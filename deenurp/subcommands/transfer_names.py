"""
Transfer names from some sequences to best-matching sequences in a reference package.
"""

import argparse
import csv
import functools
import logging

from Bio import SeqIO
from taxtastic import refpkg, taxtable

from .. import util, uclust

def build_parser(p):
    p.add_argument('refpkg', help="""Reference package""",
        type=functools.partial(refpkg.Refpkg, create=False))
    p.add_argument('fasta_file', help="""Sequence file to augment reference
        package with""", metavar='fasta_file')
    p.add_argument('seq_info', help="""Sequence info for fasta_file""",
        type=argparse.FileType('r'))
    p.add_argument('taxtable', help="""Taxtable for sequences in fasta_file""",
        type=argparse.FileType('r'))
    p.add_argument('--conflict-action', default='warn', choices=('warn', 'replace'),
        help="""Action to take when an already-named sequence matches a
        sequence in fasta_file. 'warn': show a message, do not change
        taxonomy; 'replace': replace the tax_id for the matching sequence
        [default:%(default)s]""")

    p.add_argument('--log', type=argparse.FileType('w'), help="""Log renaming
        to file""")

    p.add_argument('-i', '--percent-id', type=float, default=0.99, help="""Minimum
        percent ID to transfer taxonomy [default: %(default).2f]""")

def add_to_taxonomy(taxonomy, tax_node):
    """Add tax_node to taxonomy, including any missing"""
    to_add = []
    node = tax_node
    while True:
        if node.tax_id in taxonomy.index:
            break
        if node.is_root:
            raise ValueError("No overlap between taxonomies")
        to_add.append((node.parent.tax_id, node))
        node = node.parent

    # Add in reverse order
    to_add = to_add[::-1]
    logging.info("Adding %s", ', '.join('{0}[{1}]'.format(node.name, node.tax_id) for _, node in to_add))
    for parent_id, node in to_add:
        taxonomy.get_node(parent_id).add_child(taxtable.TaxNode(rank=node.rank,
            name=node.name, tax_id=node.tax_id,
            sequence_ids=node.sequence_ids))

        # If not present, add rank
        if not node.rank in taxonomy.ranks:
            # Get parent rank
            par_rank_idx = taxonomy.ranks.index(node.parent.rank)
            if par_rank_idx == -1:
                raise ValueError("Rank {0} not found".format(node.parent.rank))

            taxonomy.ranks.insert(par_rank_idx+1, node.rank)
            logging.info("Added rank %s", node.rank)

def ungap(seqrecord):
    seqrecord.seq = seqrecord.seq.ungap('-')
    return seqrecord

def action(args):
    log_writer = None
    if args.log:
        log_writer = csv.DictWriter(args.log, ['seqname', 'orig_tax_id',
            'renamed_tax_id', 'renamed_tax_name', 'best_hit', 'pct_id',
            'applied'], quoting=csv.QUOTE_NONNUMERIC, lineterminator='\n')
        log_writer.writeheader()

    # Load all tax_ids
    with args.taxtable as fp:
        new_tax = taxtable.read(fp)

    with args.seq_info as fp:
        new_seq_info = {row['seqname']: row for row in csv.DictReader(fp)}

    with args.refpkg.open_resource('aln_fasta') as fp:
        ref_sequences = [ungap(i) for i in SeqIO.parse(fp, 'fasta')]
    with args.refpkg.open_resource('seq_info') as fp:
        ref_seq_info_reader = csv.DictReader(fp)
        ref_seq_info = {row['seqname']: row for row in ref_seq_info_reader}
    with args.refpkg.open_resource('taxonomy') as fp:
        ref_taxonomy = taxtable.read(fp)

    search = functools.partial(uclust.search,
            pct_id=args.percent_id,
            trunclabels=True,
            search_pct_id=0.9, quiet=True)

    # Search the sequences from the reference package against the input sequences
    with util.as_fasta(ref_sequences) as ref_fasta_path, util.ntf(prefix='uclust') as tf:
        search(args.fasta_file, ref_fasta_path, tf.name)
        input_records = uclust.parse_uclust_out(i for i in tf if i.startswith('H'))

        # Also search sequences from the reference package against themselves
        # TODO: decide if we want to use this
        #with util.ntf(prefix='uclust') as self_tf:
            #search(ref_fasta_path, ref_fasta_path, self_tf.name, maxaccepts=10)
            #ref_records = uclust.parse_uclust_out(i for i in self_tf if i.startswith('H'))
            ## Drop self-hits
            #ref_records = (i for i in ref_records if i.query_label != i.target_label)
            #grouped = itertools.groupby(ref_records, operator.attrgetter('query_label'))
            #best_hit_id = dict((g, max(i.pct_id for i in v)) for g, v in grouped)

        for record in input_records:

            ref_si = ref_seq_info[record.query_label]
            target_si = new_seq_info[record.target_label]
            #if record.pct_id > best_hit_id.get(record.query_label, 0.0):
            tax_id = target_si['tax_id']
            node = new_tax.get_node(tax_id)

            if log_writer:
                log_record = {'seqname': record.query_label,
                              'best_hit': record.target_label,
                              'pct_id': record.pct_id,
                              'orig_tax_id': ref_si['tax_id'],
                              'renamed_tax_id': node.tax_id,
                              'renamed_tax_name': node.name,
                              'applied': not ref_si['tax_id'] or args.conflict_action == 'replace'}
                log_writer.writerow(log_record)

            logging.info('Naming %s %s[%s,%s] based on %s (%.2f%%)', ref_si['seqname'],
                         node.name, node.tax_id, node.rank, record.target_label, record.pct_id)
            if ref_si['tax_id'] and ref_si['tax_id'] != tax_id:
                old_node = ref_taxonomy.get_node(ref_si['tax_id'])
                logging.warn('Already named: %s[%s,%s]%s',
                             old_node.name, old_node.tax_id, old_node.rank,
                             ' - replacing' if args.conflict_action == 'replace' else '')
            if not ref_si['tax_id'] or args.conflict_action == 'replace':
                ref_si['tax_id'] = target_si['tax_id']
                if tax_id not in ref_taxonomy.index:
                    add_to_taxonomy(ref_taxonomy, node)

    # Write updated taxtable, seqinfo
    with util.ntf(prefix='taxonomy-', suffix='.csv') as new_tax, \
         util.ntf(prefix='seq_info-', suffix='.csv') as new_seq_info:
        ref_taxonomy.write_taxtable(new_tax)
        new_tax.close()

        w = csv.DictWriter(new_seq_info, ref_seq_info_reader.fieldnames)
        w.writeheader()
        w.writerows(ref_seq_info.values())
        new_seq_info.close()

        args.refpkg.start_transaction()
        args.refpkg.update_file('taxonomy', new_tax.name)
        args.refpkg.update_file('seq_info', new_seq_info.name)
        args.refpkg.commit_transaction()

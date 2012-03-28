"""
Select reference sequences for inclusion
"""
import itertools
import json
import logging
import sqlite3

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from .wrap import cmalign, as_refpkg, as_fasta, tempdir, redupfile_of_seqs, \
                  voronoi, guppy_redup, pplacer

def seqrecord(name, residues, **annotations):
    sr = SeqRecord(Seq(residues), name)
    sr.annotations.update(annotations)
    return sr

def _sequence_extractor(rdp_con):
    """
    returns a function to extract sequences that match a sequence name, or are a
    type strain from one of the represented tax_ids
    """
    def _generate_in(coll):
        return '({0})'.format(','.join('?' for i in coll))

    rdp_con.row_factory = sqlite3.Row
    cursor = rdp_con.cursor()

    stmt = """SELECT * FROM sequences
WHERE name IN {0} OR (tax_id IN {1} AND is_type = 'type');"""

    def extract_seqs(sequence_ids, tax_ids):
        s = stmt.format(_generate_in(sequence_ids), _generate_in(tax_ids))
        cursor.execute(s, list(sequence_ids) + list(tax_ids))
        for i in cursor:
            sr = seqrecord(i['name'], i['residues'], species_name=i['species_name'],
                lineage=json.loads(i['lineage']), tax_id=i['tax_id'], is_type=i['is_type'])
            yield sr

    return extract_seqs

def select_sequences_for_cluster(ref_seqs, query_seqs, keep_leaves=5):
    """
    Given a set of reference sequences and query sequences
    """
    if len(ref_seqs) <= keep_leaves:
        return [i.id for i in ref_seqs]

    c = itertools.chain(ref_seqs, query_seqs)
    ref_ids = frozenset(i.id for i in ref_seqs)
    aligned = list(cmalign(c))
    with as_refpkg(i for i in aligned if i.id in ref_ids) as rp, \
             as_fasta(aligned) as fasta, \
             tempdir(prefix='jplace') as placedir, \
             redupfile_of_seqs(query_seqs) as redup_path:
        jplace = pplacer(rp.path, fasta, out_dir=placedir())
        # Redup
        guppy_redup(jplace, redup_path, placedir('redup.jplace'))
        prune_leaves = set(voronoi(placedir('redup.jplace'), keep_leaves))

    result = frozenset(i.id for i in ref_seqs) - prune_leaves
    assert len(result) == keep_leaves
    return result

def choose_references(deenurp_db, rdp_con, refs_per_cluster=5, threads=12):
    extractor = _sequence_extractor(rdp_con)
    #total_weight = deenurp_db.total_weight()

    for cluster_id, seqs, hits in deenurp_db.hits_by_cluster():
        if not hits:
            logging.info("No hits for cluster %d", cluster_id)
            continue

        hit_names = frozenset(i['best_hit_name'] for i in hits)
        tax_ids = frozenset(i['tax_id'] for i in hits)
        seqs = [seqrecord(i['name'], i['residues'], weight=i['weight']) for i in seqs]
        ref_seqs = list(extractor(hit_names, tax_ids))
        keep = select_sequences_for_cluster(ref_seqs, seqs, refs_per_cluster)
        refs = [i for i in ref_seqs if i.id in keep]

        assert len(refs) == len(keep)

        for i in refs:
            yield i

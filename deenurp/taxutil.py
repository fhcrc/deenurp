import csv

from taxtastic import taxtable

from deenurp import util

def add_cluster_taxids(taxonomy, seqinfo_iter, cluster_col='cluster', rank='otu', rank_above='species'):
    """
    For each unlabeled sequence info record in seqinfo_iter, add an OTU label
    """
    if rank not in taxonomy.ranks:
        rank_above_idx = taxonomy.ranks.index(rank_above)
        if rank_above_idx == -1:
            raise ValueError("'{0}' not in ranks [{1}]".format(rank_above, taxonomy.ranks))
        taxonomy.ranks.insert(rank_above_idx + 1, rank)

    records = (i for i in seqinfo_iter if not i['tax_id'])
    for sr in records:
        otu = sr[cluster_col]

        # TODO: temp check
        assert otu.startswith('otu_')
        sr['tax_id'] = otu
        assert otu
        if otu not in taxonomy.index:
            taxonomy.add_child(taxtable.TaxNode(tax_id=otu, rank=rank, name=otu))

def add_clusters_to_refpkg(refpkg, **kwargs):
    with refpkg.open_resource('taxonomy') as tax_fp:
        tax = taxtable.read(tax_fp)
    with refpkg.open_resource('seq_info') as sinfo_fp:
        reader = csv.DictReader(sinfo_fp)
        sinfo = list(reader)

    # Annotate
    add_cluster_taxids(tax, sinfo, **kwargs)

    with util.ntf(prefix='seq_info-', suffix='.csv') as seqinfo_tf, \
         util.ntf(prefix='taxonomy-', suffix='.csv') as tax_tf:
        w = csv.DictWriter(seqinfo_tf, reader.fieldnames)
        w.writeheader()
        w.writerows(sinfo)
        seqinfo_tf.close()

        tax.write_taxtable(tax_tf)
        tax_tf.close()

        refpkg.start_transaction()
        refpkg.update_file('seq_info', seqinfo_tf.name)
        refpkg.update_file('taxonomy', tax_tf.name)
        refpkg.commit_transaction()

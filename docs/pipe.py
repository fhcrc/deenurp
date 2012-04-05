import subprocess
import networkx

g = networkx.DiGraph()


SEQUENCE_DB = 'sequence database(s)'
QUERY = 'query sequences'
CLUSTERS = 'clusters'
HITS = 'best hits'
MERGED = 'merged clusters'
CANDIDATES = 'candidate references'
REFERENCES = 'reference sequences'

g.add_edge(QUERY, HITS)
g.add_edge(SEQUENCE_DB, HITS)
g.add_edge(QUERY, CLUSTERS)
g.add_edge(HITS, MERGED)
g.add_edge(CLUSTERS, MERGED)
g.add_edge(MERGED, CANDIDATES)
g.add_edge(CANDIDATES, REFERENCES, label='voronoi')
with open('pipe.dot', 'w') as fp:
    networkx.write_dot(g, fp)

with open('pipe.svg', 'w') as fp:
    subprocess.check_call('dot -T svg pipe.dot', shell=True, stdout=fp)

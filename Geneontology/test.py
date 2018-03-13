import json
import scipy as sp
import scipy.stats
import numpy as np
from collections import defaultdict
from filenames import COMPLEX_BETWEENNESS, COMPLEX_CLOSENESS, COMPLEX_DEGREE, \
    GENE_ID_CONVERSION


name_map = {
    COMPLEX_BETWEENNESS: "betweenness_centrality",
    COMPLEX_CLOSENESS: "closeness_centrality",
    COMPLEX_DEGREE: "degree_centrality"
}

def parse_json(measure):
    measure_map = defaultdict(list)
    data = json.loads(open(measure).read())
    for dp in data:
        for node in dp["nodes"]:
            measure_map[node].append(dp[name_map[measure]])

    measure_map = dict(list(map(lambda k: (k[0], max(k[1])), list(measure_map.items()))))
    return measure_map

print()

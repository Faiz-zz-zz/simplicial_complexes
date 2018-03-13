import json
import scipy as sp
import scipy.stats
import numpy as np
from collections import defaultdict
from filenames import COMPLEX_BETWEENNESS, COMPLEX_CLOSENESS, COMPLEX_DEGREE, \
    GENE_ID_CONVERSION
from go_script import generate_matrix


name_map = {
    COMPLEX_BETWEENNESS: "betweenness_centrality",
    COMPLEX_CLOSENESS: "closeness_centrality",
    COMPLEX_DEGREE: "degree_centrality"
}


def get_actual_map():
    data = open("convertor.txt").read().split('\n')
    m = {}
    for each in data[1:]:
        try:
            gene_id, gene_symbol, _ = each.split('\t')
            m[gene_symbol] = gene_id
        except:
            pass
    return m


def parse_json(measure):
    measure_map = defaultdict(list)
    data = json.loads(open(measure).read())
    for dp in data:
        for node in dp["nodes"]:
            measure_map[node].append(dp[name_map[measure]])

    measure_map = dict(list(map(lambda k: (k[0], max(k[1])), list(measure_map.items()))))
    return measure_map


def generate_measure_matrix():
    gene_id_converter = get_actual_map()
    matrix, gene_list = generate_matrix()
    gene_measure = parse_json(COMPLEX_CLOSENESS)
    gene_measure_list = []
    for gene in gene_list:
        try:
            gene_measure_list.append(gene_measure[gene_id_converter[str(gene)]])
        except:
            gene_measure_list.append(0)  # me sorry
    return matrix, gene_measure_list


def calculate_regression():
    matrix, measures = generate_measure_matrix()
    go_measures = matrix[0]
    # print(len(go_measures), len(measures), "SUP SUP SUP")
    slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(np.asarray(measures), go_measures)
    print("Coefficient of determination: {}".format(r_value ** 2))
# calculate_regression()

def what():
    gene_list_ = parse_json(COMPLEX_CLOSENESS).keys()
    print(len(gene_list_))
    gene_list = set()
    gene_id_converter = get_actual_map()
    gene_id_converter = dict(zip(gene_id_converter.values(), gene_id_converter.keys()))
    for a in gene_list_:
        if a in gene_id_converter:
            gene_list.add(gene_id_converter[a])
    convertor_list = set((get_actual_map().keys()))
    print(set(gene_list) & convertor_list, len(set(gene_list) & convertor_list))

what()

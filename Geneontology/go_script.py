import csv
import pandas as pd
import numpy
import operator
from filenames import COMPLEX_BETWEENNESS, COMPLEX_CLOSENESS, COMPLEX_DEGREE, \
    GENE_ID_CONVERSION
from collections import defaultdict
import json


measures = ["betweenness", "closeness", "degree"]

name_map = {
    COMPLEX_BETWEENNESS: "betweenness_centrality",
    COMPLEX_CLOSENESS: "closeness_centrality",
    COMPLEX_DEGREE: "degree_centrality"
}

rev_map = {
    "betweenness": COMPLEX_BETWEENNESS,
    "closeness": COMPLEX_CLOSENESS,
    "degree": COMPLEX_DEGREE
}

def convert_to_csv():
    with open('sgd_merged.txt') as fin, open('csv_format.csv', 'w') as fout:
        for line in fin:
            fout.write(line.replace('\t', ','))

def delete_column():
    with open("csv_format.csv","rb") as source:
        rdr= csv.reader( source )
        with open("result.csv","wb") as result:
            wtr= csv.writer( result )
            for r in rdr:
                wtr.writerow( (r[0], r[1], r[3], r[4]) )



def parse_json(measure):
    measure = rev_map[measure]
    measure_map = defaultdict(list)
    data = json.loads(open(measure).read())
    for dp in data:
        for node in dp["nodes"]:
            measure_map[node].append(dp[name_map[measure]])

    measure_map = dict(list(map(lambda k: (k[0], max(k[1])), list(measure_map.items()))))
    data_ppi = json.loads(open(measure.replace(".json", "_PPI.json")).read())
    for data in measure_map:
        data_ppi[data] = measure_map[data]
    return data_ppi

def get_all_genes():
    return list(parse_json("betweenness").keys())

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

def generate_matrix():
    annotation_map = json.loads(open("annotation_map.json").read())
    gene_list = annotation_map.keys()
    go_ids = set()
    for gos in annotation_map.values():
        go_ids |= set(gos)
    go_ids = list(go_ids)
    matrix = []
    for go_id in go_ids:
        matrix.append([])
        for gene in gene_list:
            matrix[-1].append(1 if go_id in annotation_map[gene] else 0)
    return matrix, gene_list, go_ids
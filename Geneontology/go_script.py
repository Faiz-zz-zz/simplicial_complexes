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

def mapping():
    df = pd.read_csv('csv_format.csv', sep='\t')
    #print(df.columns)
    no_of_rows = df.shape[0]
    gene_df = df['DB_Object_Symbol']
    annotation_df = df['GO ID']
    exp_valid_df = df['Evidence']
    annotation_map = {}
    for i in range(0, no_of_rows):
        if(exp_valid_df[i] == 'EXP' or
        exp_valid_df[i] == 'IDA' or
        exp_valid_df[i] == 'IPI' or
        exp_valid_df[i] == 'IMP' or
        exp_valid_df[i] == 'IGI' or
        exp_valid_df[i] == 'IEP'):
            if(gene_df[i] in annotation_map):
                array = annotation_map[gene_df[i]]
                if(annotation_df[i] not in array):
                    array.append(annotation_df[i])
                    annotation_map[gene_df[i]] = array
            else:
                annotation_map[gene_df[i]] = [annotation_df[i]]
    # print(gene_df[10], map[gene_df[10]])
    
    return annotation_map

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
    df = pd.read_csv('csv_format.csv', sep='\t')
    no_of_genes = df.shape[0]
    annotation_df = df['GO ID']

    all_genes_id = get_all_genes()
    all_genes_symbol = []
    m = get_actual_map()
    m = dict(zip(m.values(), m.keys()))
    for gene in all_genes_id:
        if gene in m:
            all_genes_symbol.append(m[gene])
    all_genes = all_genes_symbol
    annotation_map = mapping()
    existing_genes_ann_map = {}

    for gene in annotation_map:
        if gene in all_genes:
            existing_genes_ann_map[gene] = annotation_map[gene]
    

    go_ids_occ = {}
    for gene, go_ids in existing_genes_ann_map.items():
        for go_id in go_ids:
            if go_id in go_ids_occ:
                go_ids_occ[go_id] += 1
            else:
                go_ids_occ[go_id] = 1
    
    most_occuring_go_ids = list(map(lambda k: k[0], list(sorted(list(zip(go_ids_occ.keys(), go_ids_occ.values())), key=lambda k: k[1]))))[-10:]
    matrix = []
    for go_id in most_occuring_go_ids:
        matrix.append([])
        for gene in existing_genes_ann_map:
            matrix[-1].append(1 if go_id in existing_genes_ann_map[gene] else 0)

    return matrix, annotation_map.keys(), most_occuring_go_ids

    
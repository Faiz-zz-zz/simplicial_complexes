import csv
import pandas as pd
import numpy
import operator

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

def generate_matrix():
    df = pd.read_csv('csv_format.csv', sep='\t')
    no_of_genes = df.shape[0]
    annotation_df = df['GO ID']

    go_map = {}
    for i in range(0, no_of_genes):
        if(annotation_df[i] in go_map):
            counter = go_map[annotation_df[i]]
            counter+=1
            go_map[annotation_df[i]] = counter
        else:
            go_map[annotation_df[i]] = 1
    sorted_x = sorted(go_map.items(), key=operator.itemgetter(1))
    sorted_x = sorted_x[::-1][:100]
    #go_id_list is the list of our popular GO:IDs
    go_id_list = [x[0] for x in sorted_x]
    no_go_ids = len(sorted_x)
    annotation_map = mapping()
    annotation_map_keys = list(annotation_map.keys())

    w, h = len(annotation_map_keys), no_go_ids
    # fill it up with 0's
    # matrix = [[0 for x in range(w)] for y in range(h)]
    matrix = numpy.zeros((h, w))

    for i in range(0, h):
        for j in range(0, w):
            # get the list of annotations from each gene
            values = annotation_map[annotation_map_keys[j]]
            # if the popular GO:ID matches any of the annotations associated with a gene
            if(go_id_list[i] in values or go_id_list[i] == values):
                matrix[i][j] = 1

    return matrix, annotation_map_keys, go_id_list

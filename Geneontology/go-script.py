import csv
import pandas as pd
import numpy

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
    map = {}
    for i in range(0, no_of_rows):
        if(exp_valid_df[i] == 'EXP' or
        exp_valid_df[i] == 'IDA' or
        exp_valid_df[i] == 'IPI' or
        exp_valid_df[i] == 'IMP' or
        exp_valid_df[i] == 'IGI' or
        exp_valid_df[i] == 'IEP'):
            if(gene_df[i] in map):
                array = map[gene_df[i]]
                if(annotation_df[i] not in array):
                    array.append(annotation_df[i])
                    map[gene_df[i]] = array
            else:
                map[gene_df[i]] = [annotation_df[i]]
    # print(gene_df[10], map[gene_df[10]])
    return map

def generate_matrix():
    df = pd.read_csv('csv_format.csv', sep='\t')
    no_of_genes = df.shape[0]
    annotation_df = df['GO ID']
    myset = set()
    for i in range(0, no_of_genes):
        myset.add(annotation_df[i])
    no_of_annotations = len(myset)
    # matrix = numpy.zeros((no_of_genes, no_of_annotations))    # can't use numpy cos array too big :( )
    hashmap = mapping()
    w, h = no_of_annotations, no_of_genes;
    matrix = [[0 for x in range(w)] for y in range(h)]
    for i in range(0, h):
        for j in range(0, w):
            if(annotation_df[j] in hashmap[i]):
                matrix[i][j] = 1
    print(matrix[0][0])

# convert_to_csv()
# mapping()
generate_matrix()

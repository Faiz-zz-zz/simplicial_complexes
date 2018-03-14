import json
# import scipy as sp
# import scipy.stats
import numpy as np
from sklearn import linear_model
from collections import defaultdict
from filenames import COMPLEX_BETWEENNESS, COMPLEX_CLOSENESS, COMPLEX_DEGREE, \
    GENE_ID_CONVERSION
from go_script import generate_matrix
from itertools import chain
from itertools import combinations


measures = ["betweenness", "closeness", "degree"]

name_map = {
    COMPLEX_BETWEENNESS: "betweenness_centrality",
    COMPLEX_CLOSENESS: "closeness_centrality",
    COMPLEX_DEGREE: "degree_centrality"
}


def powerset(iterable):
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(1, len(s)+1))


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
    data_ppi = json.loads(open(measure.replace(".json", "_PPI.json")).read())
    for data in measure_map:
        data_ppi[data] = measure_map[data]
    return data_ppi


def generate_measure_matrix():
    measure_map = {
        "betweenness": COMPLEX_BETWEENNESS,
        "closeness": COMPLEX_CLOSENESS,
        "degree": COMPLEX_DEGREE
    }
    gene_id_converter = get_actual_map()
    matrix, gene_list, go_ids = generate_matrix()

    measure_matrix = []
    for measure in measure_map:
        gene_measure = parse_json(measure_map[measure])
        gene_measure_list = []
        for gene in gene_list:
            try:
                gene_measure_list.append(gene_measure[gene_id_converter[str(gene)]])
            except:
                gene_measure_list.append(0)  # me sorry
        measure_matrix.append(gene_measure_list)
    return matrix, measure_matrix, go_ids, gene_list



def calculate_regression():
    for measure in measures:
        matrix, measures, _, _ = generate_measure_matrix(measure)
        go_measures = matrix[0]
        # print(len(go_measures), len(measures), "SUP SUP SUP")
        slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(np.asarray(measures_list), go_measures)
        print("Coefficient of determination for {}: {}".format(measure, r_value ** 2))


def calculate_multifit():
    measures = ["betweenness", "closeness", "degree"]
    measure_combs = powerset(measures)
    matrix, measure_list, _ = generate_measure_matrix()
    measure_dict = {}
    for i, measure in enumerate(measures):
        measure_dict[measure] = measure_list[i]

    for measure_comb in measure_combs:
        print("Calculating coeff for {}:".format(" and ".join(measure_comb)))
        X = []
        for measure in measure_comb:
            X.append(measure_dict[measure])
        X = np.asarray(X).T
        Y = matrix[0]
        clf = linear_model.LinearRegression()
        clf.fit(X, Y)
        print("Coeff is: {}".format(clf.score(X, Y)))

def get_model(train_data):
    measure_combs = powerset([0, 1, 2])
    Y, *centralities = train_data
    models = []
    for input_nums in measure_combs:
        X = []
        for i in input_nums:
            X.append(centralities[i])
        X = np.asarray(X).T
        clf = linear_model.LinearRegression()
        clf.fit(X, Y)
        models.append(clf)
    return models

# def get_train_test_matrix(gene_test, gene_list, go_matrix, measure_matrix):
#     index_list = []
#     go_train_matrix = [i for i in go_matrix]
#     measure_train_matrix = [i for i in measure_matrix]

#     go_test_matrix = []
#     measure_test_matrix = []

#     for elem in gene_test:
#         index_list.append(gene_list.index(elem))

#     for idx in index_list:
#         go_test_matrix.append(go_matrix.index(idx))
#         go_train_matrix.remove(idx)

#         measure_test_matrix.append(measure_matrix.index(idx))
#         measure_train_matrix.remove(idx)

#     return go_train_matrix, measure_train_matrix, go_test_matrix, measure_test_matrix


def zip_measure_with_go(measure_matrix, go_matrix):

    list_measure_go = []
    for go_id in go_matrix:
        list_measure_go.append(list(zip(go_id, measure_matrix)))

    return list_measure_go

def remove_all_test_genes(index_list, go_id_to_measure_map):
    for idx in index_list:
        go_id_to_measure_map.remove(idx)

def remove_test_genes_data(measure_matrix, go_matrix, gene_list, test_data_genes):

    go_id_to_measure_map = []
    # test_gene_data = []
    index_list = []
    gene_data = []

    for elem in gene_test:
        index_list.append(gene_list.index(elem))

    for go_id in go_matrix:
        go_id_to_measure_map = list(zip(go_id, measure_matrix[0], measure_matrix[1], measure_matrix[2]))
        remove_all_test_genes(index_list, go_id_to_measure_map)
        gene_data.append(go_id_to_measure_map)

    return gene_data


def build_model_input(gene_data):
    model_input = []
    for go_id in gene_data:
        ret = []
        for i in range(4):
            ret.append(list(map(lambda k: k[i], go_id)))
        model_input.append(ret)
    return model_input

def predict():
    go_matrix, measure_matrix, go_ids, gene_list = generate_measure_matrix()
    train_data_genes, test_data_genes = train_test_split(gene_list, test_size=0.2)
    genes_data = remove_test_genes_data(measure_matrix, go_matrix, gene_list, test_data_genes)

    model_input = build_model_input(genes_data)

    models_for_each_go_id = []
    for data in model_input:
        models = get_model(data)
        models_for_each_go_id.append(models)



    # list_measure_go = zip_measure_with_go(measure_matrix, go_matrix)






#predict(["m", "a", "d", "g", "h", "b", "n"])
#calculate_regression()

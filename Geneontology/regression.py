import json
# import scipy as sp
# import scipy.stats
import numpy as np
from sklearn import linear_model
from sklearn.linear_model import LogisticRegression, LinearRegression
from sklearn.model_selection import train_test_split
from collections import defaultdict
from filenames import COMPLEX_BETWEENNESS, COMPLEX_CLOSENESS, COMPLEX_DEGREE, \
    GENE_ID_CONVERSION
from go_script import generate_matrix
from itertools import chain
from itertools import combinations
from sklearn.ensemble import RandomForestRegressor, RandomForestClassifier
from sklearn.cluster import AgglomerativeClustering


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
    matrix, gene_list, go_ids = generate_matrix()

    measure_matrix = []
    for measure in measure_map:
        gene_measure = parse_json(measure_map[measure])
        gene_measure_list = []
        for gene in gene_list:
            try:
                gene_measure_list.append(gene_measure[gene])
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

global_method = RandomForestRegressor
global_params = {}

def get_model(train_data):
    measure_combs = powerset([0, 1, 2])
    Y, *centralities = train_data
    models = []
    for input_nums in measure_combs:
        X = []
        for i in input_nums:
            X.append(centralities[i])
        X = np.asarray(X).T
        clf = global_method(**global_params)
        clf.fit(X, Y)
        models.append((clf, input_nums))
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
    test_data = []
    for idx in index_list:
        test_data.append(go_id_to_measure_map[idx])
        go_id_to_measure_map.remove(idx)
    return test_data


# def remove_test_genes_data(measure_matrix, go_matrix, gene_list, test_data_genes):

#     go_id_to_measure_map = []
#     test_gene_data = []
#     index_list = []
#     gene_data = []

#     for elem in gene_test:
#         index_list.append(gene_list.index(elem))

#     for go_id in go_matrix:
#         go_id_to_measure_map = list(zip(go_id, measure_matrix[0], measure_matrix[1], measure_matrix[2]))
#         test_data = remove_all_test_genes(index_list, go_id_to_measure_map)
#         test_gene_data.append(test_data)
#         gene_data.append(go_id_to_measure_map)

#     return gene_data, test_data


def build_model_input(gene_data):
    data = []
    for i in range(4):
        data.append(list(map(lambda k: k[i], gene_data)))
    return data

def build_prediction_data(idx, test_data):
    go_id_data = test_data.index(idx)
    prediction_data = build_model_input(go_id_data)
    return prediction_data


def predict_lin_regression():
    go_matrix, measure_matrix, go_ids, gene_list = generate_measure_matrix()
    predictions = {}
    for ind, go_id in enumerate(go_matrix):
        train_data, test_data = train_test_split(list(zip(go_id, measure_matrix[0], measure_matrix[1], measure_matrix[2])), test_size=0.2)
        input_data = build_model_input(train_data)
        pred_data = build_model_input(test_data)
        models = get_model(input_data)
        go_id_models_pred = {}
        for model in models:
            model, centralities = model
            model_cent_pred = []
            pred_input = []
            for c in centralities:
                pred_input.append(pred_data[c + 1])
            pred_input = np.asarray(pred_input).T
            pred = model.predict(pred_input)
            go_id_models_pred[tuple(centralities).__str__()] = {"prediction": np.ndarray.tolist(pred), "rms": model.score(pred_input, pred_data[0])}
        predictions[go_ids[ind]] = go_id_models_pred

    return predictions

def predict_svm():
    go_matrix, measure_matrix, go_ids, gene_list = generate_measure_matrix()
    predictions = {}
    for ind, go_id in enumerate(go_matrix):
        train_data, test_data = train_test_split(list(zip(go_id, measure_matrix[0], measure_matrix[1], measure_matrix[2])), test_size=0.2)
        input_data = build_model_input(train_data)
        pred_data = build_model_input(test_data)
        if not sum(input_data[0]): continue
        models = get_model(input_data)
        go_id_models_pred = {}
        for model in models:
            model, centralities = model
            model_cent_pred = []
            pred_input = []
            for c in centralities:
                pred_input.append(pred_data[c + 1])
            pred_input = np.asarray(pred_input).T
            pred = model.predict(pred_input)
            go_id_models_pred[tuple(centralities).__str__()] = {"prediction": pred, "actual": pred_data[0]}
            # pred = list(map(lambda k: 0 if k > 0.0 else 1, pred))
        predictions[go_ids[ind]] = go_id_models_pred

    return predictions


# This is a helper function to make verification easier for me
def zip_all_train_data(train_map):
    new_map = {}
    for key, value in train_map.items():
        new_map[key] = list(zip(*value))
    return new_map

def cluster_go_ids():
    matrix, _,  _ = generate_matrix()
    trans_matrix = [list(x) for x in zip(*matrix)]
    clustering = AgglomerativeClustering(linkage='ward', n_clusters=200)
    train_data, test_data = train_test_split(trans_matrix, test_size=0.2)
    #print("this is train data {}".format(train_data))

    print("About to run clustering...")
    clustering.fit(train_data)
    #print("these is the size {}".format(len(clustering.labels_)))
    training_map = {}
    test_map = {}
    for i, elem in enumerate(clustering.labels_):
        if elem in training_map:
            training_map[elem].append(train_data[i])
        else:
            training_map[elem] = [train_data[i]]

    # for key, val in training_map.items():
    #     print("this is the key {} and value {}".format(key, val))

    print("About to predict")
    ret = clustering.fit_predict(test_data, y=None)
    print("Prediction done")

    for i, elem in enumerate(ret):
        if elem in test_map:
            test_map[elem].append(test_data[i])
        else:
            test_map[elem] = [test_data[i]]
    

    new_map = zip_all_train_data(training_map)
    for key, val in new_map.items():
        print("this is the key {} this is the value {}".format(key, val))
        if key in test_map:
            print("this is the predicted in the same cluster {}".format(test_data[key]))
        else:
            continue

cluster_go_ids()
# methods = {
#     "random_forest_classifier": [RandomForestClassifier, {"n_estimators": 2}],
#     "ranodm_forest_regressor": [RandomForestRegressor, {"n_estimators": 2}],
#     "linear_regression": [LinearRegression, {}],
#     "logistic_regression": [LogisticRegression, {}]
# }

# total_pred = {}

# for name, method in methods.items():
#     global_method, global_params = method
#     total_pred[name] = predict_lin_regression()

# import json

# with open("all_pred_ppi.json", "w") as out:
#     json.dump(total_pred, out)



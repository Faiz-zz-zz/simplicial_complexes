import json
# import scipy as sp
# import scipy.stats
import numpy as np
import heapq 
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

def append_centrality(matrix, measure_matrix):
    degree_dist = measure_matrix[2]
    betweennes = measure_matrix[1]
    closeness = measure_matrix[0]
    for index, item in enumerate(matrix):
        item.append(degree_dist[index])
        item.append(closeness[index])
        item.append(betweennes[index])
    return matrix

def count_annotated_genes(go_id_vector):
    count = 0
    for annotation in go_id_vector:
        if annotation == 1:
            count += 1
    return count

# create the representation map to analyse over represented GO Ids
def create_rep_goID_map(train_map, go_ids_list):
    rep_map = {}
    for cluster in train_map:
        rep_map[cluster] = {}
        go_ids_vectors = train_map[cluster]
        for i, go_id_vector in enumerate(go_ids_vectors):
            num_of_genes_annotated = count_annotated_genes(go_id_vector)
            if num_of_genes_annotated:
                rep_map[cluster][go_ids_list[i]] = go_id_vector
        rep_map[cluster]["genes_list"] = train_map[cluster][-1]
    return rep_map

def find_most_dominant_clusters(clustering_model):
    freq_map = {}
    heap = []
    for cluster in clustering_model.labels_:
        if cluster in freq_map:
            freq_map[cluster] +=1
        else:
            freq_map[cluster] = 1
    
    for k, v in freq_map.items():
        heapq.heappush(heap, (v, k))

    return heapq.nlargest(10, heap)

def append_gene_names(train_data_matrix, gene_list):
    for i, data in enumerate(train_data_matrix):
        data.append(gene_list[i])
    return train_data_matrix


def build_david_map(final_map, dominant_clusters):
    david_map = {}
    for _, (_,cluster_id) in enumerate(dominant_clusters):
        david_map[cluster_id] = final_map[cluster_id]
    return david_map

def get_clusters(dominant_clusters):
    clusters = []
    for _, (_, cluster_id) in enumerate(dominant_clusters):  
        clusters.append(cluster_id)
    return clusters


# returns the top most presented go ids among the genes in a cluster
def get_annotation_map(genes_in_cluster):
    genes_annotation_data = json.loads(open("annotation_map.json").read())
    go_id_map = {}
    heap = []
    for gene in genes_in_cluster:
        annotations = genes_annotation_data[gene]
        for annot in annotations:
            if annot in go_id_map:
                go_id_map[annot] += 1
            else:
                go_id_map[annot] = 1
    
    for go_id, freq in go_id_map.items():
        heapq.heappush(heap, (freq, go_id))

    top_go_ids_tuples = heapq.nlargest(10, heap)
    top_go_ids_list = []
    for _, (_, go_id) in enumerate(top_go_ids_tuples):
        top_go_ids_list.append(go_id)
    
    return top_go_ids_list


# N is the size of the cluster (only annotated genes from the cluster are taken into account),
# X is the number of genes in the cluster that are annotated with the GO term in question,
# M is the number of all genes in the network that are annotated with any GO term, 
# K is the number of genes in the network that are annotated with the GO term in question. 
# A cluster is significantly enriched in a given GO term if the corresponding p-value is smaller than or equal to 0.05. 

def get_N_and_X(cluster, go_term, genes_annotation_data):
    genes_list = cluster[0]
    N = 0 
    X = 0
    for gene in genes_list:
        if gene in genes_annotation_data:
            N += 1
            annot_list = genes_annotation_data[gene]
            if go_term in annot_list:
                X += 1
    return N, X

def get_M_and_K(genes_list, go_term, genes_annotation_data):
    M = 0
    K = 0
    for gene in genes_list:
        if gene in genes_annotation_data:
            M +=1
            annot_list = genes_annotation_data[gene]
            if go_term in annot_list:
                K += 1
    return M, K

# Each cluster is a list of genes and go ids in it
# Genes list is the list of genes in the whole network
def calc_p_value(cluster, go_term, genes_list):
    genes_annotation_data = json.loads(open("annotation_map.json").read())
    N, X = get_N_and_X(cluster, go_term, genes_annotation_data)
    M, K = get_M_and_K(genes_list, go_term, genes_annotation_data)
    p_value =  #calculation using the variables above
    # The formula is here https://www.nature.com/articles/srep35098#supplementary-information
    return p_value
    

#  We are only considering top 10 clusters for now
# You can change the value

def calculate_p_value_for_all_clusters():
    final_data, gene_list, go_ids = cluster_based_on_centrality()
    all_clusters_p_vals = {}
    for cluster_id, cluster_data in final_data.items():
        p_val_map = {}
        for go_id in go_ids:
            p_value = calc_p_value(cluster_data, go_id, gene_list)
            p_val_map[go_id] = p_value
        all_clusters_p_vals[cluster_id] = p_val_map
    
    # This is how it should look like
    # {"cluster_id": {"GO_id":"p_value"}}
    return all_clusters_p_vals        


def cluster_based_on_centrality():
    _, measure_matrix, go_ids, gene_list  = generate_measure_matrix()
    train_data = list(zip(*measure_matrix))
    clustering = AgglomerativeClustering(linkage='ward', n_clusters=400)
    print("About to run clustering...")
    clustering.fit(train_data)
    dominant_clusters = find_most_dominant_clusters(clustering)
    clusters_list = get_clusters(dominant_clusters)

    clusters_map = {}
    gene_list = list(gene_list)
    for gene_index, cluster in enumerate(clustering.labels_):
        if cluster in clusters_list:
            if cluster in clusters_map:
                clusters_map[cluster].append(gene_list[gene_index])
            else:
                clusters_map[cluster] = [gene_list[gene_index]]
       
    final_data = {}

    for cluster_id, genes_list_in_cluster in clusters_map.items():
        genes_and_top_go_ids = []
        genes_and_top_go_ids.append(genes_list_in_cluster)
        top_go_ids_in_cluster = get_annotation_map(genes_list_in_cluster)
        genes_and_top_go_ids.append(top_go_ids_in_cluster)
        final_data[cluster_id] = genes_and_top_go_ids


    for k, v in final_data.items():
        print("This is the cluster id {} \n This is the genes list inside this cluster{} \n This is the top 10 GO ids in this cluster{} \n".format(k, v[0], v[1]))

    # the format of final_data:
    # final_data { "cluster_id": [[genes_list],[top_10_go_ids]]}
    return final_data, list(gene_list), go_ids


# def cluster_based_on_go_id_centrality():
#     matrix, measure_matrix, go_ids_list, gene_list  = generate_measure_matrix()

#     trans_matrix = [list(x) for x in zip(*matrix)]
#     train_data_matrix = append_centrality(trans_matrix, measure_matrix)

#     clustering = AgglomerativeClustering(linkage='ward', n_clusters=450) 
#     print("About to run clustering...")
#     clustering.fit(train_data_matrix)

#     append_gene_names(train_data_matrix, list(gene_list))

#     training_map = {}
#     dominant_clusters = find_most_dominant_clusters(clustering)
#     print("thesea are the most dominant ones {}".format(dominant_clusters))

#     for i, elem in enumerate(clustering.labels_):
#         if elem in training_map:
#             training_map[elem].append(train_data_matrix[i])
#         else:
#             training_map[elem] = [train_data_matrix[i]]



    # for key, val in training_map.items():
    #     print("this is the key {} and value {}".format(key, val))

    # print("About to predict")
    # ret = clustering.fit_predict(test_data, y=None)
    # print("Prediction done")

    # for i, elem in enumerate(ret):
    #     if elem in test_map:
    #         test_map[elem].append(test_data[i])
    #     else:
    #         test_map[elem] = [test_data[i]]
    



    # new_map = zip_all_train_data(training_map)
    # final_map = create_rep_goID_map(new_map, go_ids_list)
    # david_input_data_map = build_david_map(final_map, dominant_clusters)

    # new_map2 = zip_all_train_data(test_map)

    # for key, val in david_input_data_map.items():
    #     print("this is the cluster id {} \n".format(key))
    #     for i, item in val.items():
    #         print("this is the go id {} this is the value {} \n".format(i, item))
        # if key in test_map:
        #     print("this is the predicted in the same cluster\n")
        #     v = new_map2[key]
        #     for i, item in enumerate(v):
        #         print("this is the go id {} this is the value {} \n".format(i, item))
        # else:
        #     continue

cluster_based_on_centrality()
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



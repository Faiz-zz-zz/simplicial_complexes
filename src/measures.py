
import networkx as nx
import pandas as pd
from tqdm import tqdm
import numpy as np



def find_shortest_path():
    pass


def generate_PPI_network(file_path):
    #need to redo this indexing part 

    # data_set_location = "../raw_data/"
    df = pd.read_csv(file_path)
    nodes_a = df['Entrez Gene Interactor A']
    nodes_b = df['Entrez Gene Interactor B']
    # get a list of all the nodes in the network
    nodes = list(set(pd.concat([nodes_a, nodes_b])))
    # get a list of all the edges
    edges = list(
        zip(df['Entrez Gene Interactor A'], df['Entrez Gene Interactor B'])
    )
  
    G = nx.Graph()
    G.add_nodes_from(nodes)
    G.add_edges_from(edges)
  
    return G
  
def get_measures(data_set):
    Graph = generate_PPI_network(data_set)
 
    #degreeOfNodes = get_degreeOfNodes(Graph)
    # get_nodesAndInteractions(Graph)
    # clusteringCoefficient = get_clusteringCoefficient(Graph)
 
    get_degreeCentrality(Graph)
    # get_betweennessCentrality(Graph)
    # get_closenessCentrality(Graph)
    
 
 
 
def convertDictToNativeType(Centrality):
    centralityNative=dict()
    for k,v in Centrality.items():
        centralityNative[np.asscalar(k)]=v
    return centralityNative
 
 
#The following functions return a dict of nodes with their corresponding centrality values
 
def get_degreeCentrality(NXGraph):
    #The degree centrality values are normalized by dividing by the maximum possible degree in a simple graph n-1 where n is the number of nodes in G. So numbers can dodgy
    degreeCentrality = nx.degree_centrality(NXGraph)
    sum=0
    for k, v in degreeCentrality.items():
        # print("Node: {} => Degree Centrality: {}".format(k, v))
        sum+=v
    # print("Average Degree Centrality: {}".format(sum/float(len(degreeCentrality))))
    degreeCentrality=convertDictToNativeType(degreeCentrality)
    import json
    with open('degreeCentrality.json', 'w') as outfile:
        json.dump(degreeCentrality, outfile)
  
def get_closenessCentrality(NXGraph):
    closenessCentrality = nx.closeness_centrality(NXGraph)
    sum=0
    for k, v in closenessCentrality.items():
        # print("Node: {} => Degree Centrality: {}".format(k, v))
        sum+=v
    # print("Average Degree Centrality: {}".format(sum/float(len(degreeCentrality))))
    closenessCentrality=convertDictToNativeType(closenessCentrality)
    import json
    with open('closenessCentrality.json', 'w') as outfile:
        json.dump(closenessCentrality, outfile)
    
 
def get_betweennessCentrality(NXGraph):
    #Choosing second betweenness centrality parameter as a low arbitrary value makes calculation much faster, acts like a radius for node being
    #calculated for. Larger value = Longer time to wait
    sum=0
    betweennessCentrality = nx.betweenness_centrality(NXGraph)
    for k, v in betweennessCentrality.items():
        # print("Node: {} => Betweenness Centrality: {}".format(k, v))
        sum+=v
    # print("Average Betweenness Centrality: {}".format(sum/float(len(betweennessCentrality))))
    betweennessCentrality=convertDictToNativeType(betweennessCentrality)
    import json
    with open('betweennessCentrality.json', 'w') as outfile:
        json.dump(betweennessCentrality, outfile)
 
def get_degreeOfNodes(NXGraph):
    listOfNodes = NXGraph.nodes()
    degreeOfNodes=dict()
    sum=0
    for node in listOfNodes:
        # print("Node: {} => Node Degree: {}".format(node,NXGraph.degree(node)))
        sum+=NXGraph.degree(node)
        degreeOfNodes[node]=NXGraph.degree(node)
 
    # print("Average Degree Centrality: {}".format(sum/float(len(listOfNodes))));
    return degreeOfNodes
  
#returns number of Nodes and Edges
def get_nodesAndInteractions(NXGraph):
    listOfNodes = NXGraph.nodes()
    listOfEdges = NXGraph.edges()
    print("Number of Nodes: {}".format(len(listOfNodes)))
    print("Number of Edges: {}".format(len(listOfEdges)))
 
 
def get_clusteringCoefficient(NXGraph):
    clusteringCoefficient = nx.clustering(NXGraph)
    for k, v in clusteringCoefficient.items():
        print("Node: {} => clusteringCoefficient: {}".format(k, v))
    return clusteringCoefficient
 
  
def mainFunction():
    #need better parameter handling when inputting files.
   # get_measures("C:/Users/Shivam/Desktop/BIOGRID-Homosapien.csv")
   #get_measures("C:/Users/Shivam/Desktop/BIOGRID-Saccharomyces-cerevisiae-(bakers_yeast)_UPDATED.csv")
   get_measures("C:/Users/shiva/Desktop/Data_Sets_Research/BIOGRID-Homosapien_UPDATED.csv")

  
mainFunction()




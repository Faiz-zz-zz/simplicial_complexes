import networkx as nx
import pandas as pd
from tqdm import tqdm


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
   # get_degreeCentrality(Graph)
    #get_degreeOfNodes(Graph)
    #get_betweennessCentrality(Graph)
    get_closenessCentrality(Graph)



#The following return a dict of nodes with their corresponding centrality values
def get_degreeCentrality(NXGraph):
    #The degree centrality values are normalized by dividing by the maximum possible degree in a simple graph n-1 where n is the number of nodes in G. So numbers can dodgy
    degreeCentrality = nx.degree_centrality(NXGraph)
    for k, v in degreeCentrality.items():
        print("Node: {} => Degree Centrality: {}".format(k, v))

def get_closenessCentrality(NXGraph):
    #closenessCentrality = nx.closeness_centrality(NXGraph)
    # for k, v in tqdm(closenessCentrality.items()):
    #     print("Node: {} => Closeness Centrality: {}".format(k, v))

    listOfNodes = NXGraph.nodes()
    for node in tqdm(listOfNodes):
        print("Node: {} => Closeness Centrality: {}".format(node,nx.closeness_centrality(NXGraph,node)))


def get_betweennessCentrality(NXGraph):
    betweennessCentrality = nx.betweenness_centrality(NXGraph)
    for k, v in betweennessCentrality.items():
        print("Node: {} => Betweenness Centrality: {}".format(k, v))

def get_degreeOfNodes(NXGraph):
    listOfNodes = NXGraph.nodes()
    for node in listOfNodes:
        print("Node: {} => Node Degree: {}".format(node,NXGraph.degree(node)))


def mainFunctionInnitThatCallsTheOthers():
    #write your fucking functi1on calls here
    get_measures("C:/Users/shiva/Desktop/Data_Sets_Research/BIOGRID-Homosapien_UPDATED.csv")

mainFunctionInnitThatCallsTheOthers()


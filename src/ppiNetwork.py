import pandas as pd
import community
import graph as graph
from graph import Node
from graph import Edge
from graph import Graph

nodes_dict = {}

def create_nodes(node_ids):
    nodes = []
    for node_id in node_ids:
        node = graph.Node(node_id, None)
        nodes_dict[node_id] = node
        nodes.append(node)
    return nodes


def create_edges(nodes_a_id, nodes_b_id, all_nodes):
    edges = []
    for node_a_id, node_b_id in zip(nodes_a_id, nodes_b_id):
        node_a = nodes_dict[node_a_id]
        node_b = nodes_dict[node_b_id]
        edges.append(graph.Edge(node_a, node_b))
    return edges


def create_ppi_network(file_path):
    df = pd.read_csv(file_path)
    nodes_a = df['Entrez Gene Interactor A']
    nodes_b = df['Entrez Gene Interactor B']
    node_ids = list(set(pd.concat([nodes_a, nodes_b])))

    nodes = create_nodes(node_ids)
    edges = create_edges(nodes_a, nodes_b, nodes)
    ppi_network_graph = graph.Graph(nodes, edges)

    return ppi_network_graph

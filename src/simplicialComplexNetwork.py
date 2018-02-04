import networkx as nx
import pandas as pd
from collections import defaultdict
from graph import Node
from graph import Edge
from graph import Graph
from ppiNetwork import create_ppi_network


class Simplice(Node):
    def __init__(self, gene_id, complex_name):
        Node.__init__(self, gene_id)
        self.complex_name = complex_name


class SimpliceEdge(Edge):
    def __init__(self, a, b):
        Edge.__init__(self, a, b)


class Simplex(Graph):
    def __init__(self, nodes, edges):
        Graph.__init__(self, nodes, edges)


def print_degree(g):
    nodes_list = g.nodes()
    for node in nodes_list:
        print("{}: {}".format(str(node), str(g.degree(node))))


def print_simplices(simplices):
    for simplice in simplices:
        print("{} {}".format(
            simplice.id, simplice.complex_name))


def print_grouped_complexes(grouped_simplices):
    for elem in grouped_simplices:
        for simplice in elem:
            print("{} {}".format(
                simplice.id, simplice.complex_name))
        print("\n")


# This is the info from David gene conversion csv file for yeast
def get_gene_conversion_info(file_path, species_type):
    gene_ids_data = pd.read_csv(file_path)
    from_ids = gene_ids_data['From']
    to_ids = gene_ids_data['To']
    species = gene_ids_data['Species']
    gene_names = gene_ids_data['Gene Name']

    from_to_ids = {}
    for f, t, s in zip(from_ids, to_ids, species):
        if s == species_type:
            from_to_ids[f] = t

    return from_to_ids, gene_names


# This is for parsing the data in biological complexes for yeast
def get_complexes_info(file_path):
    complexes_info = pd.read_csv(file_path)
    genes = complexes_info['Name']
    complexes_names = complexes_info['Complex']
    return genes, complexes_names


def group_by_complexes(simplices):
    complex_groups = defaultdict(list)

    for simplice in simplices:
        complex_groups[simplice.complex_name].append(simplice)

    grouped_simplices = complex_groups.values()
    return grouped_simplices


# Constructing the simplices from the complexes data
def construct_simplices(
        gene_conversion_file_path,
        complexes_file_path, species):
    from_to_ids, gene_names = get_gene_conversion_info(
                                gene_conversion_file_path,
                                species
                            )
    genes, complexes_names = get_complexes_info(complexes_file_path)

    simplices = []

    for g, c in zip(genes, complexes_names):
        if g in from_to_ids:
            simplices.append(Simplice(from_to_ids[g], c))

    grouped_simplices = group_by_complexes(simplices)
    simplice_nodes = []
    ppi_network = create_ppi_network("../raw_data/yeast.csv")
    simplice_edges = []
    for simplices in grouped_simplices:
        simplice_nodes.extend(simplices)
        for node_a in simplices:
            for node_b in simplices:
                if node_a != node_b:
                    simplice_edges.append(SimpliceEdge(node_a, node_b))
    simplex = Simplex(set(simplice_nodes), set(simplice_edges))
    edges = []
    for edge in simplex.edges:
        if Edge(edge.a, edge.b) not in ppi_network.edges:
            ppi_network.edges.add(Edge(edge.a, edge.b))

    return simplex, ppi_network


simplex, ppi_network = construct_simplices(
            '../raw_data/gene_ids.csv',
            '../raw_data/CYC2008_complex_v2.csv',
            'Saccharomyces cerevisiae S288C'
        )

print("Printing Neighbours")

for node in ppi_network.nodes:
    print("Current node is :{}".format(node.id))
    print(list(map(lambda k: k.id, ppi_network.find_neighbours(node))))
    print("\n======\n")

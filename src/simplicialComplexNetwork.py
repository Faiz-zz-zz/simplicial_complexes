import networkx as nx
import pandas as pd
from collections import defaultdict
import itertools
from graph import Node
from graph import Edge
from graph import Graph
from ppiNetwork import create_ppi_network

# {(2, 3, 4) (2, 3, 1) (3, 4, 5)}
# {(2, 3) (2, 4) (2, 5) (3, 1) (3, 5) (3, 2)}


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
        self.triangles = defaultdict(set)
        self.edge_triangle_map = defaultdict(set)
        self.generate_triangles()

    def build_edge_to_triangle_map(triangle):
        for (a, b) in itertools.product(triangle, triangle):
            if a == b: continue
            self.edge_triangle_map[Edge(a, b)].add(triangle)

    def generate_triangles(self):
        triangles = set()
        for node in self.nodes:
            neighbours = self.find_neighbours(node)
            for (a, b) in itertools.product(neighbours, neighbours):
                if a == b: continue
                if Edge(a, b) in self.edges:
                    # sorting to have it unique in the set.
                    triangle = tuple([a, b, node].sort(key=lambda k: k.id))
                    if triangle not in triangles:
                        triangles.add(triangle)
                        build_edge_to_triangle_map(triangle)
        self.triangles = triangles
        return self.triangles


    def find_lower_adjacent(self, triangle):
        edges = []
        edges.append(triangle[0] + triangle[1])
        edges.append(triangle[1] + triangle[2])
        edges.append(triangle[0] + triangle[2])

        num_of_lower_adjacent = 0
        for edge in edges:
            num_of_lower_adjacent += len(self.edge_triangle_map[edge])
        return num_of_lower_adjacent


    def is_partof_k_plus_one_simplex(self, tr_a, tr_b):
        diff_a = set(tr_a) - set(tr_b)
        diff_b = set(tr_b) - set(tr_a)
        edge_a_b = (diff_a, diff_b)
        if edge_a_b in self.edge_triangle_map:
            return True
        return False


    def find_upper_adjacent(self, triangle):
        edges = []
        edges.append(triangle[0] + triangle[1])
        edges.append(triangle[1] + triangle[2])
        edges.append(triangle[0] + triangle[2])

        num_of_upper_adjacent = 0
        for edge in edges:
            triangles_edge_is_partof = self.edge_triangle_map[edge]
            for tr in triangles_edge_is_partof:
                if is_partof_k_plus_one_simplex(triangle, tr):
                    num_of_upper_adjacent += 1
        return num_of_upper_adjacent


    def find_shortest_paths(self, triangle_a, triangle_b):
        pass

    # we need to make a map of the edges to the trianbles they are part of
    # {(2, 3) -> (3, 4, 1) (2, 3, 4)}
    def degree_distribution_centrality(self, triangle):
        lower_adjacent = find_lower_adjacent(triangle)
        upper_adjacent = find_upper_adjacent(triangle)
        degree = lower_adjacent + upper_adjacent

    # list of all the triangles
    # This is implemented according to definition 15 in the paper
    def closeness_centrality(self, triangle):
        for tr in self.triangles:
            shortest_paths_sum += find_shortest_path(triangle, tr)
        return (1/shortest_paths_sum)


    def count_triangle_in_path(self, triangle, paths):
        in_path = 0
        for path in paths:
            if triangle in paths:
                in_path += 1
        return in_path


    def betweenness_centrality(self, triangle):
        betweenness = 0
        for tr_a in self.triangles:
            for tr_b in self.triangles:
                if tr_a != tr_b:
                    shortest_paths = find_shortest_paths(tr_a, tr_b)
                    num_of_shortest_paths = len(shortest_paths)
                    triangle_in_path = count_triangle_in_path(triangle, shortest_paths)
                    betweenness += triangle_in_path/num_of_shortest_paths
        return betweenness


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

for triangle in simplex.triangles:
    print("Triangle: ({}, {}, {})".format(triangle[0].id, triangle[1].id, triangle[0].id))
    print("Betweennes Centrality: {}".format(simlex.betweenness_centrality(triangle)))
    print("Closeness Centrality: {}".format(simplex.closeness_centrality(triangle)))

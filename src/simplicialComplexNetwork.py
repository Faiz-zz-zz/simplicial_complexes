# import networkx as nx
import pandas as pd
from collections import defaultdict
import itertools
from graph import Node, Edge, Graph
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

    #def __eq__(self, other):
    #    return other.a == self.a and other.b == self.b


class Simplex(Graph):
    def __init__(self, nodes, edges):
        Graph.__init__(self, nodes, edges)
        self.triangles = defaultdict(set)
        self.edge_triangle_map = defaultdict(set)
        self.edge_map = {}
        self.generate_edge_map()
        self.generate_triangles()

    def generate_edge_map(self):
        for edge in self.edges:
            self.edge_map[(edge.a, edge.b)] = True

    def build_edge_to_triangle_map(self, triangle):
        for (a, b) in itertools.product(triangle, triangle):
            if a == b: continue
            self.edge_triangle_map[Edge(a, b)].add(triangle)

    def generate_triangles(self):
        triangles = set()
        for node in self.nodes:
            neighbours = self.find_neighbours(node)
            for (a, b) in itertools.product(neighbours, neighbours):
                if a == b: continue
                if (a, b) in self.edge_map or (b, a) in self.edge_map:
                    # sorting to have it unique in the set.

                    triangle = [a, b, node]
                    triangle.sort(key=lambda k: k.id)
                    triangle = tuple(triangle)
                    if triangle not in triangles:
                        triangles.add(triangle)
                        self.build_edge_to_triangle_map(triangle)
        self.triangles = triangles
        return self.triangles


    def find_lower_adjacent(self, triangle):
        edges = []
        edges.append(Edge(triangle[0], triangle[1]))
        edges.append(Edge(triangle[1],  triangle[2]))
        edges.append(Edge(triangle[0], triangle[2]))

        num_of_lower_adjacent = 0
        neighbours = set()
        for edge in edges:
            num_of_lower_adjacent += len(self.edge_triangle_map[edge])
            neighbours |= (self.edge_triangle_map[edge])
        return num_of_lower_adjacent, list(neighbours)


    def is_partof_k_plus_one_simplex(self, tr_a, tr_b):
        diff_a = set(tr_a) - set(tr_b)
        diff_b = set(tr_b) - set(tr_a)
        edge_a_b = (diff_a, diff_b)
        if edge_a_b in self.edge_triangle_map:
            return True
        return False


    def find_upper_adjacent(self, triangle):
        edges = []
        edges.append((triangle[0], triangle[1]))
        edges.append((triangle[1],triangle[2]))
        edges.append((triangle[0],triangle[2]))

        num_of_upper_adjacent = 0
        for edge in edges:
            triangles_edge_is_partof = self.edge_triangle_map[edge]
            for tr in triangles_edge_is_partof:
                if self.is_partof_k_plus_one_simplex(triangle, tr):
                    num_of_upper_adjacent += 1
        return num_of_upper_adjacent


    def dfs_all(self, start, end):
        """
        Iterative dfs that generates all paths from start to end (keeping track of the node ordering in our own stack)
        """
        stack = [(start, [start])]
        paths = []
        while stack:
            simplex_node, path = list(stack.pop())
            if simplex_node == end:
                paths.append(path)
            _, simplices_connected = self.find_lower_adjacent(simplex_node)
            for adj in (set(simplices_connected) - set(path)):
                for i in range(len(simplices_connected)):
                    _, new_simplices_connected = self.find_lower_adjacent(simplices_connected[i])
                    if(len(set(adj) - set(new_simplices_connected)) == 1):
                        stack.append((adj, path + [adj]))
        return paths

    def find_shortest_path(self, start, end):
        """
        Given a 3 tuple start [(a,b,c)] and 3 tuple end [(x,y,z)], and a list of 3 tuples of all 3-simplicies (all_simplicies)
        This tries to find the shortest path among all the paths returned by the dfs function
        """
        all_paths = self.dfs_all(start, end)
        shortest = 0
        for i in range(len(all_paths)):
            if len(all_paths[i]) > shortest:
                all_paths.pop(i)
            else:
                shortest = len(all_paths[i])
        return all_paths



    # we need to make a map of the edges to the trianbles they are part of
    # {(2, 3) -> (3, 4, 1) (2, 3, 4)}
    def degree_distribution_centrality(self, triangle):
        lower_adjacent, neighbours = self.find_lower_adjacent(triangle)
        upper_adjacent = self.find_upper_adjacent(triangle)
        degree = lower_adjacent + upper_adjacent
        return degree/len(self.triangles)

    # list of all the triangles
    # This is implemented according to definition 15 in the paper
    def closeness_centrality(self, triangle):
        shortest_paths_sum = 0
        # paths_file = open("paths.txt", "w")
        # paths_file.write("This is the node {} {} {}".format(traingle))
        for tr in self.triangles:
            paths = self.find_shortest_path(triangle, tr)
            # paths_file.write("These are the lengths of paths " + str(len(paths)))
            if len(paths):
                shortest_paths_sum += len(paths[0])
        paths_file.close()
        return (1/shortest_paths_sum)


    def give_me_random_list(self, n):
        from random import choice
        if n < len(self.triangles): return self.triangles
        ret = set()
        while (len(ret)) < n:
            _item = random.choice(self.triangles)
            ret.add(_item)
        return ret


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
                    shortest_paths = self.find_shortest_path(tr_a, tr_b)
                    num_of_shortest_paths = len(shortest_paths)
                    triangle_in_path = self.count_triangle_in_path(triangle, shortest_paths)
                    if num_of_shortest_paths:
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
    simplex = Simplex(simplice_nodes, simplice_edges)
    edges = []
    for edge in simplex.edges:
        if Edge(edge.a, edge.b) not in ppi_network.edges:
            ppi_network.edges.add(Edge(edge.a, edge.b))

    return simplex, ppi_network

simplex.generate_triangles()
for triangle in simplex.triangles:
    #degree_dist = simplex.degree_distribution_centrality(triangle)
    closeness = simplex.closeness_centrality(triangle)
    print (closeness)


def generate_metrics(methods, data_set):

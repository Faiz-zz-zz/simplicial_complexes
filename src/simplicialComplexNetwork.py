import networkx as nx
import pandas as pd
from collections import defaultdict
import itertools
from graph import Node, Edge, Graph
from ppiNetwork import create_ppi_network
from queue import Queue
from tqdm import tqdm
import json
from filenames import RAW_HUMAN_PPI, RAW_YEAST_PPI, RAW_HUMAN_COMPLEX, \
    RAW_HUMAN_COMPLEX_JSON, RAW_YEAST_COMPLEX, CLEANED_LOCATION, GENE_ID_CSV

# {(2, 3, 4) (2, 3, 1) (3, 4, 5)}
# {(2, 3) (2, 4) (2, 5) (3, 1) (3, 5) (3, 2)}
>>>>>>> master

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
        self.model_simplicial_complex_as_binary_network()

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
        # assumed they share an edge

        diff_a = set(tr_a) - set(tr_b)
        diff_b = set(tr_b) - set(tr_a)
        if diff_a and diff_b:
            edge_a_b = Edge(Node(diff_a.pop()), Node(diff_b.pop()))
            if edge_a_b in self.edge_triangle_map:
                return True
        return False


    def find_upper_adjacent(self, triangle):
        edges = []
        edges.append(Edge(triangle[0], triangle[1]))
        edges.append(Edge(triangle[1],triangle[2]))
        edges.append(Edge(triangle[0],triangle[2]))

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


    # def find_shortest_path(self, start, end):
    #     """
    #     Given a 3 tuple start [(a,b,c)] and 3 tuple end [(x,y,z)], and a list of 3 tuples of all 3-simplicies (all_simplicies)
    #     This tries to find the shortest path among all the paths returned by the dfs function
    #     """
    #     all_paths = self.dfs_all(start, end)
    #     shortest = 0
    #     for i in range(len(all_paths)):
    #         if len(all_paths[i]) > shortest:
    #             all_paths.pop(i)
    #         else:
    #             shortest = len(all_paths[i])
    #     return all_paths

    def create_edges_from_triangles(self, tr_list):
        for tr1 in tr_list:
            for tr2 in tr_list:
                if tr1 != tr2:
                    node_a = self.triangle_to_id[tr1]
                    node_b = self.triangle_to_id[tr2]
                    self.binary_network_edges.append((node_a, node_b))

    def map_triangles_to_ids(self):
        num_of_triangles = len(self.triangles)
        self.triangle_to_id = {}
        self.id_to_triangle = {}
        triangles_list = list(self.triangles)
        for i in range(num_of_triangles):
            self.triangle_to_id[triangles_list[i]] = i
            self.id_to_triangle[i] = triangles_list[i]

    def model_simplicial_complex_as_binary_network(self):
        self.binary_network_nodes = list(range(0, len(self.triangles)))
        self.binary_network_edges = []
        self.binary_network_model = nx.Graph()
        self.map_triangles_to_ids()
        for edge, tr_list in self.edge_triangle_map.items():
            self.create_edges_from_triangles(tr_list)

        self.binary_network_model.add_nodes_from(self.binary_network_nodes)
        self.binary_network_model.add_edges_from(self.binary_network_edges)


    def get_edges(self, triangle):
        edge_one = Edge(triangle[0], triangle[1])
        edge_two = Edge(triangle[1], triangle[2])
        edge_three = Edge(triangle[2], triangle[0])
        return edge_one, edge_two, edge_three


    def add_to_path(self, paths, triangle, last_node):
        for path in paths:
            if path[-1] == last_node:
                new_path = path
                new_path.append(triangle)
                paths.append(new_path)


    def remove_path_with_last_node(self, paths, last_node):
        for path in paths:
            if path[-1] == last_node:
                paths.remove(path)
                break


    def find_shortest_paths(self, paths, end_triangle):
        shortest_paths = []
        for path in paths:
            if path[-1] == end_triangle:
                shortest_paths.append(path)
        return shortest_paths


    def shortest_path(self, start_triangle, end_triangle):
        nodes_queue = Queue()
        visited_map = {}
        nodes_queue.put(start_triangle)
        paths = []


        while not nodes_queue.empty():
            node = nodes_queue.get()
            visited_map[node] = 1
            # print("This is the node {} {} {} ".format(node[0].id, node[1].id, node[2].id))
            if node == end_triangle:
                break

            edge_one, edge_two, edge_three = self.get_edges(node)
            triangles_edge_one = self.edge_triangle_map[edge_one]
            triangles_edge_two = self.edge_triangle_map[edge_two]
            triangles_edge_three = self.edge_triangle_map[edge_three]

            triangles_to_add = triangles_edge_one | triangles_edge_two | triangles_edge_three
            triangles_to_add = triangles_to_add - {node}

            for tr in triangles_to_add:
                print("These are the neighbours for node {} {} {} , {} {} {} ".format(node[0].id, node[1].id, node[2].id, tr[0].id, tr[1].id, tr[2].id))

            for triangle in triangles_to_add:

                if triangle not in visited_map and not self.is_partof_k_plus_one_simplex(triangle, node):
                    # print("This is the nodes to add to the queue {} {} {} ".format(triangle[0].id, triangle[1].id, triangle[2].id))
                    nodes_queue.put(triangle)

                    if not paths:
                        paths.append([node, triangle])
                    else:
                        self.add_to_path(paths, triangle, node)
                        if paths:
                            self.remove_path_with_last_node(paths, node)
                        print("This is number of paths {}".format(len(paths)))

        # print("I reached to the end")
        shortest_paths = self.find_shortest_paths(paths, end_triangle)
        return shortest_paths

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
        for tr in self.triangles:
            if tr != triangle:
                source = self.triangle_to_id[triangle]
                dest = self.triangle_to_id[tr]
                paths = []
                try:
                    for p in nx.all_shortest_paths(self.binary_network_model, source, dest):
                        paths.append(p)
                except:
                    continue
                if len(paths):
                    shortest_paths_sum += len(paths[0])
        return 1 / shortest_paths_sum if shortest_paths_sum else print("shortest path sum is zero")


    def count_triangle_in_path(self, triangle, paths):
        in_path = 0
        for path in paths:
            if triangle in path:
                in_path += 1
        return in_path


    def betweenness_centrality(self, triangle):
        betweenness = 0
        print("This is the node to find betweenness {}".format(self.triangle_to_id[triangle]))
        for tr_a in self.triangles:
            for tr_b in self.triangles:
                if tr_a != tr_b:
                    source = self.triangle_to_id[tr_a]
                    dest = self.triangle_to_id[tr_b]
                    shortest_paths = []
                    try:
                        for p in nx.all_shortest_paths(self.binary_network_model, source, dest):
                            # print("This is the path {}".format(p))
                            shortest_paths.append(p)
                    except:
                        continue
                    # shortest_paths = self.shortest_path(tr_a, tr_b)
                    num_of_shortest_paths = len(shortest_paths)
                    # print("this is the num of paths {}".format(len(num_of_shortest_paths)))
                    triangle_in_path = self.count_triangle_in_path(self.triangle_to_id[triangle], shortest_paths)
                    # print("this is the number of paths triangle exist in {}".format(triangle_in_path))
                    if num_of_shortest_paths:
                        betweenness += (triangle_in_path/num_of_shortest_paths)
                        print("This is betweenness until now {}".format(betweenness))

        return betweenness/(((len(self.triangles)-1)*(len(self.triangles)-2))/2)



    def write_to_file_degree_centralities(self):
        degree_dict = nx.degree_centrality(self.binary_network_model)
        degree_centralities = {}
        data = []
        for key, value in degree_dict.items():
            # print("this is for node {} {} {}".format(str(triangle[0].id), str(triangle[1].id), str(triangle[2].id)))
            triangle = self.id_to_triangle[key]
            print("this is for node {} {} {}".format(str(triangle[0].id), str(triangle[1].id), str(triangle[2].id)))
            data.append({"nodes":[str(triangle[0].id), str(triangle[1].id), str(triangle[2].id)], "degree_centrality": value})

        with open('degree_centrality.json', 'w') as outfile:
            #json.dumps(closeness_centralities)
            json.dump(data, outfile)

    def write_to_file_closeness_centralities(self):
        closeness_dict = nx.closeness_centrality(self.binary_network_model)
        closeness_centralities = {}
        data = []
        for key, value in closeness_dict.items():
            # print("this is for node {} {} {}".format(str(triangle[0].id), str(triangle[1].id), str(triangle[2].id)))
            triangle = self.id_to_triangle[key]
            data.append({"nodes":[str(triangle[0].id), str(triangle[1].id), str(triangle[2].id)], "closeness_centrality": value})

        with open('closeness_centrality.json', 'w') as outfile:
            #json.dumps(closeness_centralities)
            json.dump(data, outfile)


    def write_to_file_betweenness_centralities(self):
        betweenness_dict = nx.betweenness_centrality(self.binary_network_model)
        betweenness_centralities = {}
        data = []
        for key, value in betweenness_dict.items():
            triangle = self.id_to_triangle[key]
            data.append({"nodes":[str(triangle[0].id), str(triangle[1].id), str(triangle[2].id)], "betweenness_centrality": value})


        # print("length of data for betweenness is  {}".format(len(data)))
        # print("length of dic for betweenness is  {}".format(len(betweenness_dict)))
        with open('betweenness_centrality.json', 'w') as outfile:
            #json.dumps(closeness_centralities)
            json.dump(data, outfile)


    def closeness_centrality_all(self):
        for triangle in self.traingles:
            self.closeness_centrality(traingle)

    def degree_distribution_centrality_all(self):
        for triangle in self.triangles:
            self.degree_distribution_centrality(triangle)

    def closeness_centrality_all(self):
        for triangle in self.triangles:
            self.closeness_centrality(triangle)

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


# simplex, ppi_network = construct_simplices(
#             '../raw_data/gene_ids.csv',
#             '../raw_data/CYC2008_complex_v2.csv',
#             'Saccharomyces cerevisiae S288C'
#         )


def generate_metrics(methods, data_set):
    simplex, ppi_network = None, None
    if data_set == "yeast_complex":
        simplex, ppi_network = construct_simplices(
            CLEANED_LOCATION + GENE_ID_CSV,
            CLEANED_LOCATION + RAW_YEAST_COMPLEX,
            'Saccharomyces cerevisiae S288C'
        )
    else:
        # Make human
        pass

    available_methods = {
        "betweenness": simplex.betweenness_centrality_all(),
        "closeness": simplex.closeness_centrality_all(),
        "degree": simplex.degree_distribution_centrality_all()
    }

    for method in methods:
        if method in available_methods:
            available_methods[method]()

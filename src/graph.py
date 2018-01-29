from collections import defaultdict


class Node:
    """
    Base class node.
    Each node in Simplicial complex will inherit from this.
    This is so that we don't need to repeat the code and can use the same
    property of nodes in every case.
    """
    def __init__(self, id, description, node_type):
        self.id = id
        self.description = description
        self.type = node_type  # Simplicial or PPI


class Edge:
    """
    Edge base class incase we need weighed graphs.
    """
    def __init__(self, a, b, weight=1):
        self.a = a
        self.b = b
        self.weight = weight


class Graph:
    """
    Parent class of all the graphs. PPI will be based on a normal graph,
    where as Simplicial Complex will inherit from parent Graph class.
    All the basic measurements in (measures(PPI) - measures(simplicial))
    will be in this class as class methods.
    """
    def __init__(self, node, edges):
        self.nodes = node
        self.edges = edges
        # Only set if find_shortest_path was called once
        self.adj_matrix = None
        self.next_matrix = None

    def find_shortest_path(self):
        """
        Find shortest path between all the pairs of nodes.
        Input: nodes, edges
        Output: {(src, dest): dist} => assumes directed.
        """
        dist_dict = defaultdict(float('inf'))
        next_dict = defaultdict(None)

        for edge in self.edges:
            dist_dict[(edge.b, edge.a)] = edge.weight
            next_dict[edge.a] = edge.b

        for k in range(len(nodes)):
            for i in range(len(nodes)):
                for j in range(len(nodes)):
                    n_k, n_i, n_j = nodes[k], nodes[i], nodes[j]
                    inter_dist = dist_dict[(n_i, n_k)] + dist_dict[(n_k, n_j)]
                    if dist_dict[(n_i, n_j)] > inter_dist:
                        dist_dict[(n_i, n_j)] = inter_dist
                        next_dict[(n_i, n_j)] = next_dict[(n_i, n_k)]

        self.next_dict = next_dict
        self.adj_matrix = dist_dict
        return self.adj_matrix

    def dijkstra(self, src, dest):
        pass

    def find_distance(self, src, dist):
        if adj_matrix:
            mat = self.adj_matrix
            return mat[(src, dist)] if (src, dist) in mat else mat[(dist, src)]
        # apply dijkstra if we are using distance for just one pair.

    def find_neighbout(self, node):
        """
        Given a node, return all the neighbours.
        Input: node
        Output: set(nodes)
        """
        return set()

    def get_path(self, src, dest):
        """
        Given a node, return the path from src -> dest.
        Input: Node
        Output: List(nodes)
        """
        return []

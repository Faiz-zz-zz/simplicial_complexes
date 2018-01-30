from collections import defaultdict
from queue import PriorityQueue


class Node:
    """
    Base class node.
    Each node in Simplicial complex will inherit from this.
    This is so that we don't need to repeat the code and can use the same
    property of nodes in every case.
    """
    def __init__(self, id, description):
        self.id = id
        self.description = description


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
        self.nodes = set(node)
        self.edges = set(edges)
        self.distance_map = {}
        # Only set if find_shortest_path was called once
        self.adj_matrix = None
        self.next_node_dict = None
        # self.neighbour_map = defaultdict(set())
        # self.build_neighbour()
        # self.build_distance_map()

    def build_distance_map(self):
        self.distance_map = {(k.a, k.b): k.weight for k in edges}

    def build_neighbour(self):
        for edge in self.edges:
            self.neighbour_map[edge.a].add(edge.b)

    def find_shortest_path(self):
        """
        Find shortest path between all the pairs of nodes.
        Input: nodes, edges
        Output: {(src, dest): dist} => assumes directed.
        """
        dist_dict = defaultdict(float('inf'))
        next_node_dict = defaultdict(None)

        for edge in self.edges:
            dist_dict[(edge.b, edge.a)] = edge.weight
            next_node_dict[edge.a] = edge.b

        nodes = list(self.nodes)
        for k in range(len(nodes)):
            for i in range(len(nodes)):
                for j in range(len(nodes)):
                    n_k, n_i, n_j = nodes[k], nodes[i], nodes[j]
                    inter_dist = dist_dict[(n_i, n_k)] + dist_dict[(n_k, n_j)]
                    if dist_dict[(n_i, n_j)] > inter_dist:
                        dist_dict[(n_i, n_j)] = inter_dist
                        next_node_dict[(n_i, n_j)] = next_node_dict[(n_i, n_k)]

        self.next_node_dict = next_node_dict
        self.adj_matrix = dist_dict
        return self.adj_matrix

    def dijkstra(self, src, dest):
        dist = {}
        prev = {}
        dist[src] = 0

        vertex_hq = PriorityQueue()
        vertex_hq.put((0, src))

        for node in self.nodes:
            if node != src:
                dist[node] = float('inf')
                prev[node] = None

        while vertex_pq:
            curr_node = vertex_pq.get()
            if curr_node == dest:
                path = [src]
                while path[-1] != dest:
                    path.append(prev[path[-1]])
                return dist[dest], prev

            for nb in self.find_neighbours(curr_node):
                alt_dist = dist[curr_node] + distance_map[(curr_node, nb)]
                if alt_dist < dist[node]:
                    dist[nb] = alt_dist
                    prev[nb] = curr_node
                    vertex_hq.put((dist[nb], nb))

    def find_distance(self, src, dist):
        if self.adj_matrix:
            mat = self.adj_matrix
            return mat[(src, dist)] if (src, dist) in mat else mat[(dist, src)]
        # apply dijkstra if we are using distance for just one pair.
        dist, _ = self.dijkstra(src, dest)
        return dist

    def find_neighbours(self, node):
        """
        Given a node, return all the neighbours.
        Input: node
        Output: set(nodes)
        """
        return self.neighbour_map[node]

    def get_path(self, src, dest):
        """
        Given a node, return the path from src -> dest.
        Input: Node
        Output: List(nodes)
        """
        if self.next_node_dict:
            if not self.next_node_dict[(src, dest)]:
                return []
            path = [src]
            while path[-1] != dest:
                path.append(self.next_node_dict[(path[-1], dest)])
            return path
        _, path = self.dijkstra(src, dest)
        return path

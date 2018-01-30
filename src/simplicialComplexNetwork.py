import networkx as nx
import pandas as pd
import community
import math as m
from collections import defaultdict
import graph as graph
from graph import Node as Node
from graph import Edge as Edge
from graph import Graph as Graph


class Simplice(Node):
    def __init__(self, gene_id, complex_name):
        Node.__init__(self, gene_id)
        self.complex_name = complex_name


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
def construct_simplices(gene_conversion_file_path, complexes_file_path, species):
    from_to_ids, gene_names = get_gene_conversion_info(gene_conversion_file_path, species)
    genes, complexes_names = get_complexes_info(complexes_file_path)

    simplices = []

    for g, c in zip(genes, complexes_names):
        if g in from_to_ids:
            simplices.append(Simplice(from_to_ids[g], c))

    grouped_simplices = group_by_complexes(simplices)
    print_grouped_complexes(grouped_simplices)
    
    return grouped_simplices


# simplices = construct_simplices('/Users/matin/desktop/gene_ids.csv','/Users/matin/desktop/CYC2008_complex_v2.csv','Saccharomyces cerevisiae S288C')


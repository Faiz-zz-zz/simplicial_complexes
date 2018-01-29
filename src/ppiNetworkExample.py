import networkx as nx
import pandas as pd
import community
import math as m


class Simplice:
    def __init__(self, gene, complex_genes, complex_name):
        self.gene = gene
        self.complex_genes = complex_genes
        self.complex_name = complex_name


def print_degree(g):
    nodes_list = g.nodes()
    for node in nodes_list:
        print("{}: {}".format(str(node), str(g.degree(node))))


def print_simplices(simplices):
    for simplice in simplices:
        print("{} {} {}".format(
            simplice.gene, simplice.complex, simplice.complex_name))


# This is the info from David gene conversion csv file for yeast
def get_gene_conversion_info():
    gene_ids_data = pd.read_csv("gene_ids.csv")
    from_ids = gene_ids_data['From']
    to_ids = gene_ids_data['To']
    species = gene_ids_data['Species']
    gene_names = gene_ids_data['Gene Name']

    from_to_ids = {}
    for f, t, s in zip(from_ids, to_ids, species):
        if s == 'Saccharomyces cerevisiae S288C':
            from_to_ids[f] = t

    return from_to_ids, species, gene_names


# This is for parsing the data in biological complexes for yeast
def get_complexes_info():
    complexes_info = pd.read_csv("Annotated_YHTP2008_complex.csv")
    genes = complexes_info['Gene']
    other_involved_genes = complexes_info['Complete known complex']
    complexes_names = complexes_info['Related known complex']
    return genes, other_involved_genes, complexes_names


# Constructing the simplices from the complexes data
def construct_simplices():
    from_to_ids, species, gene_names = get_gene_conversion_info()
    genes, other_involved_genes, complexes_names = get_complexes_info()

    simplices = []

    for g, c, i in zip(genes, complexes_names, other_involved_genes):
        if g in from_to_ids:
            simplices.append(Simplice(from_to_ids[g], i, c))

    return simplices


def basicNetwork():

    df = pd.read_csv("BIOGRID-Homosapien.csv")
    nodesA = df['Entrez Gene Interactor A']
    nodesB = df['Entrez Gene Interactor B']
    nodes = pd.concat([nodesA, nodesB]) #get a list of all the nodes in the network
    allNodes = list(set(nodes)) #remove any duplicates, for the case if a gene is repeated
    #print(allNodes)

    edges = list(zip(df['Entrez Gene Interactor A'], df['Entrez Gene Interactor B'])) #get a list of all the edges
    #print(edges)

    ##Also possible to add weighted edges as well adding nodes and edges one-at-a-time
    ##https://networkx.github.io/documentation/stable/tutorial.html
    G = nx.Graph()
    G.add_nodes_from(allNodes)
    G.add_edges_from(edges)

    print(nx.info(G)) #Gives basic info of the graph, number of nodes, edges and the average degree
    print_degree(G)

basicNetwork()
simplices = construct_simplices()
print_simplices(simplices)

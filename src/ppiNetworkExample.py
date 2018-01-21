import networkx as nx 
import pandas as pd
import community 

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

basicNetwork()

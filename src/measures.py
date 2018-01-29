def find_shortest_path():
    pass


def generate_PPI_network(data_set):
    data_set_location = "../raw_data/"
    df = pd.read_csv(data_set_location + data_set + ".csv")
    nodes_a = df['Entrez Gene Interactor A']
    nodes_b = df['Entrez Gene Interactor B']
    # get a list of all the nodes in the network
    nodes = list(set(pd.concat([nodes_a, nodes_b])))
    # get a list of all the edges
    edges = list(
        zip(df['Entrez Gene Interactor A'], df['Entrez Gene Interactor B'])
    )

    G = nx.Graph()
    G.add_nodes_from(nodes)
    G.add_edges_from(edges)

    return G


def get_measures(operations, data_set):
    Graph = generate_PPI_network(data_set)
    op_map = {
        "shortest_path": find_shortest_path
    }

    for operation in operations:

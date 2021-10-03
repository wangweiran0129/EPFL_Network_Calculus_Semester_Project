# This code is written by Karim Hadidane .
# For any questions or problems, please contact the author of the code at (karim.hadidane@epfl.ch)

import pbzlib
from graph_transformer import *
import torch
from torch_geometric.data import Data, DataLoader
import pickle


# The inputs are the graph G, unique identifiers for each node, prolongation nodes if its for training
def graph2torch(G, node_ids, prolongation=None, test=False):
    # node features : the paper used 11 input dimension (corrected after asking the author). We use 9 features for
    # training but we add one more case in the feature vector to store the flow corresponding to each prolongation node
    x = torch.zeros((len(G.nodes), 10))

    # outcome of node prediction : NEED TO Fill for training: 1 or 0 for prolongation nodes
    y = torch.zeros(G.number_of_nodes(), dtype=torch.float)

    # As the loss function depends only on prolongation nodes, we use this to mask the other nodes
    prolongation_mask = torch.zeros(G.number_of_nodes(), dtype=torch.bool)

    for node, data in G.nodes(data=True):
        node_idx = node_ids[node]
        node_type = data["ntype"]

        # Adding node features
        # Adding the type in 1-hot encoding
        x[node_idx, node_type - 1] = 1

        # ADDED: A flow of interest is also a flow
        if node_type == NodeType.Flow_oi:
            x[node_idx, 1] = 1
            # y[node_idx] = 1 if worth_prolonging else 0

        # Adding server rates and latencies : position 4, 5
        if node_type == NodeType.Server:
            x[node_idx, 4] = data["rate_server"]
            x[node_idx, 5] = data["latency"]

        # Adding flow rates and bursts : position 6, 7
        if node_type == NodeType.Flow or node_type == NodeType.Flow_oi:
            x[node_idx, 6] = data["rate_flow"]
            x[node_idx, 7] = data["burst"]

        # Identify the prolongation nodes in the mask
        if node_type == NodeType.Prolong:
            prolongation_mask[node_idx] = 1
            x[node_idx, 8] = data["hop"]
            # index of the corresponding flow
            f = node_ids["f_" + node[1:node.index("_")]]
            x[node_idx, 9] = f

    # Each edge is encoded twice
    edge_index = torch.zeros((2, G.number_of_edges() * 2), dtype=torch.long)

    i = 0
    for src, dst in G.edges():
        # Each edge from the undirected graph G is encoded as two directed edges
        edge_index[0, i] = node_ids[src]
        edge_index[1, i] = node_ids[dst]
        i += 1
        edge_index[0, i] = node_ids[dst]
        edge_index[1, i] = node_ids[src]
        i += 1

    # Add target nodes labels : 1 if the prolongation node is the sink , 0 otherwise
    if test:
        for pro_node in prolongation:
            y[node_ids[pro_node]] = 1

    graph = Data(x=x, y=y, edge_index=edge_index, mask=prolongation_mask)

    return graph


# Returns best combination in a list of dictionaries, each dictionary is a mapping between key: flow and value: its sink
def get_best_deborahfp_combination(flow):
    min_delay_bound = flow.deborahfp.delay_bound

    # Maybe many combinations yield the least delay bound
    min_combinations = []

    # Iterate over all explored combinations for DEBORAH FP exhaustive search
    for comb in flow.deborahfp.explored_combination:

        if comb.delay_bound == min_delay_bound:
            min_combinations.append(dict(comb.flows_prolongation))

    return min_combinations


# Returns a list of Data object of the dataset, each one is a graph
def prepare_dataset(filename, train, to_pickle=True):
    graphs = []
    targets = []

    # For each network in the file
    for network in pbzlib.open_pbz(filename):
        # Get the base graph i.e server nodes, flow nodes, and links between them
        G, flow_paths = net2basegraph(network)

        for flow in network.flow:

            # If the flow has been explored using deborah FP
            if flow.HasField("deborahfp"):
                # create version of a graph where the current flow is the flow of interest (Algorithm 2 of the paper)
                G_f, pro_dict, node_ids = prolong_graph(G, flow.id, flow_paths)

                best_combinations = get_best_deborahfp_combination(flow)

                flows_that_can_be_prolonged = set(pro_dict.keys())

                grph = graph2torch(G_f, prolongation=None, node_ids=node_ids, test=False)

                # Append the created graph to the dataset
                graphs.append(grph)

                # worth prolonging if deborah FP gives tighter delay bound than Deborah
                worth_prolonging = flow.deborahfp.delay_bound < flow.deborah.delay_bound
                foi_index = node_ids["f_" + str(flow.id)]

                # Equally optimal targets can exist: this is implementing Equation 11 in the paper
                for comb in best_combinations:
                    possible_targets = []
                    # Prolongation nodes to activate
                    comb_nodes = ["p" + str(k) + "_" + str(v) for k, v in comb.items()]

                    # Some prolongation nodes that need to be activated (for flows that will not be prolonged)
                    # are not included in the mapping given in the dataset
                    # These lines are to mitigate this problem
                    prolonged = set(comb.keys())
                    k = flows_that_can_be_prolonged.difference(prolonged)
                    comb_nodes.extend(list(map(lambda x: "p" + str(x) + "_" + str(flow_paths[x][-1]), k)))

                    y = torch.zeros(G_f.number_of_nodes(), dtype=torch.float)
                    y[foi_index] = worth_prolonging
                    for pro_node in comb_nodes:
                        y[node_ids[pro_node]] = 1
                    possible_targets.append(y)

                targets.append(possible_targets)

    if to_pickle:
        file_name_graphs = "train_graphs.pickle"
        file_name_targets = "train_targets.pickle"
        if not train:
            file_name_graphs= "test_graphs.pickle"
            file_name_targets = "test_targets.pickle"

        # Saving the training graphs in a pickle format
        outfile = open(file_name_graphs, 'wb')
        pickle.dump(graphs, outfile)
        outfile.close()

        # Saving the training targets in a pickle format
        outfile = open(file_name_targets, 'wb')
        pickle.dump(targets, outfile)
        outfile.close()

    print(len(graphs))
    print(len(targets))
    print("finished hopefully with success")

    return graphs, targets


# still a problem to figure out if all combinations are added to the test dataset seperately as independent graphs
def prepare_test_dataset(filename):
    dataset = []

    # For each network in the file
    for network in pbzlib.open_pbz(filename):

        # Get the base graph i.e server nodes, flow nodes, and links between them
        G, flow_paths = net2basegraph(network)

        for flow in network.flow:

            if flow.HasField("deborahfp"):
                flow_id = flow.id
                worth_prolonging = flow.deborahfp.delay_bound < flow.deborah.delay_bound

                G_f, pro_dict, node_ids = prolong_graph(G, flow_id, flow_paths)
                foi_index = node_ids["f_" + str(flow.id)]
                best_combinations = get_best_deborahfp_combination(flow)

                flows_that_can_be_prolonged = set(pro_dict.keys())

                for comb in best_combinations:
                    # Prolongation nodes to activate
                    comb_nodes = ["p" + str(k) + "_" + str(v) for k, v in comb.items()]

                    prolonged = set(comb.keys())

                    k = flows_that_can_be_prolonged.difference(prolonged)
                    comb_nodes.extend(list(map(lambda x: "p" + str(x) + "_" + str(flow_paths[x][-1]), k)))
                    graph = graph2torch(G_f, prolongation=comb_nodes, node_ids=node_ids, test=True)
                    graph.y[foi_index] = worth_prolonging

                    dataset.append(graph)

    print(len(dataset))

    return dataset


def main(filename):
    prepare_dataset(filename)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("input")
    args = p.parse_args()
    main(args.input)

# This code is written by Karim Hadidane.
# For any questions or problems, please contact the author of the code at (karim.hadidane@epfl.ch)
import random

import torch
import torch_geometric.nn as pyg_nn
from torch_geometric.data import Data, DataLoader
from graph_transformer import *
import torch.nn as nn
from prepare_dataset import *
from scipy.sparse import coo_matrix
import numpy as np
import pickle
from tqdm import tqdm
import output_pb2
from pbzlib import write_pbz, open_pbz


# This class implements the gated graph neural network model that will be used for training and prediction
class GGNN(nn.Module):

    # unrolls is the number of recurrence unrolls to perform
    def __init__(self, input_dim, hidden_dim, unrolls, dropout):
        super(GGNN, self).__init__()

        self.unrolls = unrolls

        # Init
        self.init_layer = nn.Sequential(
            nn.Linear(input_dim, hidden_dim, bias=True),
            nn.LeakyReLU())

        # # Memory Unit
        self.gru = pyg_nn.GatedGraphConv(hidden_dim, 2, aggr= "mean")
        #
        # # Edge attention
        self.edge_attention = nn.Sequential(
            nn.Linear(hidden_dim * 2, hidden_dim, bias=True),
            nn.Sigmoid(),
            nn.Linear(hidden_dim, 1, bias=True),
            nn.Sigmoid()
        )
        #
        # # Output Hidden layers
        self.out = nn.Sequential(
            nn.Linear(hidden_dim, hidden_dim, bias=True),
            nn.LeakyReLU(),
            nn.Linear(hidden_dim, hidden_dim, bias=True),
            nn.LeakyReLU()
        )

        # Output layer 1 : for node classification to predict flow prolongations
        self.out1 = nn.Sequential(
            nn.Linear(hidden_dim, 1, bias=True),
            nn.Sigmoid()
        )

        # Output layer 2 : for node classification to predict if worth prolonging
        self.out2 = nn.Sequential(
            nn.Linear(hidden_dim, 1, bias=True),
            nn.Sigmoid()
        )

    def forward(self, data, adj):
        # unpacking the data object
        x, edge_index = data.x[:, :-1], data.edge_index
        
        #print("edge index shape")
        #print(edge_index.shape)
        #print(edge_index)
        #print(edge_index[:,0])
        
        edge_weight= torch.zeros(edge_index.shape[1])
        #print(edge_weight)
        


        x = self.init_layer(x)
        
        #x = self.gru(x, edge_index)
        for i in range(self.unrolls):
            x = self.gru(x, edge_index)
        # Prepare the input for the edge attention mechanism
            d = prepare_edge_attention_input(x)

        # Run the edge attention layer and squeeze to remove the last dimension
            coeffs = self.edge_attention(d).squeeze()
            zero_vec = torch.zeros_like(coeffs)
            coeffs = torch.where(adj == 1, coeffs, zero_vec)
            x = coeffs @ x
        
        
        
        #for i in range(edge_index.shape[1]):
        #    u= edge_index[0][i]
        #    v= edge_index[1][i]
        #    weight= coeffs[u][v]
        #    edge_weight[i]= weight
            
            
        
        
        # Matrix multiplication of the edge coefficients with the feature matrix
        
        

        x = self.out(x)
        out1 = self.out1(x)
        out2 = self.out2(x)
        return out1, out2


# This method is used to train the model
# needs parameters like the learning rate, number of epochs, number of unrolls for the GRU cell
# it needs the training dataset path (pickled after running prepare_dataset)
def train_model(train_graphs, train_targets, test_graphs, test_targets, learning_rate, epochs, dropout=0.2):
    # train_data = DataLoader(train_dataset, batch_size=batch_size)
    outfile = open(train_graphs, 'rb')
    train_graphs = pickle.load(outfile)
    outfile.close()
    outfile = open(train_targets, 'rb')
    train_targets = pickle.load(outfile)
    outfile.close()

    # train_data = DataLoader(train_data, batch_size=batch_size, shuffle=True)

    # Added: test data
    outfile = open(test_graphs, 'rb')
    test_graphs = pickle.load(outfile)
    outfile.close()

    outfile = open(test_targets, 'rb')
    test_targets = pickle.load(outfile)
    outfile.close()

    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")

    # Build the model
    model = GGNN(9, 96, 2, dropout=dropout)

    opt = torch.optim.SGD(model.parameters(), lr=learning_rate)
    model = model.to(device)

    losses_per_epoch = []
    accuracies_per_epoch = []

    # Tells the model that we are training
    model.train()
    for epoch in tqdm(range(epochs)):

        # Reshuffle
        train_data = random.sample(list(zip(train_graphs, train_targets)), len(train_graphs))
        test_data = random.sample(list(zip(test_graphs, test_targets)), len(test_graphs))

        correct_predicted_deepfp = 0
        correct_predicted_deepfp4 = 0

        losses_train = []
        for graph, targets in train_data:
            # Set gradients to zero
            opt.zero_grad()

            # create adjacency matrix
            adj = prepare_adjacency_matrix(graph)

            pred1, pred2 = model(graph.to(device), adj.to(device))

            # Identify nodes of flow of interests
            foi_idx = torch.where(graph.x[:, 2])[0]

            losses = []

            # Implementing equation 11 of the paper: dynamically choosing the target vector from equally good
            # combinations ===> we backward only the minimum loss
            for t in targets:
                # Define the loss function to use: here we use Cross entropy for both classification problems
                criterion1 = nn.BCELoss()
                criterion2 = nn.BCELoss()

                target_foi = torch.index_select(t.to(device), 0, foi_idx)
                output_foi = torch.index_select(pred1.view(-1), 0, foi_idx)

                idxmask = torch.where(graph.mask)[0]

                target_prolongations = torch.index_select(t.to(device), 0, idxmask)
                output_prolongations = torch.index_select(pred2.view(-1), 0, idxmask)

                loss = criterion1(output_prolongations, target_prolongations) + criterion2(output_foi,
                                                                                           target_foi)

                losses.append(loss)

            min_loss = losses[np.argmin(list(map(lambda x: x.item(), losses)))]
            min_loss.backward()
            losses_train.append(min_loss.item())

            # update the parameters
            opt.step()
            
            

        losses_per_epoch.append(np.mean(losses_train))

        # training_acc_deepfp, training_acc_deepfp4 = evaluate(model, train_data, device)

        acc_deepfp, acc_deepfp4 = evaluate(model, test_data, device)
        # accuracies_per_epoch.append(acc)
        print("********* Lossses in training ******")
        

        tqdm.write(
            f"{epoch:3d} | loss={np.mean(losses_train):.2e}  | test_accuracy_deepfp= {acc_deepfp} | test_accuracy_deepfp4= {acc_deepfp4}")
        

    return losses_per_epoch, accuracies_per_epoch, model


def evaluate(model, test_data, device, k=1):
    # Tells the model that we are evaluating
    model.eval()

    correct_predicted_deepfp = 0
    correct_predicted_deepfp4 = 0

    worth_prolonging_graphs = 0

    for graph, targets in test_data:
        with torch.no_grad():
            adj = prepare_adjacency_matrix(graph)
            out1, out2 = model(graph.to(device), adj.to(device))

        accurate_deepfp, accurate_deepfp4 = accurate_deepFP_deepFP4(out1, out2, graph, targets)
        correct_predicted_deepfp = correct_predicted_deepfp + accurate_deepfp
        correct_predicted_deepfp4 = correct_predicted_deepfp4 + accurate_deepfp4
        
        

    # for each graph in the dataset, a correct classification means that the model outputs the correct flow
    # prolongations for all cross flows
    accuracy_deepfp = correct_predicted_deepfp / len(test_data)
    accuracy_deepfp4 = correct_predicted_deepfp4 / len(test_data)

    return accuracy_deepfp, accuracy_deepfp4


def accurate_deepFP_deepFP4(out1, out2, graph, targets):
    # Doing every thing on cpu
    out1 = out1.cpu()
    out2 = out2.cpu()
    graph = graph.cpu()

    foi_idx = torch.where(graph.x[:, 2])[0]
    output_foi = torch.index_select(out1.view(-1), 0, foi_idx)
    predicted_label = 1 if output_foi.item() >= 0.5 else 0

    idxmask = torch.where(graph.mask)[0]
    output_prolongations = torch.index_select(out2.view(-1), 0, idxmask)

    pro_nodes = graph.x[graph.mask]

    prolongations_deepfp = create_output_vector(pro_nodes, output_prolongations, 1)
    prolongations_deepfp4 = create_output_vector(pro_nodes, output_prolongations, 4)

    worth_prolonging = torch.index_select(targets[0], 0, foi_idx).item() == 1

    if worth_prolonging and predicted_label == 1:
        correct_predicted_targets_deepfp = 0
        correct_predicted_targets_deepfp4 = 0

        for t in targets:
            target_prolongations = torch.index_select(t, 0, idxmask).cpu()
            correct_predicted_targets_deepfp = correct_predicted_targets_deepfp + compute_correct_solutions(
                prolongations_deepfp, target_prolongations)
            correct_predicted_targets_deepfp4 = correct_predicted_targets_deepfp4 + compute_correct_solutions(
                prolongations_deepfp4, target_prolongations)

        return correct_predicted_targets_deepfp > 0, correct_predicted_targets_deepfp4 > 0
    elif (not worth_prolonging) and predicted_label == 0:
        return True, True

    return False, False


# This function will compute the number of correct solutions from the flow prolongation combinations returned by
# deepFP or deepFP_k
def compute_correct_solutions(output_prolongations, target_prolongations):
    s = 0

    for alternative in output_prolongations:
        s = s + (torch.sum(target_prolongations == alternative.cpu()).item() == len(target_prolongations))

    return s


# pro_nodes is the prolonagtion nodes : tensor (P, D) where P is the number of prolon
# out_prol is the output of prolongation nodes : tensor (P, 1)
# number_of_combinations is the number of predictions to generate (see DeepFP_2, DeepFP_4)
def create_output_vector(pro_nodes, out_prol, number_of_combinations=1):
    target_prolo = torch.zeros(out_prol.shape[0])
    combinations = []

    flows = set(map(lambda x: int(x.item()), pro_nodes[:, 9]))

    if number_of_combinations == 1:
        for f in flows:
            # Select prolongation nodes related to the flow
            node_indx = torch.nonzero(pro_nodes[:, 9] == f)
            # Select the node with the highest score to be the sink of the prolongation
            r = torch.index_select(out_prol, 0, node_indx.squeeze())

            ind = node_indx[torch.argmax(r, dim=0)].item()

            # Create the output vector
            target_prolo[ind] = 1

        combinations.append(target_prolo)

    else:
        for i in range(number_of_combinations):
            target_prolo = torch.zeros(out_prol.shape[0])
            for f in flows:
                # Select prolongation nodes related to the flow
                node_indx = torch.nonzero(pro_nodes[:, 9] == f)

                # Select the node with the highest score to be the sink of the prolongation
                r = torch.index_select(out_prol, 0, node_indx.squeeze())

                # Generating probabilities from prediction vector output of the model (use softmax here)
                # probs = np.array(nn.Softmax(dim=0)(r).view(len(r)).cpu())
                probs = np.array(r.view(len(r)).cpu())
                probs /= np.sum(probs)

                ind = np.random.choice(a=list(range(len(probs))), p=probs)
                # Create the output vector
                target_prolo[ind + node_indx[0].item()] = 1

            combinations.append(target_prolo)

    # Combinations are the predictions that we generate from the vector of gnn predictions
    return combinations


# Adjacency matrix for the edge attention layer
def prepare_adjacency_matrix(batch):
    nodes_number = batch.x.shape[0]
    adj = coo_matrix(
        (np.ones(batch.edge_index.shape[1]), (batch.edge_index[0].cpu(), batch.edge_index[1].cpu())),
        shape=(nodes_number, nodes_number))
    adj = torch.FloatTensor(np.array(adj.todense())) + torch.eye(nodes_number)

    return adj


def predict_network(network, foi_id, model, output_file="output.pbz"):
    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")

    # Create a base graph
    G, flows_path = net2basegraph(network)

    # prolong the graph with respect to the foi
    G_f, pro_dict, node_ids = prolong_graph(G, foi_id, flows_path)

    graph = graph2torch(G_f, node_ids=node_ids)

    adj = prepare_adjacency_matrix(graph)

    out1, out2 = model(graph.to(device), adj.to(device))

    # HERE a dictionary flow: (start, sink)
    start_sink_dict = {k: [flows_path[k][0], flows_path[k][-1]] for k in flows_path.keys()}

    foi_idx = torch.where(graph.x[:, 2])[0]
    output_foi = torch.index_select(out1.view(-1), 0, foi_idx)
    predicted_label = 1 if output_foi.item() >= 0.5 else 0

    # If the prediction is that FP is not worth then write the same network
    if not predicted_label:
        write_network(network, start_sink_dict, output_file)
        return

    idxmask = torch.where(graph.mask)[0]
    output_prolongations = torch.index_select(out2.view(-1), 0, idxmask)

    pro_nodes = graph.x[graph.mask]

    prolongations_deepfp = create_output_vector(pro_nodes, output_prolongations, 1)
    sinks = idxmask[torch.where(prolongations_deepfp[0])[0]]

    # print(graph.x)
    # Create a dictionary nodeids and flow id
    inv_map = {v: k for k, v in node_ids.items()}

    z = [inv_map[x] for x in np.array(sinks)]
    to_be_prolonged = {get_flowid_from_prolongation_node_name(k): get_serverid_from_prolongation_node_name(k) for k in
                       z}

    for flow, server in to_be_prolonged.items():
        start_sink_dict[flow][1] = server

    write_network(network, start_sink_dict, output_file)

    return graph


def write_network(network, flows_start_sink, filename):
    objs = [output_pb2.Network(id=1)]

    for s in network.server:
        p = objs[0].server.add()
        p.id = s.id
        p.rate = s.rate
        p.latency = s.latency

    for f in network.flow:
        p = objs[0].flow.add()
        p.id = f.id
        p.rate = f.rate
        p.burst = f.burst
        p.start = flows_start_sink[f.id][0]
        p.sink = flows_start_sink[f.id][1]

    with write_pbz(filename, "output.descr") as w:
        for obj in objs:
            w.write(obj)

    ## Start reading
    print("************************************")
    for network in open_pbz("output.pbz"):
        print(network)


# Notes:
# check if adding self loops is needed

# Method to prepare the input for the edge attention layer by concatenating neighboring nodes feature vectors
def prepare_edge_attention_input(feature_matrix):
    # Number of nodes
    N = feature_matrix.shape[0]

    # Number of features
    D = feature_matrix.shape[1]

    # This is a way to create every possible combination of nodes in the graph
    # Selecting the right nodes i.e neighbors is done later using the adjacency matrix
    t1 = feature_matrix.repeat_interleave(N, dim=0)
    t2 = feature_matrix.repeat(N, 1)
    new = torch.cat([t1, t2], dim=1)

    return new.view(N, N, 2 * D)


def get_flowid_from_prolongation_node_name(s):
    flow = int(s[s.index("_") - 1])
    return flow


def get_serverid_from_prolongation_node_name(s):
    server = int(s[s.index("_") + 1])
    return server


def main(args):
    losses = train_model(train_dataset=args.dataset_path, learning_rate=args.learning_rate, epochs=args.epochs,
                         batch_size=args.batch_size)


if __name__ == "__main__":
    p = argparse.ArgumentParser(
        description="Train a gated graph neural network descriped in \"Tightening network calculus delay bounds by predicting flow prolongations in the FIFO Analysis\" by Fabien Geyer, Alexander Scheffler and Steffen Bondorf")
    p.add_argument("--dataset-path", type=str, help="Training dataset path", default="train_data.pickle")
    p.add_argument("--epochs", type=int, default=5, help="Number of epochs used for training")
    p.add_argument("--learning-rate", type=float, default=5e-4, help="Learning rate for Adam Optimizer")
    p.add_argument("--batch-size", type=int, default=32, help="Batch size")

    args = p.parse_args()
    main(args)

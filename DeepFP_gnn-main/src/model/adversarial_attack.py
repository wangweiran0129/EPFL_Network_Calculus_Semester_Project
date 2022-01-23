# This code is written by Weiran Wang
# For any questions or problems, please contact the author of code at (weiran.wang@epfl.ch)

from importlib.metadata import requires
import pickle
import random
import torch
import torchvision
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from train_model import *

def fgsm_update(graph_feature_matrix, graph_feature_matrix_grad, update_max_norm):
    """
    Use Fast Gradient Sign Attack (FGSM) to add adversarial attack on graph feature matrix
    It mainly follows the function: x = x + eps * sign(x.grad)
    :param graph_feature_matrix: the feature matrix of the graph, i.e., graph.x
    :param graph_feature_matrix_grad: the gradient of graph feature matrix
    :param update_max_norm: epsilons value in the attack
    :return: the graph feature matrix after the adversarial attack
    """
    # graph: <class 'torch_geometric.data.data.Data'>

    # Collect the element-wise sign of the data gradient
    sign_graph_feature_matrix_grad = graph_feature_matrix_grad.sign()

    # Create the perturbed graph
    perturbed_graph_feature_matrix = graph_feature_matrix + update_max_norm * sign_graph_feature_matrix_grad

    # Adding clipping to maintain [0,1] graph
    perturbed_graph_feature_matrix = torch.clamp(perturbed_graph_feature_matrix, 0, 1)

    # graph["x"] = perturbed_graph_feature_matrix

    return perturbed_graph_feature_matrix


def evaluate_attack(model, test_graphs, test_targets, update_max_norm):
    """
    Evaluate the adversarial attack
    :param model: the GGNN model we trained before
    :param test_graphs: the test graphs pickle
    :param test_targets: the test targets pickle
    :param update_max_norm: the epsilons list
    :return: the accuracies of deepfp and deepfp4
    """
    # Added: test data
    outfile = open(test_graphs, 'rb')
    test_graphs = pickle.load(outfile)
    outfile.close()

    outfile = open(test_targets, 'rb')
    test_targets = pickle.load(outfile)
    outfile.close()

    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    
    # shuffle the testset
    test_data = random.sample(list(zip(test_graphs, test_targets)), len(test_graphs))

    correct_predicted_deepfp = 0
    correct_predicted_deepfp4 = 0

    # Loop over all examples in test set in batches
    for graph, targets in test_data:
        # print("x before attack : ", graph["x"])
        # x is the node feature matrix in torch_geometric.data.data.Data
        node_feature_matrix = graph.x
        # 4 -> rate server, 5 -> latency
        # 6 -> rate flow, 7 -> burst
        # node_feature_matrix = graph["x"][:,4:8]
        # store the data except the indices of 4,5,6,7
        node_feature_matrix_first_part = graph.x[:,:4]
        node_feature_matrix_second_part = graph.x[:,8:]
        node_feature_matrix.requires_grad = True

        # create adjacency matrix
        adjacency = prepare_adjacency_matrix(graph)
        pred1, pred2 = model(graph.to(device), adjacency.to(device))

        # Identify nodes of flow of interests
        foi_idx = torch.where(graph.x[:,2])[0]

        losses = []

        # forward pass
        for t in targets:
            # Define the loss function to use: CrossEntropy
            criterion1 = nn.BCELoss()
            criterion2 = nn.BCELoss()

            target_foi = torch.index_select(t.to(device), 0, foi_idx)
            output_foi = torch.index_select(pred1.view(-1), 0, foi_idx)

            # indices of prolongation nodes
            idxmask = torch.where(graph.mask)[0]

            # prolongations nodes:
            target_prolongations = torch.index_select(t.to(device), 0, idxmask)
            output_prolongations = torch.index_select(pred2.view(-1), 0, idxmask)

            loss = criterion1(output_prolongations, target_prolongations) + criterion2(output_foi, target_foi)
            losses.append(loss)

        min_loss = losses[np.argmin(list(map(lambda x: x.item(), losses)))]
        # zero all existing gradients
        model.zero_grad()
        min_loss.backward()

        # add adversarial attacks to node feature by using fgsm
        perturbed_matrix = fgsm_update(node_feature_matrix, node_feature_matrix.grad, update_max_norm)
        graph["x"] = perturbed_matrix
        # restore the origin number except for indices 4,5,6,7
        graph.x[:,:4] = node_feature_matrix_first_part
        graph.x[:,8:] = node_feature_matrix_second_part
        # graph["x"][:,4:8] = perturbed_matrix
    
        # create adjacency matrix
        perturbed_adjacency = prepare_adjacency_matrix(graph)
        
        # re-classify the perturbed batch
        # now the graph node features have been attacked 
        perturbed_pred1, perturbed_pred2 = model(graph.to(device), perturbed_adjacency.to(device))

        accurate_deepfp, accurate_deepfp4 = accurate_deepFP_deepFP4(perturbed_pred1, perturbed_pred2, graph, targets)
        correct_predicted_deepfp = correct_predicted_deepfp + accurate_deepfp
        correct_predicted_deepfp4 = correct_predicted_deepfp4 + accurate_deepfp4
    
    accuracy_deepfp = correct_predicted_deepfp / len(test_data)
    accuracy_deepfp4 = correct_predicted_deepfp4 / len(test_data)
    
    # return original_accuracy
    return accuracy_deepfp, accuracy_deepfp4


def accuracy_visualization(epsilons, accuracy_deepfp, accuracy_deepfp4):
    """
    visualize the accurcy vs epsilons
    :param: epsilons list
    :param: accuracies value
    """
    fig, ax = plt.subplots()
    ax.plot(epsilons, accuracy_deepfp, marker='o', label='deepfp')
    ax.plot(epsilons, accuracy_deepfp4, marker='*', label='deepfp4')
    plt.gca().xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.set_xlabel('epsilons')
    ax.set_ylabel('accuracy')
    ax.set_title('Accuracy vs Epsilons in FGSM')
    ax.legend()

    plt.show()


def main():
    # load model
    model = torch.load("/Users/wangweiran/Desktop/SemesterProject/EPFL_Network_Calculus_Semester_Project/DeepFP_gnn-main/src/model/ggnn_learning_rate_5e4.pt", map_location=torch.device("cpu"))

    accuracy_deepfp = []
    accuracy_deepfp4 = []
    # 0 represents the model performance on the original test set.
    epsilons = [0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3]

    test_graphs = "/Users/wangweiran/Desktop/SemesterProject/EPFL_Network_Calculus_Semester_Project/DeepFP_gnn-main/src/data/processed_serialized/test_graphs.pickle"
    test_targets = "/Users/wangweiran/Desktop/SemesterProject/EPFL_Network_Calculus_Semester_Project/DeepFP_gnn-main/src/data/processed_serialized/test_targets.pickle"

    # Run test for each epsilon
    for eps in tqdm(epsilons):
        print("eps : ", eps)
        acc_deepfp, acc_deepfp4 = evaluate_attack(model, test_graphs, test_targets, eps)
        accuracy_deepfp.append(acc_deepfp)
        accuracy_deepfp4.append(acc_deepfp4)
    
    print("deepfp accuracy : ", accuracy_deepfp)
    print("deepfp4 accuracy : ", accuracy_deepfp4)

    # visualize the result
    accuracy_visualization(epsilons, accuracy_deepfp, accuracy_deepfp4)


if __name__ == "__main__":
    main()
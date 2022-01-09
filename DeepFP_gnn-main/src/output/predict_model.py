import sys
sys.path.append(r'/Users/wangweiran/Desktop/SemesterProject/EPFL_Network_Calculus_Semester_Project/DeepFP_gnn-main')
sys.path.insert(0, "../src")
from src.data.graph_transformer import *
from src.data.prepare_dataset import *
from src.model.train_model import *
from src.output.output_pb2 import *
from src.model.gnn import *
import re
import torch
import torch.nn as nn
from pbzlib import write_pbz, open_pbz


def predict_network(network, foi_id, model, output_file="output.pbz"):
    """
    A method that uses the model trained to predict new network configuration and write it in proto file format
    :param network: network parameters
    :param foi_id: the flow of interest
    :param model: the model
    :param output_file: output file to generate
    :return: 
    """
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
    # print("predicted label : ", predicted_label)

    # The graph with a single server don not need to prolong
    if len(network.server) == 1:
        print("There is only one server in this topology, so the flow cannot be prolonged")
        write_network(network, start_sink_dict, output_file)
        return

    # If the prediction is that FP is not worth then write the same network
    if not predicted_label:
        write_network(network, start_sink_dict, output_file)
        print("The network is not worth to prolong")
        return

    print("This network can be prolonged")
    idxmask = torch.where(graph.mask)[0]
    output_prolongations = torch.index_select(out2.view(-1), 0, idxmask)

    pro_nodes = graph.x[graph.mask]

    prolongations_deepfp = create_output_vector(pro_nodes, output_prolongations, 1)
    sinks = idxmask[torch.where(prolongations_deepfp[0])[0]]

    # Create a dictionary nodeids and flow id
    inv_map = {v: k for k, v in node_ids.items()}

    z = [inv_map[x] for x in np.array(sinks)]
    print("what's z?", z)
    to_be_prolonged = {get_flowid_from_prolongation_node_name(k): get_serverid_from_prolongation_node_name(k) for k in z}
    print("to_be_prolonged : ", to_be_prolonged)

    for flow, server in to_be_prolonged.items():
        start_sink_dict[flow][1] = server

    write_network(network, start_sink_dict, output_file)

    print("start_sink_dic : ", start_sink_dict)

    return graph


def write_network(network, flows_start_sink, filename):
    """
    A method that writes the network generated into a protobuf file according to the output.descr description file
    :param network: the network parameters
    :param flows_start_sink: the flows path
    :param filename: output filename
    """
    objs = [Network(id=1)]

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

    with write_pbz(filename, "/Users/wangweiran/Desktop/SemesterProject/EPFL_Network_Calculus_Semester_Project/DeepFP_gnn-main/src/output/output.descr") as w:
        for obj in objs:
            w.write(obj)


def get_flowid_from_prolongation_node_name(s):
    # flow = int(re.search(r"\d+", s).group())
    flow = int(s[s.index("_") - 1])
    return flow


def get_serverid_from_prolongation_node_name(s):
    # server_temp = re.search(r"_\d+", s).group()
    # server = int(re.search(r"\d+", server_temp).group())
    server = int(s[s.index("_") + 1])
    return server

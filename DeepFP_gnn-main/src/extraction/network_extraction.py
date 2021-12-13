# This code is written by Weiran Wang
# For any questions or problems, please contact the author of code at (weiran.wang@epfl.ch)

import pbzlib
import csv
import sys
import re
sys.path.append(r'/Users/wangweiran/Desktop/SemesterProject/EPFL_Network_Calculus_Semester_Project/DeepFP_gnn-main')

from src.data.graph_transformer import *

def info_extraction(i, network):
    
    '''
    For the server_info.csv files, the column names are
    topology_id, server_id, server_rate, server_latency

    :param i : the topology id
    :network: the network information
    :return: server & flow csv files
    '''

    # create the graph from the dataset
    G, flows_paths = net2basegraph(network)

    # i will also be the topology_id
    topology_id = str(i)

    # extract server information
    servers_id = []
    servers_rates = []
    servers_latency = []

    for node in G.nodes:
        #exclude the flow node
        if 'f_' in node:
            break
        servers_id.append(re.search(r'\d+', node).group())
        servers_rates.append(G.nodes.data()[node]['rate_server'])
        servers_latency.append(G.nodes.data()[node]['latency'])
    
    server_size = len(servers_id)

    # extract flow information
    flows_id = []
    flows_rates = []
    flows_bursts = []
    flows_src = []
    flows_dest = []
    flow_of_interest = []

    for node in G.nodes:
        # exclude the server node
        if 's_' in node:
            continue
        # ('f_0', {'ntype': <NodeType.Flow: 2>, 'rate_flow': 0.00031810534748519687, 'burst': 0.968261581700407})
        flows_id.append(re.search(r'\d+', node).group())
        flows_rates.append(G.nodes.data()[node]['rate_flow'])
        flows_bursts.append(G.nodes.data()[node]['burst'])
    
    flow_size = len(flows_id)
    
    # src and dest may change according to different flow of interest
    # sometimes, the index of srcs are larger than dests
    # sometimes, the opposite

    # for the basic graph
    for path in flows_paths:
        flow_src = flows_paths[path][0]
        flow_dest = flows_paths[path][-1]
        '''
        if flow_src > flow_dest:
            flow_src, flow_dest = flow_dest, flow_src
        '''
        flows_src.append(flow_src)
        flows_dest.append(flow_dest)
        # for the original graph, foi is -1
        flow_of_interest.append(-1)
    
    # for the prolonged graph
    # calculate the number of foi in one network
    # foi is the flow of interest e.g.[0, 5, 6, 9, 12, 13, 14, 16]
    foi = []
    for flow in network.flow:
        if flow.HasField('deborahfp'):
            foi.append(flow.id)
    
    # the basic information of a network is the same
    # the dest will change according to different foi
    for foi_index in foi:
        G_f, pro_dict, node_ids = prolong_graph(G, foi_index, flows_paths)
        # print("foi_index : ", foi_index)
        # This part will be changed later
        for path in flows_paths:
            # the src doesn't change
            flow_src = flows_paths[path][0]
            # flow in pro_dict should change their destinations
            if path in pro_dict:
                flow_dest = re.match('.*?([0-9]+)$', pro_dict[path][-1]).group(1)
            else:
                flow_dest = flows_paths[path][-1]
            flows_src.append(flow_src)
            flows_dest.append(flow_dest)
            # for the foi, it will be labeled as 1
            if path == foi_index:
                flow_of_interest.append(1)
            else:
                flow_of_interest.append(-1)

    # create a csv file for server
    path_server = '/Users/wangweiran/Desktop/SemesterProject/EPFL_Network_Calculus_Semester_Project/server_info/'
    filename_server = 'topology' + topology_id + '_server.csv'
    print('filename = ', path_server + filename_server)

    # write server data into csv file
    with open(path_server + filename_server, 'w') as csvfile:
        server_csv = csv.writer(csvfile)

        # write the columns names
        server_csv.writerow(['topology_id', 'server_id', 'server_rate', 'server_latency'])
        for j in range(server_size):
            server_csv.writerow([topology_id, servers_id[j], servers_rates[j], servers_latency[j]])

    # create a csv file for flow
    path_flow = '/Users/wangweiran/Desktop/SemesterProject/EPFL_Network_Calculus_Semester_Project/flow_info/'
    filename_flow = 'topology' + topology_id + '_flow.csv'
    print('file_name = ', path_flow + filename_flow)

    # write flow data into csv file
    with open(path_flow + filename_flow, 'w') as csvfile:
        flow_csv = csv.writer(csvfile)

        # write the columns names
        flow_csv.writerow(['topology_id', 'network_id', 'flow_id', 'flow_rate', 'flow_burst', 'flow_src', 'flow_dest', 'flow_of_interest'])
        for k in range(len(foi)+1):
            for j in range(flow_size):
                flow_csv.writerow([topology_id, k, flows_id[j], flows_rates[j], flows_bursts[j], flows_src[j + flow_size * k], flows_dest[j + flow_size * k], flow_of_interest[j + flow_size * k]])


def test(path):
    i = 0
    for i, network in enumerate(pbzlib.open_pbz(path)):
        G, flows_paths = net2basegraph(network)
        print(i)
        if (i == 4):
            print(G)
            print(flows_paths)
            break


def main(path):

    # define the name of the csv file
    i = 0
    for i, network in enumerate(pbzlib.open_pbz(path)):
        print("topology : ", i)
        info_extraction(i, network)
    # test(path)


if __name__ == "__main__":
    path = '/Users/wangweiran/Desktop/SemesterProject/EPFL_Network_Calculus_Semester_Project/DeepFP_gnn-main/src/extraction/source-sink.pbz'
    main(path);

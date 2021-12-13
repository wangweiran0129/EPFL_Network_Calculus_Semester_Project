# This code is written by Weiran Wang
# For any questions or problems, please contact the author of code at (weiran.wang@epfl.ch)

import csv
import networkx as nx
import enum
import sys
from itertools import islice
sys.path.append(r'/Users/wangweiran/Desktop/SemesterProject/EPFL_Network_Calculus_Semester_Project/DeepFP_gnn-main')

from networkx.algorithms import flow

from src.data.graph_transformer import prolong_graph

def source_sink_info(servers):
    
    '''
    Construct a source-sink network servers and flows information
    The shape of only depends on the number of server(s)
    servers rates, R = 10 Mb/s
    servers latency, T = 0.001 s
    flows bursts, b = 1kb
    flow rates depneds on the load r = (R*0.5)/servers
    flow of interest (foi) crosses all servers
    :param servers : the index of servers size
    '''
    topology_id = str(servers)
    server_rate = 10
    server_latency = 0.001
    servers_size = servers + 1

    # creaet a csv file for server
    path_server = '/Users/wangweiran/Desktop/SemesterProject/EPFL_Network_Calculus_Semester_Project/server_info/'
    filename_server = 'source-sink' + topology_id + '_server.csv'
    print('filename = ', path_server + filename_server)

    # write server data into csv file
    with open(path_server + filename_server, 'w') as csvfile:
        server_csv = csv.writer(csvfile)

        # write the columns names
        server_csv.writerow(['topology_id', 'server_id', 'server_rate', 'server_latency'])
        for i in range(servers_size):
            server_csv.writerow([topology_id, i, server_rate, server_latency])

    #topology_id = str(servers)
    flows_number = 2 * servers_size - 1
    network_id = []
    flows_bursts = 10
    flows_rates = 5 / servers_size
    flows_src = []
    flows_dest = []
    foi = []

    # define the src & dest of a flow
    # 1st category : source is the first server and the destination changes
    for dest in range(servers_size):
        flows_src.append(0)
        flows_dest.append(dest)
    
    # 2nd category : destination is the last server and sources change
    # exclude src = 0, dest = servers_size-1 which has been calculated before in the first category
    for src in range(1,servers_size):
        flows_src.append(src)
        flows_dest.append(servers_size-1)

    # set -1 for foi in the base graph
    for i in range(flows_number):
        foi.append(-1)

    # define the flow of interest
    for networkID in range(flows_number):
        network_id.append(networkID)
        for flow in range(flows_number):
            if flow == networkID:
                foi.append(1)
            else:
                foi.append(-1)

    # create a csv file for flow
    path_flow = '/Users/wangweiran/Desktop/SemesterProject/EPFL_Network_Calculus_Semester_Project/flow_info/'
    filename_flow = 'source-sink' + topology_id + '_flow.csv'
    print('file_name = ', path_flow + filename_flow)

    # write flow data into csv file
    with open(path_flow + filename_flow, 'w') as csvfile:
        flow_csv = csv.writer(csvfile)

        # write the columns names
        flow_csv.writerow(['topology_id', 'network_id', 'flow_id', 'flow_rate', 'flow_burst', 'flow_src', 'flow_dest', 'flow_of_interest'])
        # loop for networks
        for i in range(flows_number+1):
            # loop for flows
            for j in range(flows_number):
                flow_csv.writerow([topology_id, i, j, flows_rates, flows_bursts, flows_src[j], flows_dest[j], foi[j+flows_number*i]])


NodeType = enum.IntEnum('NodeType', [
    'Server', 'Flow', 'Flow_oi', 'Prolong'
])


def source_sink_construction(serverfile, flowfile, server_size):
    '''
    A method to read the network information from a csv file and convert it into a network
    :param server: the name of the csv file which stores the server information
    :param flow; the name of the csv file which stores the flow information
    :param server_size: the size of servers in this topology
    :return: the base graph of the form of networkx graph object, a dictionary of flow paths of the form flow_id: [servers]
    '''
    server_size = server_size + 1
    flow_size = 2 * server_size - 1

    # Initialize empty graph
    G = nx.Graph()

    servers_list = []
    servers_edges = []

    # Initialize nodes features dictionaries
    servers_rates = {}
    servers_latencies = {}
    flows_rates = {}
    flows_bursts = {}
    flows_src = {}
    flows_dest = {}
    flows_paths = {}
    flow_paths_list = []

    # Dictionary to save the path of each flow
    flow_paths = {}
    
    # Add server nodes to the grpah
    with open(serverfile, 'r') as csvServer:
        reader = csv.reader(csvServer)
        for line in islice(reader, 1, server_size+1):
            server_id = 's_' + line[1]
            servers_list.append(server_id)
            servers_rates[server_id] = line[2]
            servers_latencies[server_id] = line[3]
            G.add_node(server_id, ntype=NodeType.Server)
    
    # Add edges between servers node
    for i in range(server_size-1):
        # there is no edge in one single server
        if(server_size == 1):
            continue
        else:
            servers_edges.append((servers_list[i], servers_list[i+1]))
    
    # Add flow nodes to the graph
    with open(flowfile, 'r') as csvFlow:
        reader = csv.reader(csvFlow)
        for line in islice(reader, 1, flow_size+1):
            flow_paths_list.clear()
            flow_id = 'f_' + str(line[2])
            G.add_node(flow_id, ntype=NodeType.Flow)
            flows_rates[flow_id] = line[3]
            flows_bursts[flow_id] = line[4]
            flows_src[flow_id] = int(line[5])
            flows_dest[flow_id] = int(line[6])
            for path in range(flows_src[flow_id], flows_dest[flow_id]+1):
                flow_paths_list.append(path)
            flow_paths[int(line[2])] = flow_paths_list
            servers_in_path = list(map(lambda x: (flow_id, 's_' + str(x)), list(flow_paths[int(line[2])])))
            G.add_edges_from(servers_in_path)
            flow_paths[int(line[2])] = list(flow_paths[int(line[2])])
    
    # Add attributes of the nodes
    nx.set_node_attributes(G, servers_rates, name="rate_server")
    nx.set_node_attributes(G, flows_rates, name='rate_flow')
    nx.set_node_attributes(G, servers_latencies, name='latency')
    nx.set_node_attributes(G, flows_bursts, name='burst')

    return G, flow_paths

    
def main():
    server_num = 25
    for i in range(server_num):
        # source_sink_info(server)
        serverfile = '/Users/wangweiran/Desktop/SemesterProject/EPFL_Network_Calculus_Semester_Project/server_info/source-sink' + str(i) + '_server.csv'
        flowfile = '/Users/wangweiran/Desktop/SemesterProject/EPFL_Network_Calculus_Semester_Project/flow_info/source-sink' + str(i) + '_flow.csv'
        G, flow_paths = source_sink_construction(serverfile, flowfile, i)
        print('G = ', G)
        print('flow_paths = ', flow_paths)
        if i == 1:
            G_f, pro_dict, node_ids = prolong_graph(G, 0, flow_paths)
            print('G_f = ', G_f)
            print('pro_dict = ', pro_dict)
            print('node_ids = ', node_ids)
            break


if __name__ == "__main__":
    main()

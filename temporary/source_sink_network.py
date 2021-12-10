# This code is written by Weiran Wang
# For any questions or problems, please contact the author of code at (weiran.wang@epfl.ch)

import csv
import networkx as nx

def source_sink_server(servers):
    
    '''
    Construct a source-sink network servers
    The shape of only depends on the number of server(s)
    
    servers rates, R = 10 Mb/s
    servers latency, T = 0.001 s

    :servers : the number of server(s)
    '''
    topology_id = str(servers)
    server_rate = 10
    server_latency = 0.001

    # creaet a csv file for server
    path_server = '/Users/wangweiran/Desktop/SemesterProject/EPFL_Network_Calculus_Semester_Project/server_info/'
    filename_server = 'source-sink' + topology_id + '_server.csv'
    print('filename = ', path_server + filename_server)

    # write server data into csv file
    with open(path_server + filename_server, 'w') as csvfile:
        server_csv = csv.writer(csvfile)

        # write the columns names
        server_csv.writerow(['topology_id', 'server_id', 'server_rate', 'server_latency'])
        for i in range(servers):
            server_csv.writerow([topology_id, i+1, server_rate, server_latency])


def source_sink_flow(servers):
    '''
    flows bursts, b = 1kb
    flow rates depneds on the load r = (R*0.5)/servers
    flow of interest (foi) crosses all servers
    '''
    topology_id = str(servers)
    flows_number = 2 * servers - 1
    network_id = []
    flows_bursts = 10
    flows_rates = 5 / servers
    flows_src = []
    flows_dest = []
    foi = []

    # define the src & dest of a flow
    for src in range(1, servers+1):
        for dest in range(src, servers+1):
            flows_src.append(src)
            flows_dest.append(dest)

    # set -1 for foi in the base graph
    for i in range(flows_number):
        foi.append(-1)

    # define the flow of interest
    for networkID in range(1, flows_number+1):
        network_id.append(networkID)
        for flow in range(1, flows_number+1):
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
                flow_csv.writerow([topology_id, i, j+1, flows_rates, flows_bursts, flows_src[j], flows_dest[j], foi[j+flows_number*i]])
 

def main():
    server_num = 25
    for server in range(server_num):
        source_sink_server(server + 1)
        source_sink_flow(server + 1)


if __name__ == "__main__":
    main()

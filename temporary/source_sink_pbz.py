from pbzlib import write_pbz, open_pbz
import sys
sys.path.append(r'/Users/wangweiran/Desktop/SemesterProject/EPFL_Network_Calculus_Semester_Project/DeepFP_gnn-main/')
from src.extraction.tests.source_sink_pb2 import Network, Server
from tests import source_sink_pb2

def source_sink_network(server_number):

    # objs store the network information
    objs = []

    for topology_id in range(server_number):
        server_size = topology_id + 1
        flow_size = 2 * server_size - 1
        flows_src = []
        flows_dest = []

        objs.append(source_sink_pb2.Network(id=topology_id))

        # add server information
        for i in range(server_size):
            server = objs[topology_id].server.add()
            server.id = i
            server.rate = 10
            server.latency = 0.001
        
        # define the src & dest of a flow
        # 1st category : source is the first server and the destination changes
        for dest in range(server_size):
            flows_src.append(0)
            flows_dest.append(dest)
            
        # 2nd category : destination is the last server and sources change
        # exclude src = 0, dest = servers_size-1 which has been calculated before in the first category
        for src in range(1,server_size):
            flows_src.append(src)
            flows_dest.append(server_size-1)
    
        # add flow information
        for i in range(flow_size):
            flow = objs[topology_id].flow.add()
            flow.id = i
            flow.rate = 5/server_size
            flow.burst = 10
            # add flow path
            source = flows_src[i]
            destination = flows_dest[i]
            for FlowPath in range(source, destination+1):
                flow.path.append(FlowPath)
            # add flow of interest
            flow.deborahfp.delay_bound = 1

    with write_pbz("source-sink.pbz", "tests/source_sink.descr") as w:
        for obj in objs:
            w.write(obj)


def main():
    server_number = 25
    for i in range(server_number):
        source_sink_network(server_number)
    
    for network in open_pbz("source-sink.pbz"):
        print(network)

if __name__ == "__main__":
    main()
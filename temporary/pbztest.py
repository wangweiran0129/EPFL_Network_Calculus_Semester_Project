from pbzlib import write_pbz, open_pbz
import sys
sys.path.append(r'/Users/wangweiran/Desktop/SemesterProject/EPFL_Network_Calculus_Semester_Project/DeepFP_gnn-main/')
from src.extraction.tests.source_sink_pb2 import Network, Server
from tests import source_sink_pb2

def source_sink_network(server_number):

    objs = [source_sink_pb2.Network(id=0)]
    server = objs[0].server.add()
    server.id = 0
    server.rate = 10
    server.latency = 0.001
    flow = objs[0].flow.add()
    flow.id = 0
    flow.rate = 5
    flow.burst = 10
    flow.path.append(0)


    objs.append(source_sink_pb2.Network(id=1))
    server = objs[1].server.add()
    server.id = 0
    server.rate = 10
    server.latency = 0.001
    server = objs[1].server.add()
    server.id = 1
    server.rate = 10
    server.latency = 0.001
    flow = objs[1].flow.add()
    flow.id = 0
    flow.rate = 5
    flow.burst = 10
    flow.path.append(0)
    flow = objs[1].flow.add()
    flow.id = 1
    flow.rate = 5
    flow.burst = 10
    flow.path.append(0)
    flow.path.append(1)
    flow = objs[1].flow.add()
    flow.id = 2
    flow.rate = 5
    flow.burst = 10
    flow.path.append(1)

    with write_pbz("output.pbz", "tests/source_sink.descr") as w:
        for obj in objs:
            w.write(obj)
    
    for network in open_pbz("output.pbz"):
        print(network)

def main():
    server_number = 2
    source_sink_network(server_number)

if __name__ == "__main__":
    main()
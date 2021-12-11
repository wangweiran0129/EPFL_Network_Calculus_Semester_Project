from pbzlib import write_pbz, open_pbz
import sys
sys.path.append(r'/Users/wangweiran/Desktop/SemesterProject/EPFL_Network_Calculus_Semester_Project/DeepFP_gnn-main/')
from src.extraction.tests.source_sink_pb2 import Network, Server
from tests import source_sink_pb2, messages_pb2

def main():
    objs = [source_sink_pb2.Network(id=0)]

    network = source_sink_pb2.Network()
    network.id = 0

    server = objs[0].server.add()
    server.id = 0
    server.rate = 10
    server.latency = 0.001

    flow = objs[0].flow.add()
    flow.id = 0
    flow.rate = 5
    flow.burst = 10
    flow.path.append(0)
    
    with write_pbz("output.pbz", "tests/source_sink.descr") as w:
        for obj in objs:
            w.write(obj)

    for network in open_pbz("output.pbz"):
        print(network)


if __name__ == "__main__":
    # example()
    main()
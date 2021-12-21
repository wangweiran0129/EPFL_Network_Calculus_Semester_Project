import sys
sys.path.append(r'/Users/wangweiran/Desktop/SemesterProject/EPFL_Network_Calculus_Semester_Project/DeepFP_gnn-main')
from src.model.train_model import *
from src.output.output_pb2 import *
from src.model.gnn import *
from src.output.predict_model import *
from pbzlib import write_pbz, open_pbz


def source_sink_prolongation():
    model = GGNN(9, 96, unrolls=2)
    # print("model : ", model)

    main_path = "/Users/wangweiran/Desktop/SemesterProject/EPFL_Network_Calculus_Semester_Project/Network_analysis/"
    pbz_path = main_path + "source_sink_prolongation_pbz/"
    for topology_id, base_network in enumerate(pbzlib.open_pbz(main_path + "extraction/source-sink.pbz")):
        print("\n----- topology_id : ", topology_id, "-----" )
        pbz_name = "source_sink_tandem_" + str(topology_id) + ".pbz"
        prolongation_pbz = pbz_path + pbz_name
        flow_of_interest = topology_id
        prolongation_prediction = predict_network(network=base_network, foi_id=flow_of_interest, model=model, output_file=prolongation_pbz)
        '''
        if topology_id == 2:
            print("base_network : ", base_network)
            print("topology_id : ", topology_id)
            break
        '''
    # print("prolongation_prediction: ", prolongation_prediction, "\n")
    # prolonged_network=next(pbzlib.open_pbz(prolongation_pbz))
    # print("prolonged_network: ", prolonged_network, "\n")

def main():
    source_sink_prolongation()


if __name__ == "__main__":
    main()
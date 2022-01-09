import torch
import sys
sys.path.append("../")
from output.predict_model import *
import pbzlib



def source_sink_prolongation():
    model = torch.load("/Users/wangweiran/Desktop/SemesterProject/EPFL_Network_Calculus_Semester_Project/DeepFP_gnn-main/src/model/ggnn.pt", map_location=torch.device('cpu'))
    main_path = "/Users/wangweiran/Desktop/SemesterProject/EPFL_Network_Calculus_Semester_Project/"
    pbz_path = main_path + "Network_Information_and_Analysis/source_sink_prolongation_pbz/"
    for topology_id, base_network in enumerate(pbzlib.open_pbz(main_path + "DeepFP_gnn-main/src/network_analysis/source_sink_tandem_generation/source-sink.pbz")):
        print("\n----- topology_id : ", topology_id, "-----" )
        pbz_name = "source_sink_tandem_" + str(topology_id) + ".pbz"
        prolongation_pbz = pbz_path + pbz_name
        flow_of_interest = topology_id
        prolongation_prediction = predict_network(network=base_network, foi_id=flow_of_interest, model=model, output_file=prolongation_pbz)


def main():
    source_sink_prolongation()


if __name__ == "__main__":
    main()
# EPFL_Network_Calculus_Semester_Project

## Project Background
Network Calculus offers a framework for the worst-cast end-to-end delay performance analysis of data commuication. Generally, a tighter bound refers to a more complex analysis method. Thus, finding a trade-off between computation effort and tightness is worth researching. In other words, a network analysis method with good tightness and fast execution should be investigated.  
An entirely different approach was recently presented with the Flow Prolongation (FP) feature, which can create a nested tandem (sequence of servers) in a different way before facing the original network cutting-problem.

## Project Target
- The first step is to use GNN, already implemented, and combine it with NetCal/DNC (DOBORAH integrated) tool to compute delay bounds for flows after predicting prolongation setting. This enables us to compute delay bounds for new networks rather than those of the [existing dataset](https://github.com/fabgeyer/dataset-rtas2021). Then, we implement the same method on source sink tandem networks and compared them to those of Polynomial-size Linear Program (PLP).  
- The second step is to find small changes in the inputs that make GNN unstable, i.e., adding adverserial attacks to modify inputs of GNN minimally, and then observe unexpected and pessimistic delay bounds.  

## Project Environment and Prerequisites
The whole project is run on my Macbook Pro (Intel Chip) and a EPFL server.  
- IDE: Visual Studio Code (Version: 1.63.1)  
- Python: 3.8.5  
- Java: 16.0.2  
- Apache Maven: 3.8.3  
- [NetCal/DNC](https://github.com/NetCal/DNC): v2.7  
- [IBM ILOG CPLEX Optimization Studio](https://www.ibm.com/docs/en/icos/20.1.0?topic=cplex-setting-up-gnulinuxmacos): 20.1.0  
(For Mac users, there will be an existing installation error if your Mac appearance is DARK)
- Python Packages Information can be found at [requirements.txt](https://github.com/wangweiran0129/EPFL_Network_Calculus_Semester_Project/blob/master/DeepFP_gnn-main/requirements.txt)

## EPFL SCITAS Server Configuration
For EPFL servers ONLY!  
Similar configuration should also be set properly. Since students are not the admin of EPFL servers, it is therefore a good idea to put Java configuration into a user-define file under the home path, i.e., add the following three lines into the .bashrc file.
```
export JAVA_HOME=/home/weirwang/jdk-16.0.2
export PATH=$JAVA_HOME/bin:$PATH
export CLASSPATH=.:$JAVA_HOME/lib/dt.jar:$JAVA_HOME/lib/tools.jar
```
For python environment, the system default version is Python 2.7.5 and the pip version is out of fashion. However, students cannot upgrade pip due to the non-admin status. Consequently, a python virtual environment is recommended. For jupyter notebook configuration, please refer to the [official website](https://scitas-data.epfl.ch/confluence/display/DOC/Using+Jupyter). After configuration, to enter the virtual environment, please use command:
```
$source opt/$SYS_TYPE/venv-gcc/bin/activation
```
Then, python package can be installed and upgraded according to specific needs. To activate jupyter notebook, please use command:
```
$jupyter notebook --ip="$(hostname -s).epfl.ch"
```

## Codes Description
`network_extraction.py`  
Extract all the necessary network information from the original topology graphs and the corresponding post-GNN prolonged topology graphs. These network information will be used in NetCal/DNC for the network analysis.  
For servers, they are topology id, server id, server rate and server latency. For flows, they are topology id, network id, flow id, flow rate, flow burst, flow source, flow destination and flow of interest. The difference between topology id and network id is that topology id refers to one big topology which include several network ids depending on different flows of interest. The flow_of_interest is defined as a binary variable where 1 refers to a flow of interest in this network, otherwise -1.  
This python script will generate csv files for servers and flows separately.  

`source_sink_network.py`  
Generate network information for source sink tandem networks analysis and output csv files similar above. The speciality of a source sink tandem network is that each flow either starts at the first server or ends at the last server, so there will be 2n-1 flows with n servers in one topology.  

`source_sink.proto`  
Define the data structure of source sink tandem network.  

`Makefile`  
Generate [source_sink_pb2.py](https://github.com/wangweiran0129/EPFL_Network_Calculus_Semester_Project/blob/master/DeepFP_gnn-main/src/extraction/tests/source_sink_pb2.py) and [source_sink.descr](https://github.com/wangweiran0129/EPFL_Network_Calculus_Semester_Project/blob/master/DeepFP_gnn-main/src/extraction/tests/source_sink.descr) which are needed for the following `source_sink_pbz.py` to write `source-sink.pbz`.  

`source_sink_pbz.py`  
Similar functions as `source_sink_network.py`, but will output network information in a .pbz format. More information about Protocol Buffers can be found [here](https://developers.google.com/protocol-buffers). This format is required for the input of GNN. The output of this python script will be `source-sink.pbz`.  

`visualization.py`  
Visualize the results after running `Source_Sink_Analysis.java`(or the DNC.jar integrated in this project). The abscissa is the number of servers and the ordinates are the delays. (It waits to be updated after the full results from source_sink_analysis.java)

`TopologyTest.java`  
For test purpose, but is written in a straightforward way, where the tutorial of using NetCal/DNC can be learned quickly from this java code.

`Topology1.java`  
Reads topology information from server and flow csv files and calculates Total Flow Analysis (TFA), Separated Flow Analysis (SFA), PMOO Analysis, Tandem Matching Analysis (TMA) and Least Upper Delay Bound (LUDB) analysis.  
It is worth mentioning here that the stand-alone tool [DEBORAH](http://cng1.iet.unipi.it/wiki/index.php/Deborah) failed. Therefore, the NetCal/DNC developers only published their implementation of LUDB. In order to run LUDB-FF, [IBM ILOG CPLEX Optimization Studio](https://www.ibm.com/docs/en/icos/20.1.0?topic=cplex-setting-up-gnulinuxmacos) is needed. Then, the variables path_to_cplex and path_to_lp_dir in `AnalysisConfig.java` should be set. The former should be the path of the CPLEX binary path, e.g., ```path_to_cplex="/Applications/CPLEX_Studio201/cplex/bin/x86-64_osx/"```, and the directory is the LP files will be stored temporarily. It's the same for the EPFL server setting if LUDB-LL is needed.  

`Source_Sink_Analysis.java`  
Focuses on the source sink tandem network analysis compared to `Topology1.java` but the working pattern of these two java codes are the same.

`DNC.jar`  
The final package file used to aggregate all the Java classes files which can be ran directly by the command ```java -jar DNC.jar``` in the terminal. This command can be used both locally and on a server.  

## Usage Specification
1. In source_sink folder,a `Makefile` can be found. It makes preparation for source sink tandem network by compile source_sink.proto. It will output two files: source_sink_pb2.py and source_sink.descr which will be used in the generation of source sink tandem network in the next step.
```
make
```
2. In `extraction` folder, run source_sink_pbz.py script to generate the source sink tandem networks. The default number of server is 10 while it will be changed later to let users define the number of servers.
```
python3 source_sink_pbz.py
```
3. Train the GNN model. This part of code is written by Hadidane Karim and modified by me. To see the process of training model. Please refer to this [jupyter notebook](https://github.com/wangweiran0129/EPFL_Network_Calculus_Semester_Project/blob/master/DeepFP_gnn-main/notebooks/train_model.ipynb).
4. Get the information of prolonged network by running:
```
python3 source_sink_prolongation.py
```
5. Extract the information of network, which is a necessary step for the following NetCal/DNC analysis:
```
python3 network_extraction.py
```
6. For user convenience, I integrate all the java code into a DNC.jar format. It can be run by:
```
java -jar DNC.jar
```
7. Visualization: to see the comparasion of delay bound before and after the flow prolongation, run:
```
python3 visualization.py
```

## Attention
This repo is still under testing case. The final version will come out at the beginning of February.  
The path in codes are based on my own Mac absolute path. If there is an error referring to unfound files. Please check the path of the files closely. It will also be fixed in the final version.


## Disclaimer and Special Acknowledgement
- Project Student: Weiran Wang (weiran.wang@epfl.ch)
- Ph.D. Advisor: Tabatabaee Hossein (hossein.tabatabaee@epfl.ch)
- Supervisor: Prof. Le boudec Jean-Yves (jean-yves.leboudec@epfl.ch)
- Special Acknowledgement to:  
Hadidane Karim (karim.hadidane@epfl.ch)  
Bondorf Steffen (Bondorf.Steffen@ruhr-uni-bochum.de)  
Alexander Scheffler (Alexander.Scheffler@ruhr-uni-bochum.de)
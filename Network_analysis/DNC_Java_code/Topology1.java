// This code is written by Weiran Wang
// For any questions or problems, please contact the author of the code at (weiran.wang@epfl.ch)

package org.networkcalculus.dnc.semesterproject;

// network calculus analysis
import java.io.*;
import java.util.*;
import org.networkcalculus.dnc.AnalysisConfig;
import org.networkcalculus.dnc.curves.ArrivalCurve;
import org.networkcalculus.dnc.curves.Curve;
import org.networkcalculus.dnc.curves.ServiceCurve;
import org.networkcalculus.dnc.network.server_graph.Flow;
import org.networkcalculus.dnc.network.server_graph.Server;
import org.networkcalculus.dnc.network.server_graph.ServerGraph;
import org.networkcalculus.dnc.tandem.analyses.PmooAnalysis;
import org.networkcalculus.dnc.tandem.analyses.SeparateFlowAnalysis;
import org.networkcalculus.dnc.tandem.analyses.TandemMatchingAnalysis;
import org.networkcalculus.dnc.tandem.analyses.TotalFlowAnalysis;
import org.networkcalculus.num.Num;
import java.net.URLEncoder;

// csv package

public class Topology1 {

    public Topology1() {
    }

    public static void main(String[] args) {

        System.out.println("Hello to my first prolonged topology");
        System.out.println();

        Topology1 topology = new Topology1();

        // The number of server files and flow files info are the same
        String path = "/Users/wangweiran/Desktop/SemesterProject/EPFL_Network_Calculus_Semester_Project/server_info/";
        int filenumber = get_file_number(path);

        System.out.println("How many files in server folder : " + filenumber);

        // Generate the network and analyze the network
        try {
            for (int topology_id = 0; topology_id < filenumber; topology_id++) {
                topology.run(topology_id);
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
    }


    public static int get_file_number(String path) {
        int fileCount = 0;
        File d = new File(path);
        File list[] = d.listFiles();
        for (int i = 0; i<list.length; i++) {
            if(list[i].isFile()){
                fileCount++;
            }
        }
        return fileCount;
    }


    public void run(int topology_id) throws Exception {

        // define lists for server information
        ArrayList<String> server_id = new ArrayList<String>();
        ArrayList<String> server_rate = new ArrayList<String>();
        ArrayList<String> server_latency = new ArrayList<String>();
        ArrayList<ServiceCurve> service_curve = new ArrayList<ServiceCurve>();

        // define lists for flow information
        ArrayList<String> flow_id = new ArrayList<String>();
        ArrayList<String> arrival_rate = new ArrayList<String>();
        ArrayList<String> flow_burst = new ArrayList<String>();
        ArrayList<String> flow_src = new ArrayList<String>(); // flow_src stores all the sources in one network
        ArrayList<String> flow_dest = new ArrayList<String>(); // flow_dest stores all the destinations in one network
        ArrayList<String> flow_of_interest_temp = new ArrayList<String>(); // foi ids are stored here
        ArrayList<ArrivalCurve> arrival_curve = new ArrayList<ArrivalCurve>();

        // Define the list to store the analysis result;
        ArrayList<Num> TFA_DB = new ArrayList<Num>();
        ArrayList<Num> TFA_BB = new ArrayList<Num>();
        ArrayList<Num> SFA_DB = new ArrayList<Num>();
        ArrayList<Num> SFA_BB = new ArrayList<Num>();
        ArrayList<Num> PMOO_DB = new ArrayList<Num>();
        ArrayList<Num> PMOO_BB = new ArrayList<Num>();
        ArrayList<Num> TMA_DB = new ArrayList<Num>();
        ArrayList<Num> TMA_BB = new ArrayList<Num>();

        // Java reading csv file
        String path_server = "/Users/wangweiran/Desktop/SemesterProject/EPFL_Network_Calculus_Semester_Project/server_info/";
        String filename_server = "topology" + topology_id + "_server.csv";
        String csvServer = path_server + filename_server;
        System.out.println("server information stored in : " + csvServer);
        BufferedReader brServer = null;
        String line = "";
        String csvSplitBy = ",";

        // basic server information
        try {
            brServer = new BufferedReader(new FileReader(csvServer));
            brServer.readLine(); // read the title line
            while ((line = brServer.readLine()) != null) {
                // use comma as separator
                String[] server = line.split(csvSplitBy);
                server_id.add(server[1]);
                server_rate.add(server[2]);
                server_latency.add(server[3]);
                double serverRate = Double.parseDouble(server[2]);
                double serverLatency = Double.parseDouble(server[3]);
                service_curve.add(Curve.getFactory().createRateLatency(serverRate, serverLatency));
            }
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        } finally {
            if (brServer != null) {
                try {
                    brServer.close();
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
        }

        String path_flow = "/Users/wangweiran/Desktop/SemesterProject/EPFL_Network_Calculus_Semester_Project/flow_info/";
        String filename_flow = "topology" + topology_id + "_flow.csv";
        String csvFlow = path_flow + filename_flow;
        System.out.println("flow information stored in : " + csvFlow);
        BufferedReader brFlow = null;
        int flow_counter = 0;
        int topology_size = 1;
        int flow_size = 0;
        // basic flow information before flow prolongation
        // need to get the # flow in the topology
        // and how many topologies I should have in total
        try {
            brFlow = new BufferedReader(new FileReader(csvFlow));
            brFlow.readLine(); // read the title line
            while ((line = brFlow.readLine()) != null) {
                // use comma as separator
                String[] flow = line.split(csvSplitBy);
                // flow_beginning is the first flow in the network, i.e., f0
                int flow_beginning = Integer.parseInt(flow[2]);
                if (flow_beginning != flow_counter) {
                    topology_size = topology_size + 1;
                    flow_size = flow_counter;
                    flow_counter = 0;
                }
                flow_id.add(flow[2]);
                arrival_rate.add(flow[3]);
                flow_burst.add(flow[4]);
                double arrivalRate = Double.parseDouble(flow[3]);
                double flowBurst = Double.parseDouble(flow[4]);
                arrival_curve.add(Curve.getFactory().createTokenBucket(arrivalRate, flowBurst));
                int foi = Integer.parseInt(flow[7]);
                if (foi == 1) {
                    flow_of_interest_temp.add(flow[2]);
                }
                flow_counter = flow_counter + 1;
            }
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        } finally {
            if (brFlow != null) {
                try {
                    brFlow.close();
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
        }
        System.out.println("topology size : " + topology_size);
        System.out.println("flow size : " + flow_size);

        // print server information
        System.out.println();
        System.out.println("--- SERVER INFORMATION --- ");
        for (int i = 0; i < server_id.size(); i++) {
            System.out.println("server id : " + server_id.get(i));
            System.out.println("server rate : " + server_rate.get(i));
            System.out.println("server latency : " + server_latency.get(i));
            System.out.println();
        }

        // print flow information
        System.out.println();
        System.out.println("--- FLOW INFORMATION ---");
        for (int i = 0; i < flow_size; i++) {
            System.out.println("flow id : " + flow_id.get(i));
            System.out.println("arrival rate : " + arrival_rate.get(i));
            System.out.println("flow burst : " + flow_burst.get(i));
            System.out.println();
        }

        // we need to prolong the topology and connect the flow into the topology
        try {
            brFlow = new BufferedReader(new FileReader(csvFlow));
            brFlow.readLine();
            // Skip the original topology information
            for (int i = 0; i < flow_size; i++) {
                brFlow.readLine();
            }
            while ((line = brFlow.readLine()) != null) {
                // use comma as separator
                String[] flow = line.split(csvSplitBy);
                // dest_id vary according to different foi(s)
                int src_id = Integer.parseInt(flow[5]);
                flow_src.add(server_id.get(src_id));
                int dest_id = Integer.parseInt(flow[6]);
                flow_dest.add(server_id.get(dest_id));
            }
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        } finally {
            if (brFlow != null) {
                try {
                    brFlow.close();
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
        }

        // System.out.println("flow src : " + flow_src);
        // System.out.println("flow dest : " + flow_dest);
        
        // src > dest -> backward -> 0
        // src < dest -> forward -> 1
        // set a flag to detect the whether it's feed-forward or back-forward
        int flag = 0;
        for (int i = 0; i < flow_size; i++) {
            int src = Integer.parseInt(flow_src.get(i));
            int dest = Integer.parseInt(flow_dest.get(i));
            if (src == dest) {
                continue;
            }
            // src > dest -> backward -> 0
            else if (src > dest) {
                flag = 0;
                break;
            }
            // src < dest -> forward -> 1
            else {
                flag = 1;
                break;
            }
        }

        // src & dest may differ according to the prolonged topology and foi
        try {
            
            // j is the order of topology_size
            for (int j = 1; j < topology_size; j++) {

                System.out.println();
                System.out.println("----- The " + j + " Network ----- ");
                System.out.println();

                // Create a Network Topology & Add Server Curve
                ServerGraph sg = new ServerGraph();
                ArrayList<Server> serverId = new ArrayList<Server>();
                AnalysisConfig configuration = new AnalysisConfig();
                configuration.enforceMaxSC(AnalysisConfig.MaxScEnforcement.GLOBALLY_ON);
                configuration.enforceMaxScOutputRate(AnalysisConfig.MaxScEnforcement.GLOBALLY_ON);

                // add service curve into server
                for (int k = 0; k < server_id.size(); k++) {
                    serverId.add(sg.addServer(service_curve.get(k)));
                }

                if (flag == 1) {
                    System.out.println("The servers are connected by forward");
                    // addTurn : Connect the Servers
                    for (int k = 0; k < server_id.size()-1; k++) {
                        sg.addTurn(serverId.get(k), serverId.get(k+1));
                    }
                }

                if (flag == 0) {
                    System.out.println("The servers are connected by backward");
                    // addTurn : Connect the Servers
                    for (int k = server_id.size()-1; k > 0; k--) {
                        sg.addTurn(serverId.get(k), serverId.get(k-1));
                    }
                }

                // addFlow : Connect the Flow with Servers
                for (int i = 0; i < flow_size; i++) {
                    String fi = "f" + flow_id.get(i);
                    int index = i + flow_size*(j-1);
                    int src = Integer.parseInt(flow_src.get(index));
                    int dest = Integer.parseInt(flow_dest.get(index));
                    if (src == dest){
                        // System.out.println(fi + "; src : " + src + "; index : " + index);
                        sg.addFlow(fi, arrival_curve.get(i), serverId.get(src));
                    }
                    else {
                        // System.out.println(fi + "; src : " + src + "; dest : " + dest + "; index : " + index);
                        sg.addFlow(fi, arrival_curve.get(i), serverId.get(src), serverId.get(dest));
                    }
                }

                // when i reaches the flow size, we will do the analysis
                // Set the Flow of Interest
                int foi_no = Integer.parseInt(flow_of_interest_temp.get(j - 1));
                Flow flow_of_interest = sg.getFlow(foi_no);
                System.out.println("Flow of interest : " + flow_of_interest.toString());
                System.out.println();

                // Analyze the Network
                // TFA
                System.out.println("--- Total Flow Analysis ---");
                System.out.println();
                TotalFlowAnalysis tfa = new TotalFlowAnalysis(sg, configuration);
                try {
                    tfa.performAnalysis(flow_of_interest);
                    Num tfa_db = tfa.getDelayBound();
                    Num tfa_bb = tfa.getBacklogBound();
                    TFA_DB.add(tfa_db);
                    TFA_BB.add(tfa_bb);
                    System.out.println("delay bound : " + tfa_db);
                    System.out.println(" per server : " + tfa.getServerDelayBoundMapString());
                    System.out.println("backlog bound : " + tfa_bb);
                    System.out.println(" per server : " + tfa.getServerBacklogBoundMapString());
                    System.out.println("alpha per server: " + tfa.getServerAlphasMapString());
                } catch (Exception e) {
                    System.out.println("TFA analysis failed");
                    e.printStackTrace();
                }

                System.out.println();

                // SFA
                System.out.println("--- Separated Flow Analysis ---");
                System.out.println();
                SeparateFlowAnalysis sfa = new SeparateFlowAnalysis(sg, configuration);
                try {
                    sfa.performAnalysis(flow_of_interest);
                    System.out.println("e2e SFA SCs : " + sfa.getLeftOverServiceCurves());
                    System.out.println(" per server : " + sfa.getServerLeftOverBetasMapString());
                    System.out.println("xtx per server : " + sfa.getServerAlphasMapString());
                    Num sfa_db = sfa.getDelayBound();
                    Num sfa_bb = sfa.getBacklogBound();
                    SFA_DB.add(sfa_db);
                    SFA_BB.add(sfa_bb);
                    System.out.println("delay bound : " + sfa_db);
                    System.out.println("backlog bound : " + sfa_bb);
                } catch (Exception e) {
                    System.out.println("SFA analysis failed");
                    e.printStackTrace();
                }

                System.out.println();

                // PMOO
                System.out.println("--- PMOO Analysis ---");
                System.out.println();
                PmooAnalysis pmoo = new PmooAnalysis(sg, configuration);
                try {
                    pmoo.performAnalysis(flow_of_interest);
                    System.out.println("e2e PMOO SCs : " + pmoo.getLeftOverServiceCurves());
                    System.out.println("xtx per server : " + pmoo.getServerAlphasMapString());
                    Num pmoo_db = pmoo.getDelayBound();
                    Num pmoo_bb = pmoo.getBacklogBound();
                    PMOO_DB.add(pmoo_db);
                    PMOO_BB.add(pmoo_bb);
                    System.out.println("delay bound : " + pmoo_db);
                    System.out.println("backlog bound : " + pmoo_bb);
                } catch (Exception e) {
                    System.out.println("PMOO analysis failed");
                    e.printStackTrace();
                }

                System.out.println();

                // TMA
                System.out.println("--- Tandem Matching Analysis ---");
                System.out.println();
                TandemMatchingAnalysis tma = new TandemMatchingAnalysis(sg, configuration);
                try {
                    tma.performAnalysis(flow_of_interest);
                    System.out.println("e2e TMA SCs : " + tma.getLeftOverServiceCurves());
                    System.out.println("xtx per server : " + tma.getServerAlphasMapString());
                    Num tma_db = tma.getDelayBound();
                    Num tma_bb = tma.getBacklogBound();
                    TMA_DB.add(tma_db);
                    TMA_BB.add(tma_bb);
                    System.out.println("delay bound : " + tma_db);
                    System.out.println("backlog bound : " + tma_bb);
                } catch (Exception e) {
                    System.out.println("TMA analysis failed");
                    e.printStackTrace();
                }

                System.out.println();
                System.out.println();

            }

            // Write the analysis results into a csv file

            String fileName = "topology" + topology_id + "_analysis.csv";
            String filePath = "/Users/wangweiran/Desktop/SemesterProject/EPFL_Network_Calculus_Semester_Project/network_analysis/";
            File writename = new File(filePath + fileName);
            writename.createNewFile();
            BufferedWriter out = new BufferedWriter(new FileWriter(writename));
            out.write("topology_id, flow of interest, TFA_DB, TFA_BB, SFA_DB, SFA_BB, PMOO_DB, PMOO_BB, TMA_DB, TMA_BB \r\n");
            for (int temp = 0; temp < flow_of_interest_temp.size(); temp++) {
                out.write(topology_id + ", " + flow_of_interest_temp.get(temp) 
                + ", " + TFA_DB.get(temp).toString() + ", " + TFA_BB.get(temp).toString()
                + ", " + SFA_DB.get(temp).toString() + ", " + SFA_BB.get(temp).toString()
                + ", " + PMOO_DB.get(temp).toString() + ", " + PMOO_BB.get(temp).toString()
                + ", " + TMA_DB.get(temp).toString() + ", " + TMA_BB.get(temp).toString() + "\r\n");
            }
            out.flush();
            out.close();

        } catch (Exception e) {
            e.printStackTrace();
            throw new RuntimeException(e);
        }

    }

}
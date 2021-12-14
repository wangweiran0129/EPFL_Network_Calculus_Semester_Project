// This code is written by Weiran Wang
// For any questions or problems, please contact the author of the code at (weiran.wang@epfl.ch)
package org.networkcalculus.dnc.semesterproject;

import org.networkcalculus.dnc.AnalysisConfig;
import org.networkcalculus.dnc.curves.ArrivalCurve;
import org.networkcalculus.dnc.curves.Curve;
import org.networkcalculus.dnc.curves.ServiceCurve;
import org.networkcalculus.dnc.network.server_graph.Flow;
import org.networkcalculus.dnc.network.server_graph.Server;
import org.networkcalculus.dnc.network.server_graph.ServerGraph;
import org.networkcalculus.dnc.tandem.analyses.FIFOTandemAnalysis;
import org.networkcalculus.dnc.tandem.analyses.PmooAnalysis;
import org.networkcalculus.dnc.tandem.analyses.SeparateFlowAnalysis;
import org.networkcalculus.dnc.tandem.analyses.TandemMatchingAnalysis;
import org.networkcalculus.dnc.tandem.analyses.TotalFlowAnalysis;

// network calculus analysis
public class TopologyTest {
    
    public TopologyTest() {
    }

    public static void main(String[] args) {

        TopologyTest topology = new TopologyTest();

        try {
            topology.run();
        } catch  (Exception e) {
            e.printStackTrace();
        }
    }

    public void run() throws Exception {
        
        double serverRate0 = 0.652298881523439;
        double serverLatency0 = 0.5384246880613892;

        double serverRate1 = 0.1981014874337499;
        double serverLatency1 = 0.14038693814170045;
        
        double serverRate2 = 0.5586898242041254;
        double serverLatency2 = 0.4173048057735682;
        
        double serverRate3 = 0.6704675095373693;
        double serverLatency3 = 0.027387597179620027;
        
        double serverRate4 = 0.8781174306508484;
        double serverLatency4 = 0.20445225303549885;

        double serverRate5 = 0.6852195043049478;
        double serverLatency5 = 0.4191945092867837;

        double serverRate6 = 0.9168613286559002;
        double serverLatency6 = 0.0836230021541644;

        double serverRate7 = 0.49305959663062193;
        double serverLatency7 = 0.9499381437063239;

        double flowRate0 = 0.00031810534748519687;
        double flowBurst0 = 0.968261581700407;

        double flowRate1 = 0.0002356023629769934;
        double flowBurst1 = 0.6716541045589675;

        double flowRate2 = 0.0002147600101331183;
        double flowBurst2 = 0.7654851009830297;

        double flowRate3 = 4.265214627621127e-05;
        double flowBurst3 = 0.5098102725477676;

        double flowRate4 = 5.4208469518974745e-05;
        double flowBurst4 = 0.05991769065467234;

        double flowRate5 = 0.00019615591506642638;
        double flowBurst5 = 0.619955718245456;

        double flowRate6 = 0.0002990411363998393;
        double flowBurst6 = 0.06653647655939565;

        double flowRate7 = 0.0002076633728669982;
        double flowBurst7 = 0.7368114725919939;

        double flowRate8 = 0.00036800388295322284;
        double flowBurst8 = 0.29451154896603016;

        double flowRate9 = 0.0002122877824810558;
        double flowBurst9 = 0.2814767649785799;

        double flowRate10 = 7.598042900864155e-05;
        double flowBurst10 = 0.5888399423961517;

        double flowRate11 = 0.00024023069238016657;
        double flowBurst11 = 0.8288458020909559;

        double flowRate12 = 0.00012818888142791194;
        double flowBurst12 = 0.6707887868192388;

        double flowRate13 = 7.888425808616545e-05;
        double flowBurst13 = 0.4267010094225365;

        double flowRate14 = 0.0003186680735496977;
        double flowBurst14 = 0.5724885178140751;

        double flowRate15 = 9.463213208324156e-05;
        double flowBurst15 = 0.04441793533582761;

        double flowRate16 = 0.00032486826267223533;
        double flowBurst16 = 0.6563158575668357;

        ServiceCurve service_curve0 = Curve.getFactory().createRateLatency(serverRate0, serverLatency0);
        ServiceCurve service_curve1 = Curve.getFactory().createRateLatency(serverRate1, serverLatency1);
        ServiceCurve service_curve2 = Curve.getFactory().createRateLatency(serverRate2, serverLatency2);
        ServiceCurve service_curve3 = Curve.getFactory().createRateLatency(serverRate3, serverLatency3);
        ServiceCurve service_curve4 = Curve.getFactory().createRateLatency(serverRate4, serverLatency4);
        ServiceCurve service_curve5 = Curve.getFactory().createRateLatency(serverRate5, serverLatency5);
        ServiceCurve service_curve6 = Curve.getFactory().createRateLatency(serverRate6, serverLatency6);
        ServiceCurve service_curve7 = Curve.getFactory().createRateLatency(serverRate7, serverLatency7);

        ServerGraph sg = new ServerGraph();
        AnalysisConfig configuration = new AnalysisConfig();
        

        Server s0 = sg.addServer(service_curve0);
        Server s1 = sg.addServer(service_curve1);
        Server s2 = sg.addServer(service_curve2);
        Server s3 = sg.addServer(service_curve3);
        Server s4 = sg.addServer(service_curve4);
        Server s5 = sg.addServer(service_curve5);
        Server s6 = sg.addServer(service_curve6);
        Server s7 = sg.addServer(service_curve7);

        configuration.enforceMaxSC(AnalysisConfig.MaxScEnforcement.GLOBALLY_ON);
        configuration.enforceMaxScOutputRate(AnalysisConfig.MaxScEnforcement.GLOBALLY_ON);

        sg.addTurn(s7, s6);
        sg.addTurn(s6, s5);
        sg.addTurn(s5, s4);
        sg.addTurn(s4, s3);
        sg.addTurn(s3, s2);
        sg.addTurn(s2, s1);
        sg.addTurn(s1, s0);

        ArrivalCurve arrival_curve0 = Curve.getFactory().createTokenBucket(flowRate0, flowBurst0);
        ArrivalCurve arrival_curve1 = Curve.getFactory().createTokenBucket(flowRate1, flowBurst1);
        ArrivalCurve arrival_curve2 = Curve.getFactory().createTokenBucket(flowRate2, flowBurst2);
        ArrivalCurve arrival_curve3 = Curve.getFactory().createTokenBucket(flowRate3, flowBurst3);
        ArrivalCurve arrival_curve4 = Curve.getFactory().createTokenBucket(flowRate4, flowBurst4);
        ArrivalCurve arrival_curve5 = Curve.getFactory().createTokenBucket(flowRate5, flowBurst5);
        ArrivalCurve arrival_curve6 = Curve.getFactory().createTokenBucket(flowRate6, flowBurst6);
        ArrivalCurve arrival_curve7 = Curve.getFactory().createTokenBucket(flowRate7, flowBurst7);
        ArrivalCurve arrival_curve8 = Curve.getFactory().createTokenBucket(flowRate8, flowBurst8);
        ArrivalCurve arrival_curve9 = Curve.getFactory().createTokenBucket(flowRate9, flowBurst9);
        ArrivalCurve arrival_curve10 = Curve.getFactory().createTokenBucket(flowRate10, flowBurst10);
        ArrivalCurve arrival_curve11 = Curve.getFactory().createTokenBucket(flowRate11, flowBurst11);
        ArrivalCurve arrival_curve12 = Curve.getFactory().createTokenBucket(flowRate12, flowBurst12);
        ArrivalCurve arrival_curve13 = Curve.getFactory().createTokenBucket(flowRate13, flowBurst13);
        ArrivalCurve arrival_curve14 = Curve.getFactory().createTokenBucket(flowRate14, flowBurst14);
        ArrivalCurve arrival_curve15 = Curve.getFactory().createTokenBucket(flowRate15, flowBurst15);
        ArrivalCurve arrival_curve16 = Curve.getFactory().createTokenBucket(flowRate16, flowBurst16);

        sg.addFlow(arrival_curve0, s5, s1);
        sg.addFlow(arrival_curve1, s2);
        sg.addFlow(arrival_curve2, s7, s3);
        sg.addFlow(arrival_curve3, s1, s0);
        sg.addFlow(arrival_curve4, s2);
        sg.addFlow(arrival_curve5, s4, s0);
        sg.addFlow(arrival_curve6, s3, s1);
        sg.addFlow(arrival_curve7, s2);
        sg.addFlow(arrival_curve8, s1, s0);
        sg.addFlow(arrival_curve9, s5, s0);
        sg.addFlow(arrival_curve10, s7, s0);
        sg.addFlow(arrival_curve11, s0);
        sg.addFlow(arrival_curve12, s4, s2);
        sg.addFlow(arrival_curve13, s3, s1);
        sg.addFlow(arrival_curve14, s2, s0);
        sg.addFlow(arrival_curve15, s5);
        sg.addFlow(arrival_curve16, s3, s0);

        Flow flow_of_interest = sg.getFlow(0);
        System.out.println("Flow of interest : " + flow_of_interest.toString());
        System.out.println("FOI getpath : " + flow_of_interest.getPath());
        System.out.println();

        // Analyze the network
            // TFA
            System.out.println();
            System.out.println("--- Total Flow Analysis ---");
            // If no analysis configuration is given, the defaults are used
            TotalFlowAnalysis tfa = new TotalFlowAnalysis(sg, configuration);

            try {
                tfa.performAnalysis(flow_of_interest);
                System.out.println("delay bound     : " + tfa.getDelayBound());
                // System.out.println("     per server : " + tfa.getServerDelayBoundMapString());
                System.out.println("backlog bound   : " + tfa.getBacklogBound());
                // System.out.println("     per server : " + tfa.getServerBacklogBoundMapString());
                // System.out.println("alpha per server: " + tfa.getServerAlphasMapString());
            } catch (Exception e) {
                System.out.println("TFA analysis failed");
                e.printStackTrace();
            }

            System.out.println();

            // SFA
            System.out.println();
            System.out.println("--- Separated Flow Analysis ---");
            // If no analysis configuration is given, the defaults are used
            SeparateFlowAnalysis sfa = new SeparateFlowAnalysis(sg, configuration);

            try {
                sfa.performAnalysis(flow_of_interest);
                // System.out.println("e2e SFA SCs     : " + sfa.getLeftOverServiceCurves());
                // System.out.println("     per server : " + sfa.getServerLeftOverBetasMapString());
                // System.out.println("xtx per server  : " + sfa.getServerAlphasMapString());
                System.out.println("delay bound     : " + sfa.getDelayBound());
                System.out.println("backlog bound   : " + sfa.getBacklogBound());
            } catch (Exception e) {
                System.out.println("SFA analysis failed");
                e.printStackTrace();
            }

            System.out.println();

            // PMOO
            System.out.println();
            System.out.println("--- PMOO Analysis ---");
            // If no analysis configuration is given, the defaults are used
            PmooAnalysis pmoo = new PmooAnalysis(sg, configuration);

            try {
                pmoo.performAnalysis(flow_of_interest);
                // System.out.println("e2e PMOO SCs    : " + pmoo.getLeftOverServiceCurves());
                // System.out.println("xtx per server  : " + pmoo.getServerAlphasMapString());
                System.out.println("delay bound     : " + pmoo.getDelayBound());
                System.out.println("backlog bound   : " + pmoo.getBacklogBound());
            } catch (Exception e) {
                System.out.println("PMOO analysis failed");
                e.printStackTrace();
            }
            
            System.out.println();
            
            // TMA
            System.out.println();
            System.out.println("--- Tandem Matching Analysis ---");
            // If no analysis configuration is given, the defaults are used
            TandemMatchingAnalysis tma = new TandemMatchingAnalysis(sg, configuration);

            try {
                tma.performAnalysis(flow_of_interest);
                // System.out.println("e2e TMA SCs     : " + tma.getLeftOverServiceCurves());
                // System.out.println("xtx per server  : " + tma.getServerAlphasMapString());
                System.out.println("delay bound     : " + tma.getDelayBound());
                System.out.println("backlog bound   : " + tma.getBacklogBound());
            } catch (Exception e) {
                System.out.println("TMA analysis failed");
                e.printStackTrace();
            }

            // LUDB-FF
            System.out.println();
            System.out.println("--- LUDB-FF Analysis ---");
            FIFOTandemAnalysis ludb_ff = new FIFOTandemAnalysis(sg, configuration);

            try {
                ludb_ff.performAnalysis(flow_of_interest);
                System.out.println("delay bound     : " + ludb_ff.getDelayBound());
            } catch (Exception e) {
                System.out.println("LUDB-FF analysis failed");
                e.printStackTrace();
            }

            System.out.println();
            System.out.println();

    }
}

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package gephitk;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import org.gephi.graph.api.*;
import org.gephi.io.importer.api.Container;
import org.gephi.io.importer.api.EdgeDirectionDefault;
import org.gephi.io.importer.api.EdgeMergeStrategy;
import org.gephi.io.importer.api.ImportController;
import org.gephi.project.api.ProjectController;
import org.gephi.project.api.Workspace;
import org.gephi.io.processor.plugin.AppendProcessor;
import org.gephi.io.processor.plugin.DefaultProcessor;
import org.gephi.statistics.plugin.*;
import org.openide.util.Lookup;

/**
 *
 * @author Libra
 */
public class GephiTK {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {

        GephiTK gtk = new GephiTK();
        for (String m : new String[]{"s1", "s2", "non", "any"}) {
            double[][] congru = new double[116][];
            for (int i = 1; i < 117; i++) {
                String edgePath = String.format("K:\\code\\jpsth\\gephidata\\%s_edge_%03d.csv", m, i);
                String nodePath = String.format("K:\\code\\jpsth\\gephidata\\%s_node_%03d.csv", m, i);
                congru[i-1] = gtk.processFile(edgePath, nodePath);
            }

            try (FileWriter writer = new FileWriter(String.format("K:\\code\\jpsth\\gephidata\\%s_gephi_graph_sums.csv", m))) {
                for (double[] onesess : congru) {
                    writer.append(Arrays.toString(onesess).replaceAll("[\\[\\]]", ""));
                    writer.append("\n");
                }
            } catch (IOException ioe) {
                System.out.println(ioe.toString());
            }
        }
    }

    public double[] processFile(String edgePath, String nodePath) {
        double[] rtn = new double[]{-1, -1, -1, -1, -1, -1, -1};

        ProjectController pc = Lookup.getDefault().lookup(ProjectController.class);
        pc.newProject();
        Workspace workspace = pc.getCurrentWorkspace();
        ImportController importController = Lookup.getDefault().lookup(ImportController.class);
        GraphModel graphModel = Lookup.getDefault().lookup(GraphController.class).getGraphModel();
        Container nodeContainer, edgeContainer;
        try {
            File edgeFile = new File(edgePath);
            File nodeFile = new File(nodePath);
            nodeContainer = importController.importFile(nodeFile);
            nodeContainer.getLoader().setEdgeDefault(EdgeDirectionDefault.DIRECTED);   //Force DIRECTED
            nodeContainer.getLoader().setAllowAutoNode(true);  //create missing nodes
            nodeContainer.getLoader().setEdgesMergeStrategy(EdgeMergeStrategy.LAST);
            nodeContainer.getLoader().setAutoScale(true);

            edgeContainer = importController.importFile(edgeFile);
            edgeContainer.getLoader().setEdgeDefault(EdgeDirectionDefault.DIRECTED);   //Force DIRECTED
            edgeContainer.getLoader().setAllowAutoNode(true);  //create missing nodes
            edgeContainer.getLoader().setEdgesMergeStrategy(EdgeMergeStrategy.LAST);
            edgeContainer.getLoader().setAutoScale(true);
            importController.process(nodeContainer, new DefaultProcessor(), workspace);
            importController.process(edgeContainer, new AppendProcessor(), workspace);

            DirectedGraph graph = graphModel.getDirectedGraph();
            rtn[0] = graph.getNodeCount();
            rtn[1] = graph.getEdgeCount();

            ConnectedComponents cmpo = new ConnectedComponents();
            cmpo.setDirected(true);
            cmpo.execute(graphModel);
            rtn[2] = cmpo.getConnectedComponentsCount();

            ClusteringCoefficient cc = new ClusteringCoefficient();
            cc.setDirected(true);
            cc.execute(graph);
            rtn[3] = cc.getAverageClusteringCoefficient();

            GraphDensity gdens = new GraphDensity();
            rtn[4] = gdens.calculateDensity(graph, true);

            Degree dg = new Degree();
            dg.execute(graph);
            rtn[5] = dg.getAverageDegree();

            GraphDistance gdist = new GraphDistance();
            gdist.execute(graph);
            rtn[6] = gdist.getPathLength();

        } catch (FileNotFoundException oe) {
            System.out.println("Error Processing File");
            System.out.println(edgePath);
            System.out.println(nodePath);
        }
        return rtn;
    }
}

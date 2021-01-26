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
import org.gephi.graph.GraphControllerImpl;
import org.gephi.graph.api.*;
import org.gephi.io.importer.api.Container;
import org.gephi.io.importer.api.ContainerLoader;
import org.gephi.io.importer.api.EdgeMergeStrategy;
import org.gephi.io.importer.api.ImportController;
import org.gephi.io.importer.impl.ImportControllerImpl;
import org.gephi.project.api.ProjectController;
import org.gephi.project.api.Workspace;
import org.gephi.io.importer.plugin.file.spreadsheet.ImporterSpreadsheetCSV;
import org.gephi.io.processor.plugin.DefaultProcessor;
import org.gephi.project.impl.ProjectControllerImpl;
import org.gephi.statistics.plugin.*;

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
        double[][] congru = new double[115][];
        for (int i = 1; i < 115; i++) {
            String edgePath = String.format("K:\\code\\jpsth\\congru_conn_sc_inter_%03d.csv", i);
            String nodePath = String.format("K:\\code\\jpsth\\congru_conn_sc_node_coord_%03d.csv", i);
//            congru[i] = gtk.processFile(edgePath, nodePath);
            congru[i] = gtk.processFile(edgePath, nodePath,
                    edgePath.replace("congru", "incongru"),
                    nodePath.replace("congru", "incongru"));
//            System.out.println(Arrays.toString(k));
        }

        try (FileWriter writer = new FileWriter("K:\\code\\jpsth\\comb_conn_inter_sums.csv")) {
            for (int j = 1; j < congru.length; j++) {
                writer.append(Arrays.toString(congru[j]));
                writer.append("\n");
            }
        } catch (IOException ioe) {
            System.out.println(ioe.toString());
        }

    }

    public double[] processFile(String edgePath, String nodePath, String... addPath) {
        double[] rtn = new double[]{-1, -1, -1, -1, -1, -1, -1};
        try {
            //Init a project - and therefore a workspace
            ProjectController pc = new ProjectControllerImpl();
            pc.newProject();
            Workspace workspace = pc.getCurrentWorkspace();

            //Get models and controllers for this new workspace - will be useful later
            GraphModel graphModel = new GraphControllerImpl().getGraphModel(workspace);
            ImportController importController = new ImportControllerImpl();

            File edgeFile = new File(edgePath);
            File nodeFile = new File(nodePath);
            ImporterSpreadsheetCSV importer = new ImporterSpreadsheetCSV();
            importer.setFile(edgeFile);

            Container edgeContainer = importController.importFile(
                    edgeFile, importer
            );
            ContainerLoader cl = edgeContainer.getLoader();
            cl.setEdgesMergeStrategy(EdgeMergeStrategy.LAST);

            importController.process(edgeContainer, new DefaultProcessor(), workspace);
            importer.setFile(nodeFile);
            Container nodeContainer = importController.importFile(
                    nodeFile, importer
            );
            importController.process(nodeContainer, new DefaultProcessor(), workspace);

            if (addPath.length > 0) {
                edgeFile = new File(addPath[0]);
                nodeFile = new File(addPath[1]);
                importer.setFile(edgeFile);

                edgeContainer = importController.importFile(
                        edgeFile, importer
                );
                cl = edgeContainer.getLoader();
                cl.setEdgesMergeStrategy(EdgeMergeStrategy.LAST);

                importController.process(edgeContainer, new DefaultProcessor(), workspace);
                importer.setFile(nodeFile);
                nodeContainer = importController.importFile(
                        nodeFile, importer
                );
                importController.process(nodeContainer, new DefaultProcessor(), workspace);

            }

            DirectedGraph graph = graphModel.getDirectedGraph();
//            System.out.println("Nodes: " + graph.getNodeCount());
//            System.out.println("Edges: " + graph.getEdgeCount());
            rtn[0] = graph.getNodeCount();
            rtn[1] = graph.getEdgeCount();

            ConnectedComponents cmpo = new ConnectedComponents();
            cmpo.setDirected(true);
            cmpo.execute(graphModel);
//            System.out.println("Components");
//            System.out.println(cmpo.getConnectedComponentsCount());
            rtn[2] = cmpo.getConnectedComponentsCount();

            ClusteringCoefficient cc = new ClusteringCoefficient();
            cc.setDirected(true);
            cc.execute(graph);
//            System.out.println("CC");
//            System.out.println(cc.getAverageClusteringCoefficient());
            rtn[3] = cc.getAverageClusteringCoefficient();

            GraphDensity gdens = new GraphDensity();
//            System.out.println("Graph Density");
//            System.out.println(gdens.calculateDensity(graph, true));
            rtn[4] = gdens.calculateDensity(graph, true);

            Degree dg = new Degree();
            dg.execute(graph);
//            System.out.println("Avg Degree");
//            System.out.println(dg.getAverageDegree());
            rtn[5] = dg.getAverageDegree();

            GraphDistance gdist = new GraphDistance();
            gdist.execute(graph);
//            System.out.println("Path length");
//            System.out.println(gdist.getPathLength());
            rtn[6] = gdist.getPathLength();
        } catch (FileNotFoundException oe) {
            System.out.println("Error Processing File");
            System.out.println(edgePath);
            System.out.println(nodePath);
        }
        return rtn;

    }

}

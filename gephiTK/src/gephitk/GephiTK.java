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
import java.util.HashSet;
import java.util.LinkedList;
import java.util.Set;
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
import org.gephi.filters.api.FilterController;
import org.gephi.filters.api.Query;
import org.gephi.filters.plugin.attribute.AttributeEqualBuilder;
import org.gephi.filters.plugin.attribute.AttributeEqualBuilder.EqualNumberFilter;

public class GephiTK {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        process_stats();
        process_shufs();
    }

    public static void process_stats() {
        GephiTK gtk = new GephiTK();
        for (String r : new String[]{"within", "cross", "all"}) {
            for (String m : new String[]{"s1", "s2", "incong", "non", "any"}) {
                double[][] congru = new double[116][];
                LinkedList<double[]> cmpPath = new LinkedList<>();
                for (int i = 1; i < 117; i++) {
                    String edgePath = String.format("K:\\code\\jpsth\\gephidata\\%s_edge_%03d_%s.csv", m, i, r);
                    String nodePath = String.format("K:\\code\\jpsth\\gephidata\\%s_node_%03d_%s.csv", m, i, r);
                    LinkedList rtn = gtk.processFile(edgePath, nodePath);
                    congru[i - 1] = (double[]) rtn.get(0);
                    cmpPath.addAll((LinkedList<double[]>) rtn.get(1));
                }

                try (FileWriter writer = new FileWriter(String.format("K:\\code\\jpsth\\gephidata\\%s_%s_gephi_graph_sums.csv", m, r))) {
                    for (double[] onesess : congru) {
                        writer.append(Arrays.toString(onesess).replaceAll("[\\[\\]]", ""));
                        writer.append("\n");
                    }
                } catch (IOException ioe) {
                    System.out.println(ioe.toString());
                }
                try (FileWriter writer = new FileWriter(String.format("K:\\code\\jpsth\\gephidata\\%s_%s_gephi_component_path.csv", m, r))) {
                    for (double[] onesess : cmpPath) {
                        writer.append(Arrays.toString(onesess).replaceAll("[\\[\\]]", ""));
                        writer.append("\n");
                    }
                } catch (IOException ioe) {
                    System.out.println(ioe.toString());
                }

            }
        }
    }

    public static void process_shufs() {
        GephiTK gtk = new GephiTK();
        for (String r : new String[]{"within", "cross", "all"}) {
            for (String m : new String[]{"shuf1_non", "shuf2_non"}) {
                double[][] congru = new double[116][];
                LinkedList<double[]> cmpPath = new LinkedList<>();
                for (int rpt = 0; rpt < 100; rpt++) {
                    for (int i = 1; i < 117; i++) {
                        String edgePath = String.format("K:\\code\\jpsth\\gephidata\\%s_edge_%03d_%s_%d.csv", m, i, r, rpt + 1);
                        String nodePath = String.format("K:\\code\\jpsth\\gephidata\\%s_node_%03d_%s_%d.csv", m, i, r, rpt + 1);
                        if (new File(edgePath).isFile() && new File(nodePath).isFile()) {
                            LinkedList rtn = gtk.processFile(edgePath, nodePath);
                            congru[i - 1] = (double[]) rtn.get(0);
                            cmpPath.addAll((LinkedList<double[]>) rtn.get(1));
                        }
                    }
                    try (FileWriter writer = new FileWriter(String.format("K:\\code\\jpsth\\gephidata\\%s_%s_%d_gephi_graph_sums.csv", m, r, rpt + 1))) {
                        for (double[] onesess : congru) {
                            writer.append(Arrays.toString(onesess).replaceAll("[\\[\\]]", ""));
                            writer.append("\n");
                        }
                    } catch (IOException ioe) {
                        System.out.println(ioe.toString());
                    }

                    try (FileWriter writer = new FileWriter(String.format("K:\\code\\jpsth\\gephidata\\%s_%s_%d_gephi_component_path.csv", m, r, rpt + 1))) {
                        for (double[] onesess : cmpPath) {
                            writer.append(Arrays.toString(onesess).replaceAll("[\\[\\]]", ""));
                            writer.append("\n");
                        }
                    } catch (IOException ioe) {
                        System.out.println(ioe.toString());
                    }
                }
            }
        }
    }

    public LinkedList processFile(String edgePath, String nodePath) {
        LinkedList rtn = new LinkedList();
        double[] stats = new double[]{-1, -1, -1, -1, -1, -1, -1, -1};

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
            nodeContainer.getLoader().setEdgeDefault(EdgeDirectionDefault.UNDIRECTED);   //Force UNDIRECTED
            nodeContainer.getLoader().setAllowAutoNode(true);  //create missing nodes
            nodeContainer.getLoader().setEdgesMergeStrategy(EdgeMergeStrategy.LAST);
            nodeContainer.getLoader().setAutoScale(true);

            edgeContainer = importController.importFile(edgeFile);
            edgeContainer.getLoader().setEdgeDefault(EdgeDirectionDefault.UNDIRECTED);   //Force UNDIRECTED
            edgeContainer.getLoader().setAllowAutoNode(true);  //create missing nodes
            edgeContainer.getLoader().setEdgesMergeStrategy(EdgeMergeStrategy.LAST);
            edgeContainer.getLoader().setAutoScale(true);
            importController.process(nodeContainer, new DefaultProcessor(), workspace);
            importController.process(edgeContainer, new AppendProcessor(), workspace);

            UndirectedGraph graph = graphModel.getUndirectedGraph();
            stats[0] = graph.getNodeCount();
            stats[1] = graph.getEdgeCount();

            ConnectedComponents cmpo = new ConnectedComponents();
            cmpo.setDirected(false);
            cmpo.execute(graphModel);
            stats[2] = cmpo.getConnectedComponentsCount();

            ClusteringCoefficient cc = new ClusteringCoefficient();
            cc.setDirected(false);
            cc.execute(graph);
            stats[3] = cc.getAverageClusteringCoefficient();

            GraphDensity gdens = new GraphDensity();
            gdens.setDirected(false);
            stats[4] = gdens.calculateDensity(graph, true);

            Degree dg = new Degree();
            dg.execute(graph);
            stats[5] = dg.getAverageDegree();

            GraphDistance gdist = new GraphDistance();
            gdist.setDirected(false);
            gdist.execute(graph);
            stats[6] = gdist.getPathLength();
            stats[7] = gdist.getDiameter();

            LinkedList<double[]> cmpPath = pathInSubgraph(graphModel);
            rtn.add(stats);
            rtn.add(cmpPath);
        } catch (FileNotFoundException oe) {
            System.out.println("Error Processing File");
            System.out.println(edgePath);
            System.out.println(nodePath);
        }
        return rtn;
    }

    public static LinkedList<double[]> pathInSubgraph(GraphModel graphModel) {
        LinkedList<double[]> rtn = new LinkedList<>();
        Column modColumn = graphModel.getNodeTable().getColumn("componentnumber");
        Set<Number> components = new HashSet<>();
        for (Node n : graphModel.getGraphVisible().getNodes()) {
            components.add((Number) n.getAttribute(modColumn));
        }
        FilterController filterController = Lookup.getDefault().lookup(FilterController.class);
        EqualNumberFilter enf = new AttributeEqualBuilder.EqualNumberFilter.Node(modColumn);
        for (Number c : components) {
            enf.setMatch(c);
            Query query = filterController.createQuery(enf);
            GraphView view = filterController.filter(query);
            graphModel.setVisibleView(view);
            Graph subgraph = graphModel.getUndirectedGraphVisible();
            int nc = subgraph.getNodeCount();
            if (nc > 2) {
                GraphDistance gdist = new GraphDistance();
                gdist.setDirected(false);
                gdist.execute(subgraph);
                double pl = gdist.getPathLength();
                double di = gdist.getDiameter();
                gdist.execute(subgraph);
                double radi = gdist.getRadius();
                rtn.add(new double[]{nc, pl, di, radi});
            }
        }
        return rtn;
    }
}

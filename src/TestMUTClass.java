import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.*;

public class TestMUTClass
{
    public static void main(String[] args) throws Exception
    {
        String mutationsFile="Data/BRCA_snp_gene_matrix.txt";
        //String listDriverGenes="BRCA1,BRCA2";
        String listDriverGenes="ABCE1";
        //String listDriverGenes=null;
        //String weightsFile="Data/CDS_Lengths.txt";
        String weightsFile=null;
        //int[] maxNumSolGenesPos={5,10,15,20,25,30};
        int[] maxNumSolGenesPos={5,15,30,50};
        //int[] maxNumSolGenesNeg={5,10,15,20,25,30};
        int[] maxNumSolGenesNeg={5,15,30,50};
        int numIterCV=1;
        int numFolds=5;
        double alpha=0.8;
        String outputFolder="Results/CrossValidation_BRCA";

        //System.out.println("Reading mutation matrix...");
        FileManager fm=new FileManager();
        //System.out.println(setPositives.size());
        //BufferedWriter bw=new BufferedWriter(new FileWriter(outputFolder+".txt"));

        //Read mutation data
        MutationMatrix mutMatrix;
        if(listDriverGenes==null)
            mutMatrix=fm.readMutationDataWithClasses(mutationsFile);
        else
        {
            String[] split=listDriverGenes.split(",");
            HashSet<String> setRefGenes = new HashSet<>(Arrays.asList(split));
            mutMatrix=fm.readMutationDataWithoutClasses(mutationsFile,setRefGenes);
        }
        HashMap<String,HashSet<String>> mutationData=mutMatrix.getMutationData();
        Hashtable<String,String> mapSampleClasses=mutMatrix.getMapSampleClasses();
        //System.out.println(mapSampleClasses.size());
        HashSet<String> setPositives=mutMatrix.getPositiveSet();
        HashSet<String> setNegatives=mutMatrix.getNegativeSet();
        Vector<String> listPositives=new Vector<>(setPositives);
        Vector<String> listNegatives=new Vector<>(setNegatives);
        //System.out.println(setPositives.size());
        //System.out.println(setNegatives.size());

        //Read gene weights (if specified)
        HashMap<String,Double> mapWeights;
        if(weightsFile==null)
        {
            mapWeights = null;
            alpha=1.0;
        }
        else
            mapWeights=fm.readWeightsFile(weightsFile);

        //Get maximum positive and negative panel size
        int maxPosSize=Arrays.stream(maxNumSolGenesPos).max().getAsInt();
        int maxNegSize=Arrays.stream(maxNumSolGenesNeg).max().getAsInt();

        //Initialize averages over cross-validation iterations
        HashMap<String,Double> avgAccuracy = new HashMap<>();
        HashMap<String,Double> avgPrecision = new HashMap<>();
        HashMap<String,Double> avgRecall = new HashMap<>();
        HashMap<String,Double> avgFpr = new HashMap<>();
        HashMap<String,Double> avgFnr = new HashMap<>();
        HashMap<String,Double> avgSpecificity = new HashMap<>();
        HashMap<String,Double> avgBalancedAccuracy = new HashMap<>();
        HashMap<String,Double> avgUnclassified = new HashMap<>();
        HashMap<String,Double> avgTruePositives = new HashMap<>();
        HashMap<String,Double> avgFalsePositives = new HashMap<>();
        HashMap<String,Double> avgTrueNegatives = new HashMap<>();
        HashMap<String,Double> avgFalseNegatives = new HashMap<>();
        HashMap<String,Double> avgF1 = new HashMap<>();
        for(int m=0;m<maxNumSolGenesPos.length;m++)
        {
            for (int n = 0; n < maxNumSolGenesNeg.length; n++)
            {
                avgAccuracy.put(m+"-"+n,0.0);
                avgPrecision.put(m+"-"+n,0.0);
                avgRecall.put(m+"-"+n,0.0);
                avgFpr.put(m+"-"+n,0.0);
                avgFnr.put(m+"-"+n,0.0);
                avgSpecificity.put(m+"-"+n,0.0);
                avgBalancedAccuracy.put(m+"-"+n,0.0);
                avgUnclassified.put(m+"-"+n,0.0);
                avgTruePositives.put(m+"-"+n,0.0);
                avgFalsePositives.put(m+"-"+n,0.0);
                avgTrueNegatives.put(m+"-"+n,0.0);
                avgFalseNegatives.put(m+"-"+n,0.0);
                avgF1.put(m+"-"+n,0.0);
            }
        }

        for(int i=0;i<numIterCV;i++)
        {
            //System.out.println("CV iteration: "+(i+1)+"/"+numIterCV);

            //Create folds
            Collections.shuffle(listPositives);
            Collections.shuffle(listNegatives);
            Vector<List<String>> posFolds = new Vector<>();
            Vector<List<String>> negFolds = new Vector<>();
            int foldPosSize = (int) (((double) listPositives.size()) / numFolds);
            int foldNegSize = (int) (((double) listNegatives.size()) / numFolds);
            int startPos = 0;
            int startNeg = 0;
            for (int j = 0; j < numFolds - 1; j++)
            {
                posFolds.add(listPositives.subList(startPos, startPos + foldPosSize));
                negFolds.add(listNegatives.subList(startNeg, startNeg + foldNegSize));
                startPos = startPos + foldPosSize;
                startNeg = startNeg + foldNegSize;
            }
            posFolds.add(listPositives.subList(startPos, listPositives.size()));
            negFolds.add(listNegatives.subList(startNeg, listNegatives.size()));

            //Initialize total TP, TN, FP, FN and number of unclassified samples over all folds
            HashMap<String,Double> truePositives = new HashMap<>();
            HashMap<String,Double> trueNegatives = new HashMap<>();
            HashMap<String,Double> falsePositives = new HashMap<>();
            HashMap<String,Double> falseNegatives = new HashMap<>();
            HashMap<String,Double> unclassified = new HashMap<>();
            for(int m=0;m<maxNumSolGenesPos.length;m++)
            {
                for (int n = 0; n < maxNumSolGenesNeg.length; n++)
                {
                    truePositives.put(m+"-"+n,0.0);
                    trueNegatives.put(m+"-"+n,0.0);
                    falsePositives.put(m+"-"+n,0.0);
                    falseNegatives.put(m+"-"+n,0.0);
                    unclassified.put(m+"-"+n,0.0);
                }
            }

            for (int j = 0; j < numFolds; j++)
            {

                //Create training set
                HashSet<String> trainingSetSamples = new HashSet<>();
                for (int k = 0; k < numFolds; k++)
                {
                    if (k != j)
                    {
                        trainingSetSamples.addAll(posFolds.get(k));
                        trainingSetSamples.addAll(negFolds.get(k));
                    }
                }
                //System.out.println("FOLD "+(j+1));
                //System.out.println("Training size: "+trainingSetSamples.size());
                HashMap<String, HashSet<String>> trainingMutationData = Utility.extractSubMutData(mutationData, trainingSetSamples);
                Hashtable<String, String> trainingSampleClasses = Utility.extractSubSampleClasses(mapSampleClasses, trainingSetSamples);
                //System.out.println("SAMPLES TRAINING: "+trainingSampleClasses.size());

                //Create test set
                HashSet<String> testSetSamples = new HashSet<>();
                testSetSamples.addAll(posFolds.get(j));
                testSetSamples.addAll(negFolds.get(j));
                //System.out.println("FOLD "+(j+1));
                //System.out.println("Test size: "+testSetSamples.size());
                HashMap<String, HashSet<String>> testMutationData = Utility.extractSubMutData(mutationData, testSetSamples);
                Hashtable<String, String> testSampleClasses = Utility.extractSubSampleClasses(mapSampleClasses, testSetSamples);
                //System.out.println("SAMPLES TEST: "+testSampleClasses.size());

                //Run MUTClass with maximum positive and negative panel size
                DMGSFinderApp mc = new DMGSFinderApp(trainingMutationData, trainingSampleClasses, mapWeights, true, maxPosSize, alpha);
                Vector<String> panelPosGenes = mc.runAlgorithm();
                mc = new DMGSFinderApp(trainingMutationData, trainingSampleClasses, mapWeights, false, maxNegSize, alpha);
                Vector<String> panelNegGenes = mc.runAlgorithm();
                if (panelPosGenes == null || panelNegGenes == null)
                    continue;
                //System.out.println("Panel size for positives: "+panelPosGenes.size());
                //System.out.println("Panel size for negatives: "+panelNegGenes.size());

                for (int m = 0; m < maxNumSolGenesPos.length; m++)
                {
                    for (int n = 0; n < maxNumSolGenesNeg.length; n++)
                    {
                        //Extract sub-panels of corresponding sizes
                        List<String> subPanelPosGenes = panelPosGenes.subList(0, maxNumSolGenesPos[m]);
                        List<String> subPanelNegGenes = panelNegGenes.subList(0, maxNumSolGenesNeg[n]);
                        //System.out.println(subPanelPosGenes);

                        //Make classification
                        RuleClassifier rc = new RuleClassifier(subPanelPosGenes, subPanelNegGenes);
                        Hashtable<String, String> subPredictedClasses = rc.predictClasses2(testMutationData, testSampleClasses);
                        //System.out.println(subPredictedClasses);
                        Hashtable<String, Double> accuracyResults = rc.evaluateResults2(subPredictedClasses, testSampleClasses);
                        //System.out.println(accuracyResults);

                        //Update accuracy results
                        truePositives.put(m + "-" + n, truePositives.get(m + "-" + n) + accuracyResults.get("TP"));
                        //System.out.println(truePositives.get(m+"-"+n));
                        trueNegatives.put(m + "-" + n, trueNegatives.get(m + "-" + n) + accuracyResults.get("TN"));
                        falsePositives.put(m + "-" + n, falsePositives.get(m + "-" + n) + accuracyResults.get("FP"));
                        falseNegatives.put(m + "-" + n, falseNegatives.get(m + "-" + n) + accuracyResults.get("FN"));
                        unclassified.put(m + "-" + n, unclassified.get(m + "-" + n) + accuracyResults.get("Unclassified"));
                    }
                }
            }

            //Compute accuracy, precision, recall and other statistics for i-th iteration of CV
            for (int m = 0; m < maxNumSolGenesPos.length; m++)
            {
                for (int n = 0; n < maxNumSolGenesNeg.length; n++)
                {
                    double precision=(truePositives.get(m+"-"+n)/(truePositives.get(m+"-"+n)+falsePositives.get(m+"-"+n)))*100;
                    double recall=(truePositives.get(m+"-"+n)/(truePositives.get(m+"-"+n)+falseNegatives.get(m+"-"+n)))*100;
                    double specificity=(trueNegatives.get(m+"-"+n)/(trueNegatives.get(m+"-"+n)+falsePositives.get(m+"-"+n)))*100;
                    double fpr=(falsePositives.get(m+"-"+n)/(falsePositives.get(m+"-"+n)+trueNegatives.get(m+"-"+n)))*100;
                    double fnr=(falseNegatives.get(m+"-"+n)/(falseNegatives.get(m+"-"+n)+truePositives.get(m+"-"+n)))*100;
                    double accuracy=((truePositives.get(m+"-"+n)+trueNegatives.get(m+"-"+n))/(truePositives.get(m+"-"+n)
                            +trueNegatives.get(m+"-"+n)+falsePositives.get(m+"-"+n)+falseNegatives.get(m+"-"+n)))*100;
                    double balancedAccuracy=(recall+specificity)/2;
                    double f1=2*(precision*recall/(precision+recall));
                    /*System.out.println("Accuracy: "+accuracy);
                    System.out.println("Precision: "+precision);
                    System.out.println("Sensitivity/Recall/TPR: "+recall);
                    System.out.println("Specificity: "+specificity);
                    System.out.println("FPR: "+fpr);
                    System.out.println("Balanced accuracy: "+balancedAccuracy);*/

                    //Update average statistics
                    avgAccuracy.put(m+"-"+n, avgAccuracy.get(m+"-"+n)+accuracy);
                    avgPrecision.put(m+"-"+n, avgPrecision.get(m+"-"+n)+precision);
                    avgRecall.put(m+"-"+n, avgRecall.get(m+"-"+n)+recall);
                    avgFpr.put(m+"-"+n, avgFpr.get(m+"-"+n)+fpr);
                    avgFnr.put(m+"-"+n, avgFnr.get(m+"-"+n)+fnr);
                    avgSpecificity.put(m+"-"+n, avgSpecificity.get(m+"-"+n)+specificity);
                    avgBalancedAccuracy.put(m+"-"+n, avgBalancedAccuracy.get(m+"-"+n)+balancedAccuracy);
                    avgUnclassified.put(m+"-"+n, avgUnclassified.get(m+"-"+n)+unclassified.get(m+"-"+n)/mapSampleClasses.size()*100);
                    avgTruePositives.put(m+"-"+n, avgTruePositives.get(m+"-"+n)+truePositives.get(m+"-"+n));
                    avgFalsePositives.put(m+"-"+n, avgFalsePositives.get(m+"-"+n)+falsePositives.get(m+"-"+n));
                    avgTrueNegatives.put(m+"-"+n, avgTrueNegatives.get(m+"-"+n)+trueNegatives.get(m+"-"+n));
                    avgFalseNegatives.put(m+"-"+n, avgFalseNegatives.get(m+"-"+n)+falseNegatives.get(m+"-"+n));
                    avgF1.put(m+"-"+n, avgF1.get(m+"-"+n)+f1);

                }
            }
        }

        //Normalize over the number of CVs and get the final averages
        for(int m=0;m<maxNumSolGenesPos.length;m++)
        {
            for (int n = 0; n < maxNumSolGenesNeg.length; n++)
            {
                avgAccuracy.put(m+"-"+n, avgAccuracy.get(m+"-"+n)/numIterCV);
                avgPrecision.put(m+"-"+n, avgPrecision.get(m+"-"+n)/numIterCV);
                avgRecall.put(m+"-"+n, avgRecall.get(m+"-"+n)/numIterCV);
                avgFpr.put(m+"-"+n, avgFpr.get(m+"-"+n)/numIterCV);
                avgFnr.put(m+"-"+n, avgFnr.get(m+"-"+n)/numIterCV);
                avgSpecificity.put(m+"-"+n, avgSpecificity.get(m+"-"+n)/numIterCV);
                avgBalancedAccuracy.put(m+"-"+n, avgBalancedAccuracy.get(m+"-"+n)/numIterCV);
                avgUnclassified.put(m+"-"+n, avgUnclassified.get(m+"-"+n)/numIterCV);
                avgTruePositives.put(m+"-"+n, avgTruePositives.get(m+"-"+n)/numIterCV);
                avgTrueNegatives.put(m+"-"+n, avgTrueNegatives.get(m+"-"+n)/numIterCV);
                avgFalsePositives.put(m+"-"+n, avgFalsePositives.get(m+"-"+n)/numIterCV);
                avgFalseNegatives.put(m+"-"+n, avgFalseNegatives.get(m+"-"+n)/numIterCV);
                avgF1.put(m+"-"+n, avgF1.get(m+"-"+n)/numIterCV);

                /*bw.write("------------------------\n");
                bw.write("POSITIVE RULE SIZE: "+maxNumSolGenesPos[m]+"\n");
                bw.write("NEGATIVE RULE SIZE: "+maxNumSolGenesNeg[n]+"\n");
                bw.write("AVERAGE PRECISION: "+avgPrecision+"\n");
                bw.write("AVERAGE SENSITIVITY/RECALL/TPR: "+avgRecall+"\n");
                bw.write("AVERAGE FPR: "+avgFpr+"\n");
                bw.write("AVERAGE FNR: "+avgFnr+"\n");
                bw.write("AVERAGE SPECIFICITY: "+avgSpecificity+"\n");
                bw.write("AVERAGE ACCURACY: "+avgAccuracy+"\n");
                bw.write("AVERAGE BALANCED ACCURACY: "+avgBalancedAccuracy+"\n");
                bw.write("AVERAGE PERCENTAGE OF UNCLASSIFIED: "+avgUnclassified+"\n");
                bw.write("AVERAGE NUMBER OF TRUE POSITIVES: "+avgTruePositives+"\n");
                bw.write("AVERAGE NUMBER OF FALSE POSITIVES: "+avgFalsePositives+"\n");
                bw.write("AVERAGE NUMBER OF TRUE NEGATIVES: "+avgTrueNegatives+"\n");
                bw.write("AVERAGE NUMBER OF FALSE NEGATIVES: "+avgFalseNegatives+"\n");*/

                System.out.println("------------------------");
                System.out.println("POSITIVE RULE SIZE: "+maxNumSolGenesPos[m]);
                System.out.println("NEGATIVE RULE SIZE: "+maxNumSolGenesNeg[n]);
                System.out.println("AVERAGE PRECISION: "+avgPrecision.get(m+"-"+n));
                System.out.println("AVERAGE SENSITIVITY/RECALL/TPR: "+avgRecall.get(m+"-"+n));
                System.out.println("AVERAGE FPR: "+avgFpr.get(m+"-"+n));
                System.out.println("AVERAGE FNR: "+avgFnr.get(m+"-"+n));
                System.out.println("AVERAGE SPECIFICITY: "+avgSpecificity.get(m+"-"+n));
                System.out.println("AVERAGE ACCURACY: "+avgAccuracy.get(m+"-"+n));
                //System.out.println("AVERAGE BALANCED ACCURACY: "+avgBalancedAccuracy.get(m+"-"+n));
                System.out.println("AVERAGE F1 SCORE: "+avgF1.get(m+"-"+n));
                System.out.println("AVERAGE PERCENTAGE OF UNCLASSIFIED: "+avgUnclassified.get(m+"-"+n));
                System.out.println("AVERAGE NUMBER OF TRUE POSITIVES: "+avgTruePositives.get(m+"-"+n));
                System.out.println("AVERAGE NUMBER OF FALSE POSITIVES: "+avgFalsePositives.get(m+"-"+n));
                System.out.println("AVERAGE NUMBER OF TRUE NEGATIVES: "+avgTrueNegatives.get(m+"-"+n));
                System.out.println("AVERAGE NUMBER OF FALSE NEGATIVES: "+avgFalseNegatives.get(m+"-"+n));

            }
        }
        //bw.close();

    }

}

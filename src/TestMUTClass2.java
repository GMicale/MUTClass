import java.io.*;
import java.util.*;

public class TestMUTClass2
{
    public static void main(String[] args) throws Exception
    {
        String[] tumors={"BRCA","COAD","LUAD"};
        int[] maxNumSolGenesPos={5,15,30,50};
        int[] maxNumSolGenesNeg={5,15,30,50};
        int numIterCV=10;
        int numFolds=5;

        //Read mutation data
        //System.out.println("Reading mutation matrix...");
        FileManager fm=new FileManager();
        //System.out.println(setPositives.size());
        for(int a=0;a<tumors.length;a++)
        {
            String tumor=tumors[a];
            File dir=new File("C:/Ricerca/MUTClass/Dataset/"+tumor);
            File[] listFiles=dir.listFiles();
            for(int b=0;b<listFiles.length;b++)
            {
                String filename = listFiles[b].getName();
                String inputFile="C:/Ricerca/MUTClass/Dataset/"+tumor+"/"+filename;
                String outputFolder="Results/CrossValidation_"+tumor;
                String[] splitFile=filename.replace(".txt","").split("_");
                String classType=splitFile[4];
                String classVal=splitFile[5];
                System.out.println("Cross validation for "+tumor+" - "+classType+" - "+classVal+"...");
                BufferedWriter bw=new BufferedWriter(new FileWriter(outputFolder+"/"+tumor+"_"+classType+"_"+classVal+".txt"));
                MutationMatrix mutMatrix=fm.readMutationDataWithClasses(inputFile);
                HashMap<String,HashSet<String>> mutationData=mutMatrix.getMutationData();
                Hashtable<String,String> mapSampleClasses=mutMatrix.getMapSampleClasses();
                //System.out.println(mapSampleClasses.size());
                HashSet<String> setPositives=mutMatrix.getPositiveSet();
                HashSet<String> setNegatives=mutMatrix.getNegativeSet();
                Vector<String> listPositives=new Vector<>(setPositives);
                Vector<String> listNegatives=new Vector<>(setNegatives);
                //System.out.println(setPositives.size());
                //System.out.println(setNegatives.size());
                if(setPositives.size()==0 || setNegatives.size()==0)
                    continue;
                for(int m=0;m<maxNumSolGenesPos.length;m++)
                {
                    for(int n=0;n<maxNumSolGenesNeg.length;n++)
                    {
                        int i, j, k;
                        double avgAccuracy=0.0;
                        double avgPrecision=0.0;
                        double avgRecall=0.0;
                        double avgFpr=0.0;
                        double avgFnr=0.0;
                        double avgSpecificity=0.0;
                        double avgBalancedAccuracy=0.0;
                        double avgUnclassified=0.0;
                        double avgTruePositives=0.0;
                        double avgFalsePositives=0.0;
                        double avgTrueNegatives=0.0;
                        double avgFalseNegatives=0.0;

                        for(i=0;i<numIterCV;i++)
                        {
                            //System.out.println("CV iteration: "+(i+1)+"/"+numIterCV);
                            Collections.shuffle(listPositives);
                            Collections.shuffle(listNegatives);

                            //Create folds
                            Vector<List<String>> posFolds=new Vector<>();
                            Vector<List<String>> negFolds=new Vector<>();
                            int foldPosSize=(int)(((double)listPositives.size())/numFolds);
                            int foldNegSize=(int)(((double)listNegatives.size())/numFolds);
                            int startPos=0;
                            int startNeg=0;
                            for(j=0;j<numFolds-1;j++)
                            {
                                posFolds.add(listPositives.subList(startPos,startPos+foldPosSize));
                                negFolds.add(listNegatives.subList(startNeg,startNeg+foldNegSize));
                                startPos=startPos+foldPosSize;
                                startNeg=startNeg+foldNegSize;
                            }
                            posFolds.add(listPositives.subList(startPos,listPositives.size()));
                            negFolds.add(listNegatives.subList(startNeg,listNegatives.size()));

                            double truePositives=0.0;
                            double trueNegatives=0.0;
                            double falsePositives=0.0;
                            double falseNegatives=0.0;
                            double unclassified=0.0;

                            for(j=0;j<numFolds;j++)
                            {
                                //Create positive and negative training and testing sets
                                HashSet<String> trainingSetSamples=new HashSet<>();
                                for(k=0;k<numFolds;k++)
                                {
                                    if(k!=j)
                                    {
                                        trainingSetSamples.addAll(posFolds.get(k));
                                        trainingSetSamples.addAll(negFolds.get(k));
                                    }
                                }
                                //System.out.println("FOLD "+(j+1));
                                //System.out.println("Training size: "+trainingSetSamples.size());
                                HashMap<String,HashSet<String>> trainingMutationData=Utility.extractSubMutData(mutationData,trainingSetSamples);
                                Hashtable<String,String> trainingSampleClasses=Utility.extractSubSampleClasses(mapSampleClasses,trainingSetSamples);
                                //System.out.println("SAMPLES TRAINING: "+trainingSampleClasses.size());

                                //Run MUTClass to find panels of mutated genes
                                DMGSFinderApp mc=new DMGSFinderApp(trainingMutationData,trainingSampleClasses,null,true,maxNumSolGenesPos[m],1.0);
                                Vector<String> panelPosGenes=mc.runAlgorithm();
                                mc=new DMGSFinderApp(trainingMutationData,trainingSampleClasses,null,false,maxNumSolGenesNeg[n],1.0);
                                Vector<String> panelNegGenes=mc.runAlgorithm();
                                if(panelPosGenes==null || panelNegGenes==null)
                                    continue;
                                //System.out.println("Panel size for positives: "+panelPosGenes.size());
                                //System.out.println("Panel size for negatives: "+panelNegGenes.size());

                                //Extract test set
                                HashSet<String> testSetSamples=new HashSet<>();
                                testSetSamples.addAll(posFolds.get(j));
                                testSetSamples.addAll(negFolds.get(j));
                                //System.out.println("FOLD "+(j+1));
                                //System.out.println("Test size: "+testSetSamples.size());
                                HashMap<String,HashSet<String>> testMutationData=Utility.extractSubMutData(mutationData,testSetSamples);
                                Hashtable<String,String> testSampleClasses=Utility.extractSubSampleClasses(mapSampleClasses,testSetSamples);
                                //System.out.println("SAMPLES TEST: "+testSampleClasses.size());

                                //Make classification
                                RuleClassifier rc=new RuleClassifier(panelPosGenes,panelNegGenes);
                                Hashtable<String,String> subPredictedClasses=rc.predictClasses2(testMutationData,testSampleClasses);
                                Hashtable<String,Double> accuracyResults=rc.evaluateResults2(subPredictedClasses,testSampleClasses);

                                //Update accuracy results
                                truePositives+=accuracyResults.get("TP");
                                trueNegatives+=accuracyResults.get("TN");
                                falsePositives+=accuracyResults.get("FP");
                                falseNegatives+=accuracyResults.get("FN");
                                unclassified+=accuracyResults.get("Unclassified");

                            }
                            //System.out.println();

                            //Compute accuracy, precision and recall
                            /*System.out.println("TP: "+truePositives);
                            System.out.println("FN: "+falseNegatives);
                            System.out.println("TN: "+trueNegatives);
                            System.out.println("FP: "+falsePositives);
                            System.out.println("Unclassified: "+unclassified);*/
                            double precision=(truePositives/(truePositives+falsePositives))*100;
                            double recall=(truePositives/(truePositives+falseNegatives))*100;
                            double specificity=(trueNegatives/(trueNegatives+falsePositives))*100;
                            double fpr=(falsePositives/(falsePositives+trueNegatives))*100;
                            double fnr=(falseNegatives/(falseNegatives+truePositives))*100;
                            double accuracy=((truePositives+trueNegatives)/(truePositives+trueNegatives+falsePositives+falseNegatives))*100;
                            double balancedAccuracy=(recall+specificity)/2;
                            /*System.out.println("Accuracy: "+accuracy);
                            System.out.println("Precision: "+precision);
                            System.out.println("Sensitivity/Recall/TPR: "+recall);
                            System.out.println("Specificity: "+specificity);
                            System.out.println("FPR: "+fpr);
                            System.out.println("Balanced accuracy: "+balancedAccuracy);*/

                            //Update average statistics
                            avgAccuracy+=accuracy;
                            avgPrecision+=precision;
                            avgRecall+=recall;
                            avgFpr+=fpr;
                            avgFnr+=fnr;
                            avgSpecificity+=specificity;
                            avgBalancedAccuracy+=balancedAccuracy;
                            avgUnclassified+=(unclassified/mapSampleClasses.size()*100);
                            avgTruePositives+=truePositives;
                            avgFalsePositives+=falsePositives;
                            avgTrueNegatives+=trueNegatives;
                            avgFalseNegatives+=falseNegatives;

                        }

                        avgAccuracy/=numIterCV;
                        avgPrecision/=numIterCV;
                        avgRecall/=numIterCV;
                        avgFpr/=numIterCV;
                        avgFnr/=numIterCV;
                        avgSpecificity/=numIterCV;
                        avgBalancedAccuracy/=numIterCV;
                        avgUnclassified/=numIterCV;
                        avgTruePositives/=numIterCV;
                        avgTrueNegatives/=numIterCV;
                        avgFalsePositives/=numIterCV;
                        avgFalseNegatives/=numIterCV;

                        bw.write("------------------------\n");
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
                        bw.write("AVERAGE NUMBER OF FALSE NEGATIVES: "+avgFalseNegatives+"\n");

                        /*System.out.println("------------------------");
                        System.out.println("POSITIVE RULE SIZE: "+maxNumSolGenesPos[m]);
                        System.out.println("NEGATIVE RULE SIZE: "+maxNumSolGenesNeg[n]);
                        System.out.println("AVERAGE PRECISION: "+avgPrecision);
                        System.out.println("AVERAGE SENSITIVITY/RECALL/TPR: "+avgRecall);
                        System.out.println("AVERAGE FPR: "+avgFpr);
                        System.out.println("AVERAGE FNR: "+avgFnr);
                        System.out.println("AVERAGE SPECIFICITY: "+avgSpecificity);
                        System.out.println("AVERAGE ACCURACY: "+avgAccuracy);
                        System.out.println("AVERAGE BALANCED ACCURACY: "+avgBalancedAccuracy);
                        System.out.println("AVERAGE PERCENTAGE OF UNCLASSIFIED: "+avgUnclassified);
                        System.out.println("AVERAGE NUMBER OF TRUE POSITIVES: "+avgTruePositives);
                        System.out.println("AVERAGE NUMBER OF FALSE POSITIVES: "+avgFalsePositives);
                        System.out.println("AVERAGE NUMBER OF TRUE NEGATIVES: "+avgTrueNegatives);
                        System.out.println("AVERAGE NUMBER OF FALSE NEGATIVES: "+avgFalseNegatives);*/
                    }
                }
                bw.close();
            }
        }
    }
}

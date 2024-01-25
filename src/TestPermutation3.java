import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.*;

public class TestPermutation3
{
    public static void main(String[] args) throws Exception
    {
        String tumor="LUAD";
        String inputFile="Data/"+tumor+"_snp_gene_matrix.txt";
        String refGeneFile="Results/CrossValidation/CV_genes_"+tumor+".txt";
        int posPanelSize=5;
        int negPanelSize=50;
        int numVariants=1000;
        String outputFolder="Results/PermutationTest";

        //Read the list of reference genes
        Vector<String> listRefGenes=new Vector<>();
        BufferedReader br=new BufferedReader(new FileReader(refGeneFile));
        String str;
        while((str=br.readLine())!=null)
            listRefGenes.add(str);
        br.close();
        //System.out.println(listRefGenes);

        //System.out.println("Reading mutation matrix...");
        FileManager fm=new FileManager();

        /*System.out.print("Gene\tTP\tTN\tFP\tFN\tUnclassified\tPrecision\tRecall\tSpecificity\tFPR\tFNR\tAccuracy\tF1Score");
        System.out.print("\tTPRand\tTNRand\tFPRand\tFNRand\tUnclassifiedRand\tPrecisionRand\tRecallRand\tSpecificityRand");
        System.out.println("\tFPRRand\tFNRRand\tAccuracyRand\tF1ScoreRand\tpvalAccuracy\tpvalF1Score");*/

        BufferedWriter bw=new BufferedWriter(new FileWriter(outputFolder+"/"+tumor+"_perm.txt"));
        bw.write("Gene\tTP\tTN\tFP\tFN\tUnclassified\tPrecision\tRecall\tSpecificity\tFPR\tFNR\tAccuracy\tF1Score");
        bw.write("\tTPRand\tTNRand\tFPRand\tFNRand\tUnclassifiedRand\tPrecisionRand\tRecallRand\tSpecificityRand");
        bw.write("\tFPRRand\tFNRRand\tAccuracyRand\tF1ScoreRand\tpvalAccuracy\tpvalF1Score\n");
        bw.close();

        for (String refGene : listRefGenes)
        {
            String[] split = refGene.split(",");
            HashSet<String> setRefGenes = new HashSet<>(Arrays.asList(split));
            System.out.println("Permutation test for gene "+refGene+"...");

            //Read mutation data
            MutationMatrix mutMatrix = fm.readMutationDataWithoutClasses(inputFile, setRefGenes);
            HashMap<String, HashSet<String>> mutationData = mutMatrix.getMutationData();
            Hashtable<String, String> mapSampleClasses = mutMatrix.getMapSampleClasses();
            //System.out.println(mapSampleClasses.size());
            HashSet<String> setPositives = mutMatrix.getPositiveSet();
            HashSet<String> setNegatives = mutMatrix.getNegativeSet();
            Vector<String> listPositives = new Vector<>(setPositives);
            Vector<String> listNegatives = new Vector<>(setNegatives);
            //System.out.println(setPositives.size());
            //System.out.println(setNegatives.size());

            //Create training set using holdout method
            Random r=new Random(1250);
            Collections.shuffle(listPositives,r);
            Collections.shuffle(listNegatives,r);
            int trainingPosSize = (int) (((double) listPositives.size()) * 80/100);
            int trainingNegSize = (int) (((double) listNegatives.size()) * 80/100);
            HashSet<String> trainingSetSamples = new HashSet<>();
            trainingSetSamples.addAll(listPositives.subList(0,trainingPosSize));
            trainingSetSamples.addAll(listNegatives.subList(0,trainingNegSize));
            //System.out.println("Training size: "+trainingSetSamples.size());
            HashMap<String, HashSet<String>> trainingMutationData = Utility.extractSubMutData(mutationData, trainingSetSamples);
            Hashtable<String, String> trainingSampleClasses = Utility.extractSubSampleClasses(mapSampleClasses, trainingSetSamples);

            //Create test set
            HashSet<String> testSetSamples = new HashSet<>();
            testSetSamples.addAll(listPositives.subList(trainingPosSize,listPositives.size()));
            testSetSamples.addAll(listNegatives.subList(trainingNegSize,listNegatives.size()));
            //System.out.println("Test size: "+testSetSamples.size());
            HashMap<String, HashSet<String>> testMutationData = Utility.extractSubMutData(mutationData, testSetSamples);
            Hashtable<String, String> testSampleClasses = Utility.extractSubSampleClasses(mapSampleClasses, testSetSamples);

            //Run MUTClass
            DMGSFinderApp mc = new DMGSFinderApp(trainingMutationData, trainingSampleClasses, null, true, posPanelSize, 1.0);
            Vector<String> panelPosGenes = mc.runAlgorithm();
            mc = new DMGSFinderApp(trainingMutationData, trainingSampleClasses, null, false, negPanelSize, 1.0);
            Vector<String> panelNegGenes = mc.runAlgorithm();

            //Make classification
            RuleClassifier rc = new RuleClassifier(panelPosGenes, panelNegGenes);
            Hashtable<String, String> predictedClasses = rc.predictClasses2(testMutationData, testSampleClasses);
            //System.out.println(predictedClasses);
            Hashtable<String, Double> accuracyResults = rc.evaluateResults2(predictedClasses, testSampleClasses);
            //System.out.println(accuracyResults);

            //Compute accuracy results
            double tp=accuracyResults.get("TP");
            double tn=accuracyResults.get("TN");
            double fp=accuracyResults.get("FP");
            double fn=accuracyResults.get("FN");
            double unclassified=accuracyResults.get("Unclassified");
            double precision=(tp/(tp+fp))*100;
            double recall=(tp/(tp+fn))*100;
            double specificity=(tn/(tn+fp))*100;
            double fpr=(fp/(fp+tn))*100;
            double fnr=(fn/(fn+tp))*100;
            double accuracy=((tp+tn)/(tp+tn+fp+fn))*100;
            double f1score;
            if(precision==0.0 && recall==0.0)
                f1score=0.0;
            else if(Double.isNaN(precision) || Double.isNaN(recall))
                f1score=0.0;
            else
                f1score=2*(precision*recall/(precision+recall));

            //Initialize average accuracy results for random variants
            double avgTpRand=0.0;
            double avgTnRand=0.0;
            double avgFpRand=0.0;
            double avgFnRand=0.0;
            double avgUnclassifiedRand=0.0;
            double avgPrecisionRand=0.0;
            double avgRecallRand=0.0;
            double avgSpecificityRand=0.0;
            double avgFprRand=0.0;
            double avgFnrRand=0.0;
            double avgAccuracyRand=0.0;
            double pvalAccuracy=0.0;
            double avgF1scoreRand=0.0;
            double pvalF1score=0.0;

            for(int i=0;i<numVariants;i++)
            {
                //Permute sample classes to create a random variant of the training set
                Hashtable<String,String> randomTrainingClasses=new Hashtable<>();
                Vector<String> listClasses = new Vector<>(trainingSampleClasses.values());
                Collections.shuffle(listClasses,new Random(i));
                int j=0;
                for(String sample : trainingSampleClasses.keySet())
                {
                    randomTrainingClasses.put(sample,listClasses.get(j));
                    j++;
                }
                //System.out.println(listClasses);

                //Run MUTClass
                DMGSFinderApp mcRand = new DMGSFinderApp(trainingMutationData, randomTrainingClasses, null, true, posPanelSize, 1.0);
                Vector<String> panelPosGenesRand = mcRand.runAlgorithm();
                mcRand = new DMGSFinderApp(trainingMutationData, randomTrainingClasses, null, false, negPanelSize, 1.0);
                Vector<String> panelNegGenesRand = mcRand.runAlgorithm();

                //Make classification
                RuleClassifier rcRand = new RuleClassifier(panelPosGenesRand, panelNegGenesRand);
                Hashtable<String, String> predictedClassesRand = rcRand.predictClasses2(testMutationData, testSampleClasses);
                //System.out.println(predictedClasses);
                Hashtable<String, Double> accuracyResultsRand = rcRand.evaluateResults2(predictedClassesRand, testSampleClasses);
                //System.out.println(accuracyResults);

                //Update accuracy results for random variants
                double tpRand=accuracyResultsRand.get("TP");
                double tnRand=accuracyResultsRand.get("TN");
                double fpRand=accuracyResultsRand.get("FP");
                double fnRand=accuracyResultsRand.get("FN");
                avgTpRand+=tpRand;
                avgTnRand+=tnRand;
                avgFpRand+=fpRand;
                avgFnRand+=fnRand;
                avgUnclassifiedRand+=accuracyResultsRand.get("Unclassified");
                double precisionRand=(tpRand/(tpRand+fpRand))*100;
                avgPrecisionRand+=precisionRand;
                double recallRand=(tpRand/(tpRand+fnRand))*100;
                avgRecallRand+=recallRand;
                avgSpecificityRand+=(tnRand/(tnRand+fpRand))*100;
                avgFprRand+=(fpRand/(fpRand+tnRand))*100;
                avgFnrRand+=(fnRand/(fnRand+tpRand))*100;
                double accuracyRand=((tpRand+tnRand)/(tpRand+tnRand+fpRand+fnRand))*100;
                avgAccuracyRand+=accuracyRand;
                //System.out.println(precisionRand+"\t"+recallRand);
                double f1scoreRand;
                if(precisionRand==0.0 && recallRand==0.0)
                    f1scoreRand=0.0;
                else if(Double.isNaN(precisionRand) || Double.isNaN(recallRand))
                    f1scoreRand=0.0;
                else
                    f1scoreRand=2*(precisionRand*recallRand/(precisionRand+recallRand));
                avgF1scoreRand+=f1scoreRand;

                //Update pvalues
                if(accuracyRand>=accuracy)
                    pvalAccuracy++;
                if(f1scoreRand>=f1score)
                    pvalF1score++;

            }

            //Compute final accuracy results for random variants
            avgTpRand/=numVariants;
            avgTnRand/=numVariants;
            avgFpRand/=numVariants;
            avgFnRand/=numVariants;
            avgUnclassifiedRand/=numVariants;
            avgPrecisionRand/=numVariants;
            avgRecallRand/=numVariants;
            avgSpecificityRand/=numVariants;
            avgFprRand/=numVariants;
            avgFnrRand/=numVariants;
            avgAccuracyRand/=numVariants;
            avgF1scoreRand/=numVariants;

            //Compute final pvalues
            pvalAccuracy/=numVariants;
            pvalF1score/=numVariants;

            //Print results
            /*System.out.print(refGene+"\t"+tp+"\t"+tn+"\t"+fp+"\t"+fn+"\t"+unclassified+"\t"+precision+"\t"+recall+"\t"+specificity);
            System.out.print("\t"+fpr+"\t"+fnr+"\t"+accuracy+"\t"+f1score);
            System.out.print("\t"+avgTpRand+"\t"+avgTnRand+"\t"+avgFpRand+"\t"+avgFnRand+"\t"+avgUnclassifiedRand+"\t"+avgPrecisionRand);
            System.out.print("\t"+avgRecallRand+"\t"+avgSpecificityRand+"\t"+avgFprRand+"\t"+avgFnrRand+"\t"+avgAccuracyRand);
            System.out.println("\t"+avgF1scoreRand+"\t"+pvalAccuracy+"\t"+pvalF1score);*/

            bw=new BufferedWriter(new FileWriter(outputFolder+"/"+tumor+"_perm.txt",true));
            bw.write(refGene+"\t"+tp+"\t"+tn+"\t"+fp+"\t"+fn+"\t"+unclassified+"\t"+precision+"\t"+recall+"\t"+specificity);
            bw.write("\t"+fpr+"\t"+fnr+"\t"+accuracy+"\t"+f1score);
            bw.write("\t"+avgTpRand+"\t"+avgTnRand+"\t"+avgFpRand+"\t"+avgFnRand+"\t"+avgUnclassifiedRand+"\t"+avgPrecisionRand);
            bw.write("\t"+avgRecallRand+"\t"+avgSpecificityRand+"\t"+avgFprRand+"\t"+avgFnrRand+"\t"+avgAccuracyRand);
            bw.write("\t"+avgF1scoreRand+"\t"+pvalAccuracy+"\t"+pvalF1score+"\n");
            bw.close();

        }

    }
}

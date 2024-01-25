import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.*;

public class TestPanel
{
    public static void main(String[] args) throws Exception
    {
        //String dataset="DBGAP";
        String dataset="PCAWG-BRCA";
        String trainingFile="Data/BRCA_snp_gene_matrix.txt";
        String testFile="Data/"+dataset+"_gene_matrix.txt";
        String outputFile="Results/TestDataset/"+dataset+"_results.txt";
        //String weightsFile="Data/CDS_Lengths.txt";
        String weightsFile=null;
        //String listRefGenes=null;
        String listRefGenes="BRCA1,BRCA2";
        int[] posPanelSize={5,15,30,50};
        int[] negPanelSize={5,15,30,50};

        FileManager fm=new FileManager();

        //System.out.println("Reading training data...");
        MutationMatrix trainingMatrix;
        if(listRefGenes==null)
            trainingMatrix=fm.readMutationDataWithClasses(trainingFile);
        else
        {
            String[] split=listRefGenes.split(",");
            HashSet<String> setRefGenes = new HashSet<>(Arrays.asList(split));
            trainingMatrix=fm.readMutationDataWithoutClasses(trainingFile,setRefGenes);
        }
        HashMap<String,HashSet<String>> trainingData=trainingMatrix.getMutationData();
        Hashtable<String,String> mapTrainingSampleClasses=trainingMatrix.getMapSampleClasses();
        HashSet<String> setTrainingPositives=trainingMatrix.getPositiveSet();
        int numTrainingPositives=setTrainingPositives.size();
        System.out.println(numTrainingPositives);
        HashSet<String> setTrainingNegatives=trainingMatrix.getNegativeSet();
        int numTrainingNegatives=setTrainingNegatives.size();
        System.out.println(numTrainingNegatives);

        //System.out.println("Reading test data...");
        MutationMatrix testMatrix;
        if(listRefGenes==null)
            testMatrix=fm.readMutationDataWithClasses(testFile);
        else
        {
            String[] split=listRefGenes.split(",");
            HashSet<String> setRefGenes = new HashSet<>(Arrays.asList(split));
            testMatrix=fm.readMutationDataWithoutClasses(testFile,setRefGenes);
        }
        Hashtable<String,String> mapTestSampleClasses=testMatrix.getMapSampleClasses();
        int numTestSamples=mapTestSampleClasses.size();
        HashMap<String,HashSet<String>> testData=testMatrix.getMutationData();
        HashSet<String> setTestPositives=testMatrix.getPositiveSet();
        int numTestPositives=setTestPositives.size();
        System.out.println(numTestPositives);
        HashSet<String> setTestNegatives=testMatrix.getNegativeSet();
        int numTestNegatives=setTestNegatives.size();
        System.out.println(numTestNegatives);


        HashMap<String,Double> mapWeights;
        if(weightsFile==null)
            mapWeights = null;
        else
            mapWeights=fm.readWeightsFile(weightsFile);

        BufferedWriter bw=new BufferedWriter(new FileWriter(outputFile));
        bw.write("(PosRuleSize,NegRuleSize)\tTP\tTN\tFP\tFN\tUnclassified\tPrecision\tRecall\tSpecificity\tFPR\tFNR\tAccuracy\tF1Score\n");

        for(int i=0;i<posPanelSize.length;i++)
        {
            for(int j=0;j<negPanelSize.length;j++)
            {
                //System.out.println("Extract mutation panel for positives...");
                DMGSFinderApp mc=new DMGSFinderApp(trainingData,mapTrainingSampleClasses,mapWeights,true,posPanelSize[i],1.0);
                Vector<String> panelPosGenes=mc.runAlgorithm();

                //System.out.println("Extract mutation panel for negatives...");
                mc=new DMGSFinderApp(trainingData,mapTrainingSampleClasses,mapWeights,false,negPanelSize[j],1.0);
                Vector<String> panelNegGenes=mc.runAlgorithm();

                //System.out.println("Classifying test samples...");
                RuleClassifier rc=new RuleClassifier(panelPosGenes,panelNegGenes);
                Hashtable<String,String> predictedClasses=rc.predictClasses2(testData,mapTestSampleClasses);
                Hashtable<String,Double> accuracyResults=rc.evaluateResults2(predictedClasses,mapTestSampleClasses);

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

                //Write accuracy results to output file
                bw.write("("+posPanelSize[i]+","+negPanelSize[j]+")\t"+tp+"\t"+tn+"\t"+fp+"\t"+fn+"\t"+(unclassified/numTestSamples*100)+"\t");
                bw.write(precision+"\t"+recall+"\t"+specificity+"\t"+fpr+"\t"+fnr+"\t"+accuracy+"\t"+f1score+"\n");

                /*System.out.println("------------------------");
                System.out.println("POSITIVE RULE SIZE: "+posPanelSize[i]);
                System.out.println("NEGATIVE RULE SIZE: "+negPanelSize[j]);
                System.out.println("PRECISION: "+accuracyResults.get("Precision"));
                System.out.println("SENSITIVITY/RECALL/TPR: "+accuracyResults.get("Sensitivity"));
                System.out.println("FPR: "+accuracyResults.get("FPR"));
                System.out.println("FNR: "+accuracyResults.get("FNR"));
                System.out.println("SPECIFICITY: "+accuracyResults.get("Specificity"));
                System.out.println("ACCURACY: "+accuracyResults.get("Accuracy"));
                System.out.println("BALANCED ACCURACY: "+accuracyResults.get("Balanced accuracy"));
                System.out.println("PERCENTAGE OF UNCLASSIFIED: "+(accuracyResults.get("Unclassified")/numTestSamples*100));
                System.out.println("NUMBER OF TRUE POSITIVES: "+accuracyResults.get("TP"));
                System.out.println("NUMBER OF FALSE POSITIVES: "+accuracyResults.get("FP"));
                System.out.println("NUMBER OF TRUE NEGATIVES: "+accuracyResults.get("TN"));
                System.out.println("NUMBER OF FALSE NEGATIVES: "+accuracyResults.get("FN"));*/
            }
        }
        bw.close();

    }
}

import java.util.*;

public class TestMUTClassCV
{
    public static void main(String[] args) throws Exception
    {
        String tumor=args[0];
        String trainingFile="C:/Ricerca/Temp/"+tumor+"_trainingSet.txt";
        String testFile="C:/Ricerca/Temp/"+tumor+"_testSet.txt";
        int posPanelSize=5;
        int negPanelSize=50;
        String predictionsFile="C:/Ricerca/Temp/"+tumor+"_predictions.txt";

        //Read training and test set
        FileManager fm=new FileManager();
        MutationMatrix trainingSet=fm.readMutationDataCV(trainingFile);
        MutationMatrix testSet=fm.readMutationDataCV(testFile);

        //Run MUTClass on training set
        long inizio=System.currentTimeMillis();
        DMGSFinderApp mc = new DMGSFinderApp(trainingSet.getMutationData(), trainingSet.getMapSampleClasses(), null, true, posPanelSize, 1.0);
        Vector<String> panelPosGenes = mc.runAlgorithm();
        mc = new DMGSFinderApp(trainingSet.getMutationData(), trainingSet.getMapSampleClasses(), null, false, negPanelSize, 1.0);
        Vector<String> panelNegGenes = mc.runAlgorithm();
        double fine=System.currentTimeMillis();
        double runningTime=(fine-inizio)/1000;

        //Run MUTClass on test set
        RuleClassifier rc = new RuleClassifier(panelPosGenes, panelNegGenes);
        Hashtable<String, String> predictedClasses = rc.predictClasses2(testSet.getMutationData(), testSet.getMapSampleClasses());

        //Write predictions on output file
        fm.writePredictions(testSet.getMapSampleClasses(), predictedClasses, runningTime, predictionsFile);

    }

}

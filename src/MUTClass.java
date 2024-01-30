import java.util.*;

public class MUTClass
{
    public static void main(String[] args)
    {
        //Algorithm parameters
        String trainingFile=null;
        String testFile=null;
        String listDriverGenes=null;
        String weightsFile=null;
        int posPanelSize=5;
        int negPanelSize=50;
        String resultsFile="results.txt";

        int i;
        for (i=0;i<args.length;i++)
        {
            switch (args[i])
            {
                case "-train" -> trainingFile = args[++i];
                case "-test" -> testFile = args[++i];
                case "-d" -> listDriverGenes = args[++i];
                case "-w" -> weightsFile = args[++i];
                case "-kmax" -> posPanelSize = Integer.parseInt(args[++i]);
                case "-kmin" -> negPanelSize = Integer.parseInt(args[++i]);
                case "-o" -> resultsFile = args[++i];
                default -> {
                    System.out.println("Error! Unrecognizable command '" + args[i] + "'");
                    printHelp();
                    System.exit(1);
                }
            }
        }

        //Error in case training or test files are missing or wrong
        if(trainingFile==null || testFile==null)
        {
            System.out.println("Error! No file for mutation matrix has been specified!");
            printHelp();
            System.exit(1);
        }

        if(listDriverGenes!=null)
            System.out.println("DRIVER GENES: "+listDriverGenes);
        else
            System.out.println("DRIVER GENES: none");
        System.out.println("SIZE OF POSITIVE PANEL: "+posPanelSize);
        System.out.println("SIZE OF NEGATIVE PANEL: "+negPanelSize);

        FileManager fm=new FileManager();

        System.out.println("Read training data...");
        MutationMatrix trainingMatrix;
        if(listDriverGenes==null)
            trainingMatrix=fm.readMutationDataWithClasses(trainingFile);
        else
        {
            String[] split=listDriverGenes.split(",");
            HashSet<String> setRefGenes = new HashSet<>(Arrays.asList(split));
            trainingMatrix=fm.readMutationDataWithoutClasses(trainingFile,setRefGenes);
        }
        HashMap<String,HashSet<String>> trainingData=trainingMatrix.getMutationData();
        Hashtable<String,String> mapTrainingSampleClasses=trainingMatrix.getMapSampleClasses();
        HashMap<String,Double> mapWeights;
        if(weightsFile==null)
            mapWeights = null;
        else
            mapWeights=fm.readWeightsFile(weightsFile);

        System.out.println("Read test data...");
        MutationMatrix testMatrix=fm.readTestData(testFile);
        HashMap<String,HashSet<String>> testData=testMatrix.getMutationData();
        Hashtable<String,String> mapTestSampleClasses=testMatrix.getMapSampleClasses();

        System.out.println("Extract mutation panel for positives...");
        DMGSFinderApp mc=new DMGSFinderApp(trainingData,mapTrainingSampleClasses,mapWeights,true,posPanelSize);
        Vector<String> panelPosGenes=mc.runAlgorithm();
        System.out.println("Extract mutation panel for negatives...");
        mc=new DMGSFinderApp(trainingData,mapTrainingSampleClasses,mapWeights,false,negPanelSize);
        Vector<String> panelNegGenes=mc.runAlgorithm();
        System.out.println("Classify test samples...");
        RuleClassifier rc=new RuleClassifier(panelPosGenes,panelNegGenes);
        Hashtable<String,String> predictedClasses=rc.predictClasses(testData,mapTestSampleClasses);

        System.out.println("Write results to output file...");
        fm.writeMUTClassResults(resultsFile,panelPosGenes,panelNegGenes,predictedClasses);

    }

    private static void printHelp()
    {
        String help = "Usage: java -cp ./out MUTClass -train <trainingSetFile> -test <testSetFile> [-d <listDriverGenes> " +
                "-kmax <positivePanelSize> -kmin <negativePanelSize> -o <resultsFile>\n\n";
        help+="REQUIRED PARAMETERS:\n";
        help+="-train\tTraining mutation matrix file\n\n";
        help+="-test\tTest mutation matrix file\n\n";
        help+="OPTIONAL PARAMETERS:\n";
        help+="-d\tList of driver genes (default=training matrix file with classes provided as input)\n";
        help+="-kmax\tSize of positive gene panel (default=5)\n";
        help+="-kmin\tSize of negative gene panel (default=50)\n";
        help+="-o\tOutput result file (default='results.txt')\n\n";
        System.out.println(help);
    }

}

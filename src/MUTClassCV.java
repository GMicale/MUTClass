import java.util.*;

public class MUTClassCV
{
    public static void main(String[] args)
    {
        //Algorithm parameters
        String mutationsFile=null;
        String listDriverGenes=null;
        String weightsFile=null;
        int maxNumSolGenesPos=10;
        int maxNumSolGenesNeg=10;
        int numIterCV=10;
        int numFolds=5;

        int i;
        for (i=0;i<args.length;i++)
        {
            switch (args[i])
            {
                case "-m" -> mutationsFile = args[++i];
                case "-d" -> listDriverGenes = args[++i];
                case "-w" -> weightsFile = args[++i];
                case "-kmax" -> maxNumSolGenesPos = Integer.parseInt(args[++i]);
                case "-kmin" -> maxNumSolGenesNeg = Integer.parseInt(args[++i]);
                case "-cv" -> numIterCV = Integer.parseInt(args[++i]);
                case "-f" -> numFolds = Integer.parseInt(args[++i]);
                default -> {
                    System.out.println("Error! Unrecognizable command '" + args[i] + "'");
                    printHelp();
                    System.exit(1);
                }
            }
        }

        //Error in case mutation matrix file and/or info about positive samples or gene types are missing or wrong
        if(mutationsFile==null)
        {
            System.out.println("Error! No file for mutation matrix has been specified!");
            printHelp();
            System.exit(1);
        }

        if(listDriverGenes!=null)
            System.out.println("DRIVER GENES: "+listDriverGenes);
        System.out.println("CV ITERATIONS: "+numIterCV);
        System.out.println("CV FOLDS: "+numFolds);
        System.out.println("GENES IN POSITIVE RULE: "+maxNumSolGenesPos);
        System.out.println("GENES IN NEGATIVE RULE: "+maxNumSolGenesNeg);

        //Read mutation data
        //System.out.println("Reading mutation matrix...");
        FileManager fm=new FileManager();
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
        HashSet<String> setPositives=mutMatrix.getPositiveSet();
        HashSet<String> setNegatives=mutMatrix.getNegativeSet();
        Vector<String> listPositives=new Vector<>(setPositives);
        Vector<String> listNegatives=new Vector<>(setNegatives);
        //System.out.println(numPositivesTraining);
        //System.out.println(numNegativesTraining);

        //Read gene weights (if specified)
        HashMap<String,Double> mapWeights;
        if(weightsFile==null)
            mapWeights = null;
        else
            mapWeights=fm.readWeightsFile(weightsFile);

        int j, k;
        double avgAccuracy=0.0;
        double avgPrecision=0.0;
        double avgRecall=0.0;
        double avgFpr=0.0;
        double avgFnr=0.0;
        double avgSpecificity=0.0;
        double avgUnclassified=0.0;
        double avgTruePositives=0.0;
        double avgTrueNegatives=0.0;
        double avgFalsePositives=0.0;
        double avgFalseNegatives=0.0;
        double avgF1=0.0;

        for(i=0;i<numIterCV;i++)
        {
            //System.out.println("CV iteration: "+(i+1)+"/"+numIterCV);

            //Shuffle the sets of positive and negative samples
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
                        trainingSetSamples.addAll(posFolds.get(j));
                        trainingSetSamples.addAll(negFolds.get(j));
                    }
                }
                HashMap<String,HashSet<String>> subMutationData=Utility.extractSubMutData(mutationData,trainingSetSamples);
                Hashtable<String,String> subSampleClasses=Utility.extractSubSampleClasses(mapSampleClasses,trainingSetSamples);

                //Run MUTClass to find panels of mutated genes
                DMGSFinderApp mc=new DMGSFinderApp(subMutationData,subSampleClasses,mapWeights,true,maxNumSolGenesPos);
                Vector<String> panelPosGenes=mc.runAlgorithm();
                mc=new DMGSFinderApp(subMutationData,subSampleClasses,mapWeights,false,maxNumSolGenesNeg);
                Vector<String> panelNegGenes=mc.runAlgorithm();
                //System.out.println("Panel sizes: "+panelPosGenes.size()+"\t"+panelNegGenes.size());

                //Make classification
                RuleClassifier rc=new RuleClassifier(panelPosGenes,panelNegGenes);
                Hashtable<String,String> subPredictedClasses=rc.predictClasses(subMutationData,subSampleClasses);
                Hashtable<String,Double> accuracyResults=rc.evaluateResults(subPredictedClasses,subSampleClasses);

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
            System.out.println("FP: "+falsePositives);*/
            double precision=truePositives/(truePositives+falsePositives);
            double recall=truePositives/(truePositives+falseNegatives);
            double specificity=trueNegatives/(trueNegatives+falsePositives);
            double fpr=falsePositives/(falsePositives+trueNegatives);
            double fnr=falseNegatives/(falseNegatives+truePositives);
            double accuracy=(truePositives+trueNegatives)/(truePositives+trueNegatives+falsePositives+falseNegatives);
            double f1=2*(precision*recall/(precision+recall));

            /*System.out.println("Accuracy: "+accuracy);
            System.out.println("Precision: "+precision);
            System.out.println("Sensitivity/Recall/TPR: "+recall);
            System.out.println("Specificity: "+specificity);
            System.out.println("FPR: "+fpr);*/

            //Update average statistics
            avgAccuracy+=accuracy;
            avgPrecision+=precision;
            avgRecall+=recall;
            avgFpr+=fpr;
            avgFnr+=fnr;
            avgSpecificity+=specificity;
            avgUnclassified+=unclassified;
            avgTruePositives+=truePositives;
            avgTrueNegatives+=trueNegatives;
            avgFalsePositives+=falsePositives;
            avgFalseNegatives+=falseNegatives;
            avgF1+=f1;

        }

        avgAccuracy/=numIterCV;
        avgPrecision/=numIterCV;
        avgRecall/=numIterCV;
        avgFpr/=numIterCV;
        avgFnr/=numIterCV;
        avgSpecificity/=numIterCV;
        avgUnclassified/=numIterCV;
        avgTruePositives/=numIterCV;
        avgTrueNegatives/=numIterCV;
        avgFalsePositives/=numIterCV;
        avgFalseNegatives/=numIterCV;
        avgF1/=numIterCV;

        System.out.println("AVERAGE TRUE POSITIVES: "+((double)Math.round(avgTruePositives*100))/100);
        System.out.println("AVERAGE TRUE NEGATIVES: "+((double)Math.round(avgTrueNegatives*100))/100);
        System.out.println("AVERAGE FALSE POSITIVES: "+((double)Math.round(avgFalsePositives*100))/100);
        System.out.println("AVERAGE FALSE NEGATIVES: "+((double)Math.round(avgFalseNegatives*100))/100);
        System.out.println("AVERAGE UNCLASSIFIED: "+((double)Math.round(avgUnclassified*100))/100);
        System.out.println("AVERAGE PRECISION: "+((double)Math.round(avgPrecision*10000))/100+"%");
        System.out.println("AVERAGE SENSITIVITY/RECALL/TPR: "+((double)Math.round(avgRecall*10000))/100+"%");
        System.out.println("AVERAGE FPR: "+((double)Math.round(avgFpr*10000))/100+"%");
        System.out.println("AVERAGE FNR: "+((double)Math.round(avgFnr*10000))/100+"%");
        System.out.println("AVERAGE SPECIFICITY: "+((double)Math.round(avgSpecificity*10000))/100+"%");
        System.out.println("AVERAGE ACCURACY: "+((double)Math.round(avgAccuracy*10000))/100+"%");
        System.out.println("AVERAGE F1 SCORE: "+((double)Math.round(avgF1*10000))/100+"%");
        System.out.println("------------------------");
    }

    private static void printHelp()
    {
        String help = "Usage: java -cp ./out MUTClassCV -m <mutationsFile> [-d <listDriverGenes> " +
                "-kmax <positivePanelSize> -kmin <negativePanelSize> " +
                "-cv <crossValidationIterations> -f <crossValidationFolds>]\n\n";
        help+="REQUIRED PARAMETERS:\n";
        help+="-m\tMutation matrix file\n\n";
        help+="OPTIONAL PARAMETERS:\n";
        help+="-d\tList of driver genes (default=matrix file with classes provided as input)\n";
        help+="-kmax\tSize of positive gene panel (default=10)\n";
        help+="-kmin\tSize of negative gene panel (default=10)\n";
        help+="-cv\tNumber of iterations of cross validation (default=10)\n";
        help+="-f\tNumber of folds for cross validation (default=5)\n\n";
        System.out.println(help);
    }

}

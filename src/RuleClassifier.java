import java.util.HashMap;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.List;

public class RuleClassifier
{
    private final List<String> posPanel;
    private final List<String> negPanel;

    public RuleClassifier(List<String> posPanel, List<String> negPanel)
    {
        this.posPanel=posPanel;
        this.negPanel=negPanel;
    }

    public List<String> getPosPanel()
    {
        return posPanel;
    }

    public Hashtable<String,String> predictClasses(HashMap<String,HashSet<String>> mutData, Hashtable<String,String> mapRealClasses)
    {
        Hashtable<String,String> mapPredictedClasses=new Hashtable<>();
        for(String sample : mapRealClasses.keySet())
        {
            double numCovGenesPosRule = 0.0;
            double numCovGenesNegRule = 0.0;
            String mutClass = "P";
            for (String gene : posPanel)
            {
                if (mutData.containsKey(gene) && mutData.get(gene).contains(sample))
                    numCovGenesPosRule++;
            }
            //double posScore=(numCovGenesPosRule/posPanel.size())*Math.log(posPanel.size());
            double posScore = numCovGenesPosRule/posPanel.size();
            for (String gene : negPanel)
            {
                if (mutData.containsKey(gene) && mutData.get(gene).contains(sample))
                    numCovGenesNegRule++;
            }
            //double negScore=(numCovGenesNegRule/negPanel.size())*Math.log(negPanel.size());
            double negScore = numCovGenesNegRule/negPanel.size();
            //percCovGenesNegRule /= negPanel.size();
            //System.out.println(posScore+" --- "+negScore);
            if (negScore > posScore)
                mutClass = "N";
            else if(negScore==posScore)
                mutClass = "U";
            //System.out.println(posScore+"/"+negScore+" --> ");
            mapPredictedClasses.put(sample,mutClass);
        }
        return mapPredictedClasses;
    }

    public Hashtable<String,Double> evaluateResults(Hashtable<String,String> predictedClasses, Hashtable<String,String> trueClasses)
    {
        Hashtable<String,Double> listResults=new Hashtable<>();
        double truePositives=0.0;
        double trueNegatives=0.0;
        double falsePositives=0.0;
        double falseNegatives=0.0;
        double numUnclassified=0.0;

        for(String sample : predictedClasses.keySet())
        {
            String predClass=predictedClasses.get(sample);
            String trueClass=trueClasses.get(sample);
            if(predClass.equals("U"))
                numUnclassified++;
            else
            {
                if(predClass.equals(trueClass))
                {
                    if(predClass.equals("P"))
                        truePositives++;
                    else
                        trueNegatives++;
                }
                else
                {
                    if(predClass.equals("N"))
                        falseNegatives++;
                    else
                        falsePositives++;
                }
            }
        }
        /*System.out.println("TP: "+truePositives);
        System.out.println("FN: "+falseNegatives);
        System.out.println("TN: "+trueNegatives);
        System.out.println("FP: "+falsePositives);*/
        /*accuracy=(accuracy/numClassified)*100;
        double precision=(truePositives/(truePositives+falsePositives))*100;
        double recall=(truePositives/(truePositives+falseNegatives))*100;
        double specificity=(trueNegatives/(trueNegatives+falsePositives))*100;
        double fpr=(falsePositives/(falsePositives+trueNegatives))*100;
        double fnr=(falseNegatives/(falseNegatives+truePositives))*100;
        double balancedAccuracy=(recall+specificity)/2;*/
        listResults.put("TP",truePositives);
        listResults.put("TN",trueNegatives);
        listResults.put("FP",falsePositives);
        listResults.put("FN",falseNegatives);
        listResults.put("Unclassified",numUnclassified);
        /*listResults.put("Accuracy",accuracy);
        listResults.put("Precision",precision);
        listResults.put("Sensitivity",recall);
        listResults.put("Recall",recall);
        listResults.put("TPR",recall);
        listResults.put("Specificity",specificity);
        listResults.put("FPR",fpr);
        listResults.put("FNR",fnr);
        listResults.put("Balanced accuracy",balancedAccuracy);*/
        return listResults;
    }

}

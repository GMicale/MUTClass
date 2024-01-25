import java.util.*;

public class Utility
{

    public static HashSet<String> getCandidateGenes(HashMap<String,HashSet<String>> mutData, HashSet<String> setPositives,
                                                    HashSet<String> setNegatives, boolean maximization)
    {
        HashSet<String> setCandidateGenes=new HashSet<>();
        for(String gene : mutData.keySet())
        {
            Set<String> commonSamples = new HashSet<String>(mutData.get(gene));
            if(maximization)
                commonSamples.retainAll(setPositives);
            else
                commonSamples.retainAll(setNegatives);
            if(commonSamples.size()>0)
                setCandidateGenes.add(gene);
        }
        return setCandidateGenes;
    }

    public static HashSet<String> getCandidateGenes(HashMap<String,HashSet<String>> mutData, Hashtable<String,String> mapSampleClasses, boolean maximization)
    {
        HashSet<String> setCandidateGenes=new HashSet<>();
        for(String gene : mutData.keySet())
        {
            HashSet<String> setMutatedSamples=mutData.get(gene);
            String refClass="P";
            if(!maximization)
                refClass="N";
            int numRefSamples=0;
            for(String sample: setMutatedSamples)
            {
                if(mapSampleClasses.get(sample).equals(refClass))
                    numRefSamples++;
            }
            if(numRefSamples>0)
                setCandidateGenes.add(gene);
        }
        return setCandidateGenes;
    }

    public static Vector<String> findBestSol(HashSet<String> setUnexploredGenes, HashMap<String,HashSet<String>> mapMutations,
                                             Hashtable<String,String> mapSampleClasses, int currTotalPosCov,
                                             int currTotalNegCov, int cardinalitySol, boolean maximization, int numPositives,
                                              int numNegatives, double alpha, HashMap<String,Double> mapWeights, double totalCandWeight)
    {
        double bestWeight=0.0;
        double bestScore=Double.MAX_VALUE;
        if(maximization)
            bestScore=-Double.MAX_VALUE;
        double bestAvgPosCov=0.0;
        double bestAvgNegCov=0.0;
        String bestGene="";
        for(String gene : setUnexploredGenes)
        {
            //Compute the average coverage of the gene in the positive and negative sets
            HashSet<String> setMutatedSamples=mapMutations.get(gene);
            double newTotalPosCov=currTotalPosCov;
            double newTotalNegCov=currTotalNegCov;
            for(String sample : setMutatedSamples)
            {
                String sampleClass=mapSampleClasses.get(sample);
                if(sampleClass.equals("P"))
                    newTotalPosCov++;
                else
                    newTotalNegCov++;
            }
            double avgPosCov=newTotalPosCov/((cardinalitySol+1)*numPositives);
            double avgNegCov=newTotalNegCov/((cardinalitySol+1)*numNegatives);

            //Compute differential coverage of the gene
            double covDiff=avgPosCov-avgNegCov;

            //Get gene's weight
            double geneWeight=0.0;
            if(mapWeights!=null)
                //geneWeight=Math.log(mapWeights.get(gene))/Math.log(totalCandWeight);
                //geneWeight=mapWeights.get(gene)/totalCandWeight;
                geneWeight=1/mapWeights.get(gene);
                //geneWeight=1/Math.log(mapWeights.get(gene));


            //Compute gene score
            //double score=alpha*covDiff - (1.0-alpha)*geneWeight;
            double score=covDiff;
            if(mapWeights!=null)
                score=score*geneWeight;
            //System.out.println(gene);

            //Update best solution, if necessary
            if(maximization)
            {
                if(score>bestScore)
                {
                    bestWeight=geneWeight;
                    bestScore=score;
                    bestAvgPosCov=avgPosCov;
                    bestAvgNegCov=avgNegCov;
                    bestGene=gene;
                }
            }
            else
            {
                if(score<bestScore)
                {
                    bestWeight=geneWeight;
                    bestScore=score;
                    bestAvgPosCov=avgPosCov;
                    bestAvgNegCov=avgNegCov;
                    bestGene=gene;
                }
            }
        }
        Vector<String> bestRes=new Vector<>();
        bestRes.add(bestGene);
        bestRes.add(String.valueOf(bestAvgPosCov));
        bestRes.add(String.valueOf(bestAvgNegCov));
        bestRes.add(String.valueOf(bestWeight));
        bestRes.add(String.valueOf(bestScore));
        return bestRes;
    }

    public static HashMap<String,HashSet<String>> extractSubMutData(HashMap<String,HashSet<String>> mutData,
                                                                    HashSet<String> subsetSamples)
    {
        HashMap<String,HashSet<String>> subMutData=new HashMap<>();
        for(String gene : mutData.keySet())
        {
            HashSet<String> setSamples=mutData.get(gene);
            for(String sample : setSamples)
            {
                if(subsetSamples.contains(sample))
                {
                    if(subMutData.containsKey(gene))
                        subMutData.get(gene).add(sample);
                    else
                    {
                        HashSet<String> subSetSamples=new HashSet<>();
                        subSetSamples.add(sample);
                        subMutData.put(gene,subSetSamples);
                    }
                }
            }
        }
        return subMutData;
    }

    public static Hashtable<String,String> extractSubSampleClasses(Hashtable<String,String> mapSampleClasses,
                                                                    HashSet<String> subsetSamples)
    {
        Hashtable<String,String> subSampleClasses=new Hashtable<>();
        for(String sample : subsetSamples)
            subSampleClasses.put(sample,mapSampleClasses.get(sample));
        return subSampleClasses;
    }

    public static void normalizeWeights(HashMap<String,Double> mapWeights)
    {
        double total=0.0;
        for(double weight : mapWeights.values())
            total+=weight;
        for(String gene : mapWeights.keySet())
            mapWeights.put(gene,mapWeights.get(gene)/total);
    }

}

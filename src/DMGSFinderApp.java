import java.util.*;

public class DMGSFinderApp
{
    private final HashMap<String,HashSet<String>> mutationData;
    private final Hashtable<String,String> mapSampleClasses;
    private final HashMap<String,Double> mapWeights;
    private final boolean maximization;
    private final int maxNumSolGenes;
    private final double alpha;

    public DMGSFinderApp(HashMap<String,HashSet<String>> mutationData, Hashtable<String,String> mapSampleClasses,
                       HashMap<String,Double> mapWeights, boolean maximization, int maxNumSolGenes, double alpha)
    {
        this.mutationData=mutationData;
        this.mapSampleClasses=mapSampleClasses;
        this.mapWeights=mapWeights;
        this.maximization=maximization;
        this.maxNumSolGenes=maxNumSolGenes;
        this.alpha=alpha;
    }

    public Vector<String> runAlgorithm()
    {
        //Build the set of candidate genes
        HashSet<String> setGenes;
        double totalCandWeight=0.0;
        if(mapWeights==null)
            //setGenes=new HashSet<>(mutationData.keySet());
            setGenes=Utility.getCandidateGenes(mutationData,mapSampleClasses,maximization);
        else
        {
            setGenes=new HashSet<>();
            Set<String> mutGenes=mutationData.keySet();
            for(String gene : mutGenes)
            {
                if(mapWeights.containsKey(gene))
                {
                    setGenes.add(gene);
                    totalCandWeight+=mapWeights.get(gene);
                }

            }
        }

        //Initialize data structures
        Vector<String> currSol=new Vector<>();
        int currTotalPosCov=0;
        int currTotalNegCov=0;

        //Count the number of positives and negatives
        int numPositives=0;
        int numNegatives=0;
        for(String sample : mapSampleClasses.keySet())
        {
            if(mapSampleClasses.get(sample).equals("P"))
                numPositives++;
            else
                numNegatives++;
        }
        //System.out.println(numPositives);
        //System.out.println(numNegatives);
        if(numPositives==0 || numNegatives==0)
            return null;

        //Run greedy algorithm
        while(currSol.size()<maxNumSolGenes)
        {
            //Find the gene that, if added to current solution, maximizes (minimizes) the differential coverage
            Vector<String> bestRes=Utility.findBestSol(setGenes,mutationData,mapSampleClasses,currTotalPosCov,currTotalNegCov,
                    currSol.size(),maximization,numPositives,numNegatives,alpha,mapWeights,totalCandWeight);
            String bestGene=bestRes.get(0);

            //Update current solution
            setGenes.remove(bestGene);
            currSol.add(bestGene);
            //System.out.println("Best gene: "+bestGene);
            HashSet<String> setMutatedSamples=mutationData.get(bestGene);
            for(String sample : setMutatedSamples)
            {
                String sampleClass=mapSampleClasses.get(sample);
                if(sampleClass.equals("P"))
                    currTotalPosCov++;
                else
                    currTotalNegCov++;
            }

        }

        //System.out.println(currSol+"\t"+lastBestCov);
        return currSol;

    }
}

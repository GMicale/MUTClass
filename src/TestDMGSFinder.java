import java.util.*;

public class TestDMGSFinder
{
   public static void main(String[] args)
   {
       //String inputFile="C:/Ricerca/MUTClass/Dataset/BRCA/BRCA_snp_gene_matrix_Main_Classification_LumA.txt";
       String inputFile="Data/BRCA_snp_gene_matrix.txt";
       String listRefGenes="BRCA1,BRCA2";
       //String listRefGenes=null;
       //String weightsFile="Data/CDS_Lengths.txt";
       String weightsFile=null;
       boolean maximization=true;
       int maxNumSolGenes=15;
       double alpha=0.8;
       //String outputFile=null;
       //String outputFile="Results/BRCA_maximal_itemsets_SNP.txt";

       //Read mutation data
       System.out.println("Reading mutation matrix...");
       FileManager fm=new FileManager();
       //System.out.println(setPositives.size());
       MutationMatrix mutMatrix;
       if(listRefGenes==null)
           mutMatrix=fm.readMutationDataWithClasses(inputFile);
       else
       {
           String[] split=listRefGenes.split(",");
           HashSet<String> setRefGenes = new HashSet<>(Arrays.asList(split));
           mutMatrix=fm.readMutationDataWithoutClasses(inputFile,setRefGenes);
       }
       HashMap<String,HashSet<String>> mutationData=mutMatrix.getMutationData();
       Hashtable<String,String> mapSampleClasses=mutMatrix.getMapSampleClasses();
       HashSet<String> setPositives=mutMatrix.getPositiveSet();
       int numPositives=setPositives.size();
       System.out.println(numPositives);
       HashSet<String> setNegatives=mutMatrix.getNegativeSet();
       int numNegatives=setNegatives.size();
       System.out.println(numNegatives);

       //Read gene weights (if specified)
       HashMap<String,Double> mapWeights;
       if(weightsFile==null)
       {
           mapWeights = null;
           alpha=1.0;
       }
       else
           mapWeights=fm.readWeightsFile(weightsFile);

       //Build the set of candidate genes
       HashSet<String> setGenes;
       double totalCandWeight=0.0;
       if(mapWeights==null)
            //setGenes=mutMatrix.getAllGenes();
            setGenes=Utility.getCandidateGenes(mutationData,setPositives,setNegatives,maximization);
       else
       {
           setGenes=new HashSet<>();
           HashSet<String> mutGenes=mutMatrix.getAllGenes();
           for(String gene : mutGenes)
           {
               if(mapWeights.containsKey(gene))
               {
                   setGenes.add(gene);
                   totalCandWeight+=mapWeights.get(gene);
               }
           }
       }
       System.out.println("Number of candidate genes: "+setGenes.size());

       //Initialize data structures
       Vector<String> currSol=new Vector<>();
       int currTotalPosCov=0;
       int currTotalNegCov=0;

       //Run greedy algorithm
       double lastBestCov=1.0;
       double lastBestAvgPosCov=1.0;
       double lastBestAvgNegCov=1.0;
       double totalWeight=0.0;
       double inizio=System.currentTimeMillis();
       while(currSol.size()<maxNumSolGenes)
       {
           //Find the gene that, if added to current solution, maximizes (minimizes) the differential coverage
           Vector<String> bestRes=Utility.findBestSol(setGenes,mutationData,mapSampleClasses,currTotalPosCov,currTotalNegCov,
                   currSol.size(),maximization,numPositives,numNegatives,alpha,mapWeights,totalCandWeight);
           String bestGene=bestRes.get(0);
           double bestAvgPosCov=Double.parseDouble(bestRes.get(1));
           double bestAvgNegCov=Double.parseDouble(bestRes.get(2));
           double bestWeight=0.0;
           System.out.println("Best gene: "+bestGene);
           if(mapWeights!=null)
                bestWeight=mapWeights.get(bestGene);
           double bestDiffCov=bestAvgPosCov-bestAvgNegCov;
           System.out.println(bestGene);
           System.out.println(bestDiffCov+"\t"+bestWeight);

           //Update current solution
           setGenes.remove(bestGene);
           currSol.add(bestGene);
           HashSet<String> setMutatedSamples=mutationData.get(bestGene);
           for(String sample : setMutatedSamples)
           {
               String sampleClass=mapSampleClasses.get(sample);
               if(sampleClass.equals("P"))
                   currTotalPosCov++;
               else
                   currTotalNegCov++;
           }
           lastBestCov=bestDiffCov;
           lastBestAvgPosCov=bestAvgPosCov;
           lastBestAvgNegCov=bestAvgNegCov;
           totalWeight+=bestWeight;
       }
       System.out.println(currSol+"\t"+lastBestCov+"\t"+lastBestAvgPosCov+"\t"+lastBestAvgNegCov+"\t"+totalWeight);
       double fine=System.currentTimeMillis();
       double totalTime=(fine-inizio)/1000;
       System.out.println("Time elapsed: "+totalTime+" secs");

   }
}

import java.util.*;

public class DMGSFinder
{
   public static void main(String[] args)
   {
       String mutationsFile=null;
       String listDriverGenes=null;
       String weightsFile=null;
       int maxNumSolGenes=10;
       double alpha=0.5;
       String outputFile=null;
       boolean maximization=true;

       int i;
       for (i=0;i<args.length;i++)
       {
           switch (args[i])
           {
               case "-m" -> mutationsFile = args[++i];
               case "-d" -> listDriverGenes = args[++i];
               case "-w" -> weightsFile = args[++i];
               case "-ms" -> maxNumSolGenes = Integer.parseInt(args[++i]);
               case "-a" -> alpha = Double.parseDouble(args[++i]);
               case "-o" -> outputFile = args[++i];
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

       //Read mutation data
       System.out.println("\nReading mutation matrix...");
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
       int numPositives=setPositives.size();
       int numNegatives=setNegatives.size();

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

       //Initialize data structures
       Vector<String> currSol=new Vector<>();
       int currTotalPosCov=0;
       int currTotalNegCov=0;

       //Run greedy algorithm
       System.out.println("Start algorithm...");
       double lastPosCov=1.0;
       double lastNegCov=1.0;
       double lastDiffCov=1.0;
       double totalWeight=0.0;
       int iter=1;
       while(currSol.size()<maxNumSolGenes)
       {
           System.out.print("Step "+iter+"/"+maxNumSolGenes+"\t");

           //Find the gene that, if added to current solution, maximizes (minimizes) the differential coverage
           Vector<String> bestRes=Utility.findBestSol(setGenes,mutationData,mapSampleClasses,currTotalPosCov,currTotalNegCov,
                   currSol.size(),maximization,numPositives,numNegatives,alpha,mapWeights,totalCandWeight);
           String bestGene=bestRes.get(0);
           double bestAvgPosCov=Double.parseDouble(bestRes.get(1));
           double bestAvgNegCov=Double.parseDouble(bestRes.get(2));
           double bestWeight=Double.parseDouble(bestRes.get(3));
           double bestDiffCov=bestAvgPosCov-bestAvgNegCov;
           double posCovPerc=((double)Math.round(bestAvgPosCov*10000))/100;
           double negCovPerc=((double)Math.round(bestAvgNegCov*10000))/100;
           double diffCovPerc=((double)Math.round(bestDiffCov*10000))/100;
           //System.out.println(bestGene);
           //System.out.println(bestDiffCov);
           System.out.println("Best gene: "+bestGene+"\t"+"PosCov: "+posCovPerc+"%\tNegCov: "+negCovPerc+"%\tDiffCov: "+diffCovPerc+"%\tWeight: "+bestWeight);

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
           lastPosCov=posCovPerc;
           lastNegCov=negCovPerc;
           lastDiffCov=diffCovPerc;
           totalWeight+=bestWeight;
           iter++;
       }

       //Print or save results
       fm.writeResults(outputFile,currSol,lastPosCov,lastNegCov,lastDiffCov);

   }

    private static void printHelp()
    {
        String help = "Usage: java -cp ./out DMGSFinder -m <mutationsFile> [-d <listGenes> " +
                "-w <weightsFile> -ms <maxNumSolutionGenes> -a <alpha> -o <resultsFile>]\n\n";
        help+="REQUIRED PARAMETERS:\n";
        help+="-m\tMutation matrix file\n\n";
        help+="OPTIONAL PARAMETERS:\n";
        help+="-d\tList of driver genes\n";
        help+="-w\tFile of gene weights\n";
        help+="-ms\tMaximum number of genes in the reported solution (default=10)\n";
        help+="-a\tAlpha parameter to combine differential coverage and weight scores (default=0.5 or 1.0 if weight are not specified)\n";
        help+="-o\tResults file (default=print results to standard output)\n\n";
        System.out.println(help);
    }
}

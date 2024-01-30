import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Vector;

public class FileManager
{

    public MutationMatrix readMutationDataWithoutClasses(String inputFile, HashSet<String> listRefGenes)
    {
        MutationMatrix mutMatrix=null;
        try
        {
            BufferedReader br=new BufferedReader(new FileReader(inputFile));

            //Read samples
            String str=br.readLine();
            String[] sampleNames=str.split("\t");
            Hashtable<String,String> mapSampleClasses=new Hashtable<>();
            for (String sampleName : sampleNames)
                mapSampleClasses.put(sampleName, "N");
            mutMatrix=new MutationMatrix();

            //Read mutation values
            while((str=br.readLine())!=null)
            {
                String[] split=str.split("\t");
                String gene=split[0];
                if(listRefGenes.contains(gene))
                {
                    for(int i=1;i<split.length;i++)
                    {
                        int freq=Integer.parseInt(split[i]);
                        if(freq>0)
                            mapSampleClasses.put(sampleNames[i-1],"P");
                    }
                }
                else
                {
                    for(int i=1;i<split.length;i++)
                    {
                        int freq=Integer.parseInt(split[i]);
                        String sample=sampleNames[i-1];
                        if(freq>0)
                            mutMatrix.addMutation(gene,sample);
                    }
                }

            }
            br.close();

            //Store sample classes information
            mutMatrix.setMapSampleClasses(mapSampleClasses);
        }
        catch (Exception e){
            System.out.println(e.getMessage());
        }
        return mutMatrix;
    }

    public MutationMatrix readMutationDataWithClasses(String inputFile)
    {
        MutationMatrix mutMatrix=null;
        try
        {
            BufferedReader br=new BufferedReader(new FileReader(inputFile));

            //Read samples
            String str=br.readLine();
            String[] sampleNames=str.split("\t");
            //Read classes
            str=br.readLine();
            String[] listClasses=str.split("\t");
            mutMatrix=new MutationMatrix(sampleNames,listClasses);

            //Read mutation values
            while((str=br.readLine())!=null)
            {
                String[] split=str.split("\t");
                String gene=split[0];
                for(int i=1;i<split.length;i++)
                {
                    int freq=Integer.parseInt(split[i]);
                    String sample=sampleNames[i-1];
                    if(freq>0)
                        mutMatrix.addMutation(gene,sample);
                }
            }
            br.close();
        }
        catch (Exception e){
            System.out.println(e.getMessage());
        }
        return mutMatrix;
    }

    public MutationMatrix readTestData(String inputFile)
    {
        MutationMatrix mutMatrix=null;
        try
        {
            BufferedReader br=new BufferedReader(new FileReader(inputFile));

            //Read samples
            String str=br.readLine();
            String[] sampleNames=str.split("\t");
            Hashtable<String,String> mapSampleClasses=new Hashtable<>();
            for (String sampleName : sampleNames)
                mapSampleClasses.put(sampleName, "-");
            mutMatrix=new MutationMatrix();

            //Read mutation values
            while((str=br.readLine())!=null)
            {
                String[] split=str.split("\t");
                String gene=split[0];
                for(int i=1;i<split.length;i++)
                {
                    int freq=Integer.parseInt(split[i]);
                    String sample=sampleNames[i-1];
                    if(freq>0)
                        mutMatrix.addMutation(gene,sample);
                }

            }
            br.close();

            //Store sample classes information
            mutMatrix.setMapSampleClasses(mapSampleClasses);
        }
        catch (Exception e){
            System.out.println(e.getMessage());
        }
        return mutMatrix;
    }

    public MutationMatrix readMutationDataCV(String inputFile)
    {
        MutationMatrix mutMatrix=null;
        try
        {
            BufferedReader br=new BufferedReader(new FileReader(inputFile));
            mutMatrix=new MutationMatrix();
            Hashtable<String,String> mapSampleClasses=new Hashtable<>();

            //Read genes
            String str=br.readLine();
            String[] split=str.split("\t");
            String[] listGenes=new String[split.length-1];
            System.arraycopy(split, 0, listGenes, 0, split.length - 1);

            //Read mutation values and sample classes
            while((str=br.readLine())!=null)
            {
                split=str.split("\t");
                String sample=split[0];
                int i;
                for(i=1;i<split.length-1;i++)
                {
                    int freq=Integer.parseInt(split[i]);
                    String gene=listGenes[i-1];
                    if(freq>0)
                        mutMatrix.addMutation(gene,sample);
                }
                mapSampleClasses.put(sample,split[i]);
            }
            br.close();
            mutMatrix.setMapSampleClasses(mapSampleClasses);

        }
        catch (Exception e){
            System.out.println(e.getMessage());
        }
        return mutMatrix;
    }

    public HashMap<String,Double> readWeightsFile(String weightsFile)
    {
        HashMap<String,Double> mapWeights=new HashMap<>();
        try
        {
            BufferedReader br=new BufferedReader(new FileReader(weightsFile));
            String str;
            while((str=br.readLine())!=null)
            {
                String[] split=str.split("\t");
                String gene=split[0];
                double weight=Double.parseDouble(split[1]);
                mapWeights.put(gene,weight);
            }
            br.close();
        }
        catch(Exception e) {
            System.out.println(e.getMessage());
        }
        return mapWeights;
    }

    private HashSet<String> readPositiveSetFile(String samplesFile)
    {
        HashSet<String> positiveSet=new HashSet<>();
        try
        {
            String str;
            BufferedReader br=new BufferedReader(new FileReader(samplesFile));
            while((str=br.readLine())!=null)
                positiveSet.add(str);
            br.close();
        }
        catch(Exception e){
            System.out.println(e.getMessage());
        }
        return positiveSet;
    }

    public HashSet<String> getPositiveSet(String inputFile, HashSet<String> driverGenesSet, String positiveFile, String infoType)
    {
        HashSet<String> positiveSet;
        if(infoType.equals("genes"))
        {
            positiveSet=new HashSet<>();
            try
            {
                BufferedReader br=new BufferedReader(new FileReader(inputFile));
                String str=br.readLine();
                String[] setSamples=str.split("\t");
                while((str=br.readLine())!=null)
                {
                    String[] split=str.split("\t");
                    String gene=split[0];
                    String[] split2=gene.split("_");
                    String[] split3=split2[0].split(",");
                    for (String s : split3)
                    {
                        if (driverGenesSet.contains(s))
                        {
                            for (int i = 1; i < split.length; i++)
                            {
                                int freq = Integer.parseInt(split[i]);
                                if (freq > 0)
                                    positiveSet.add(setSamples[i - 1]);
                            }
                        }
                    }
                }
                br.close();
            }
            catch(Exception e) {
                System.out.println(e.getMessage());
            }
        }
        else
            positiveSet=readPositiveSetFile(positiveFile);
        return positiveSet;
    }

    public void writeResults(String outputFile, Vector<String> geneSet, double posCov, double negCov, double diffCov)
    {
        if(outputFile==null)
        {
            System.out.println("\nDMGS: "+geneSet);
            System.out.println("Average positive coverage: "+posCov+"%");
            System.out.println("Average negative coverage: "+negCov+"%");
            System.out.println("Differential coverage: "+diffCov+"%");
        }
        else
        {
            try{
                BufferedWriter bw=new BufferedWriter(new FileWriter(outputFile));
                bw.write("DMGS\tAverage positive coverage\tAverage negative coverage\tDifferential coverage\n");
                bw.write(geneSet+"\t"+posCov+"%\t"+negCov+"%\t"+diffCov+"%\n");
                bw.close();
            }
            catch (Exception e){
                System.out.println(e.getMessage());
            }
            System.out.println("Results written to "+outputFile);
        }
    }

    public void writePredictions(Hashtable<String,String> realClasses, Hashtable<String,String> predictedClasses, double runningTime, String outputFile)
    {
        try {
            BufferedWriter bw = new BufferedWriter(new FileWriter(outputFile));
            bw.write("real\tpredicted\n");
            for(String sample : realClasses.keySet())
                bw.write(realClasses.get(sample)+"\t"+predictedClasses.get(sample)+"\n");
            bw.write("Time\t"+runningTime+"\n");
            bw.close();
        }
        catch (Exception e){
            System.out.println(e.getMessage());
        }
    }

    public void writeMUTClassResults(String resultsFile, Vector<String> panelPosGenes, Vector<String> panelNegGenes, Hashtable<String,String> predictedClasses)
    {
        try{
            BufferedWriter bw=new BufferedWriter(new FileWriter(resultsFile));
            bw.write("Positive panel:\n");
            bw.write(panelPosGenes+"\n\n");
            bw.write("Negative panel:\n");
            bw.write(panelNegGenes+"\n\n");
            bw.write("Predicted classes:\n");
            for(String sample : predictedClasses.keySet())
                bw.write(sample+"\t"+predictedClasses.get(sample)+"\n");
            bw.close();
        }
        catch(Exception e){
            System.out.println(e.getMessage());
        }
    }
}

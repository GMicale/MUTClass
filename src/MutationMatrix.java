import java.util.HashMap;
import java.util.HashSet;
import java.util.Hashtable;

public class MutationMatrix
{
    private Hashtable<String,String> mapSampleClasses;
    private final HashMap<String,HashSet<String>> mutationData;

    public MutationMatrix()
    {
        mapSampleClasses=new Hashtable<>();
        this.mutationData=new HashMap<>();
    }

    public MutationMatrix(String[] listSamples, String[] sampleClasses)
    {
        mapSampleClasses=new Hashtable<>();
        for(int i=0;i< listSamples.length;i++)
            mapSampleClasses.put(listSamples[i],sampleClasses[i]);
        this.mutationData=new HashMap<>();
    }

    public void setMapSampleClasses(Hashtable<String,String> mapSampleClasses)
    {
        this.mapSampleClasses=mapSampleClasses;
    }

    public HashSet<String> getPositiveSet()
    {
        HashSet<String> positiveSet=new HashSet<>();
        for(String sample : mapSampleClasses.keySet())
        {
            if(mapSampleClasses.get(sample).equals("P"))
                positiveSet.add(sample);
        }
        return positiveSet;
    }

    public HashSet<String> getNegativeSet()
    {
        HashSet<String> negativeSet=new HashSet<>();
        for(String sample : mapSampleClasses.keySet())
        {
            if(mapSampleClasses.get(sample).equals("N"))
                negativeSet.add(sample);
        }
        return negativeSet;
    }

    public Hashtable<String,String> getMapSampleClasses()
    {
        return mapSampleClasses;
    }

    public HashMap<String,HashSet<String>> getMutationData()
    {
        return mutationData;
    }

    public HashSet<String> getAllGenes()
    {
        return new HashSet<>(mutationData.keySet());
    }

    public void addMutation(String gene, String sample)
    {
        if(mutationData.containsKey(gene))
            mutationData.get(gene).add(sample);
        else
        {
            HashSet<String> setSamples=new HashSet<>();
            setSamples.add(sample);
            mutationData.put(gene,setSamples);
        }
    }

}

import be.ac.ulg.montefiore.run.jahmm.Hmm;
import be.ac.ulg.montefiore.run.jahmm.Observation;
import be.ac.ulg.montefiore.run.jahmm.Opdf;
import br.ufpr.bioinfo.genbank.Entry;
import br.ufpr.bioinfo.genbank.Feature;
import br.ufpr.bioinfo.genbank.Sequence;
import br.ufpr.bioinfo.genbank.io.GenBankReader;
import br.ufpr.bioinfo.genbank.parser.GenBankParser;

import java.io.IOException;
import java.lang.String;
import java.util.ArrayList;
import java.util.List;



public class main {

		
	public static void main(String[] args) throws IOException {
		
		//buffering *.gb file
		String buffer = GenBankReader.readFile("Data/coli.gb");
		//instanciating GenPank Parsing class
		GenBankParser gp = new GenBankParser();
		//get entry
		Entry entry = gp.process(buffer);
		System.out.println(entry);
		//get all the features of the entry
		List<Feature>features = entry.getFeatures();
		//extract Sequence
		Sequence seq = entry.getSequence();
		
		//initialize List for all Gene positions
		List<int[]> genePositions = new ArrayList<int[]>();
		
		//for all features
		for(int i = 1; i < features.size(); i++){
			
		Feature f = features.get(i);
		//if feature is not on the complementary strand and is a gene add its position
		if(!f.getLocation().isComplement()){ 
			if(f.getKey().equals("gene")){
					//getting positions
					int[] positions = {f.getLocation().getMinor(),f.getLocation().getMajor()};
					genePositions.add(positions);
		
			}
			}
		}
		//success
		System.out.println("Gene positions extracted");
		buildHmm(seq, genePositions);
		

	}
	
	public static void buildHmm(Sequence seq, List<int[]> genePositions) throws IOException  
	{
		//Probability of starting in Non-Gene vs Gene. The assumptions is we always start in Non-Gene
		double[] pi = {1,0};
		//transistion probabilitys of going from state i to state j
		double[][] a = calcTransitionProbs(seq, genePositions);
		//List of Observation distributions
		List<Opdf> opdfs;
		int[][] counts;
	
	
		
		

		
	}
	
	public static double[][] calcTransitionProbs(Sequence seq, List<int[]> genePositions)
	{
		
		String seqString = seq.getSequence();
		//initialize transistion probability matrix 
		double[][] transProbs = new double[2][2];
		//look at all transitions
		for(int i = 0; i < seqString.length()-1; i++ )
		{
			if( isInAGene(i,genePositions) )
			{
				if(isInAGene(i+1,genePositions))
				{
					transProbs[1][1]++;
				}
				else
				{
					transProbs[1][0]++;
				}
			}
			else
				{
				if(isInAGene(i+1,genePositions))
				{
					transProbs[0][1]++;
				}
				else
				{
					transProbs[0][0]++;
				}
				}
						
				
		}
		//normalize over all transitions in one state
		
		double[] count = {(transProbs[0][0] + transProbs[0][1]) ,(transProbs[1][0] + transProbs[1][1])};
		for(int i = 0; i < transProbs.length; i ++)	
		{
			for(int j = 0; j < transProbs[0].length; j++)
			{

				transProbs[i][j] = transProbs[i][j]/count[i];
				System.out.println("from " + i + " to " + j);
				System.out.println(transProbs[i][j]);
			}
		}
		
		return transProbs;
	}
	
	
	
	
	//checks if a position is the given list of Gene Positions
	public static boolean isInAGene(int pos, List<int[]> genePositions)
	{
		//interates over all positions pairs
		for(int i = 0; i < genePositions.size(); i++)
		{
			int[] checkPos = genePositions.get(i);
			//if position is in range, then true
			if(pos >= checkPos[0] & pos <= checkPos[1] )
			{
				return true;
			}
		}
		return false;
	}

}

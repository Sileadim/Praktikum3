import be.ac.ulg.montefiore.run.jahmm.Hmm;
import be.ac.ulg.montefiore.run.jahmm.Observation;
import be.ac.ulg.montefiore.run.jahmm.ObservationInteger;
import be.ac.ulg.montefiore.run.jahmm.Opdf;
import be.ac.ulg.montefiore.run.jahmm.OpdfInteger;
import be.ac.ulg.montefiore.run.jahmm.learn.BaumWelchLearner;
import be.ac.ulg.montefiore.run.jahmm.toolbox.MarkovGenerator;
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
		System.out.println(baumWelch(buildHmm(seq, genePositions),100).toString());
		
//		List<int[]>test = new ArrayList<int[]>();
//		int[]test1 = {3,8};
//		int[]test2 = {12,19};
//		test.add(test1);
//		test.add(test2);		
//		System.out.println(liststring(threeStates(test, 22)));

	}
	public static String liststring(List<List<int[]>> list)
	{
		String bla = "[";
		for(int i = 0; i < (list.size()); i++ )
		{
			bla = bla+"[";
			for(int j = 0; j < (list.get(i).size()); j++ )
			{
				bla = bla+"[";
				for(int k = 0; k < (list.get(i).get(j).length); k++ )
				{
					bla = bla+ (list.get(i).get(j)[k]) + " ";
				}
				bla = bla+"]";
			}
			bla = bla+"]";
		}
		bla = bla+"]";
		return(bla);
	}

	public static List<List<int[]>> twoStates(List<int[]> genelist, int l) // l: length of seq
	{
		List<int[]> secondState = new ArrayList<int[]>();
		
		if (genelist.get(0)[0] != 0)
		{
			int[] firstpair = {0, genelist.get(0)[0] - 1};
			secondState.add(firstpair);
		}
					
		for(int i = 0; i < (genelist.size() - 1); i++ )
		{
			int anf = (genelist.get(i)[1]+1);
			int end = (genelist.get(i+1)[0]-1);
			int[] pair = {anf,end};
			secondState.add(pair);
		}
		
		if (genelist.get(genelist.size()-1)[1] != l-1)
		{
			int[] lastpair = {genelist.get(genelist.size()-1)[1],l-1};
			secondState.add(lastpair);
		}
		
		List<List<int[]>> ret = new ArrayList<List<int[]>>();
		ret.add(genelist);
		ret.add(secondState);
		return(ret);
	}	
	
	public static List<List<int[]>> threeStates(List<int[]> genelist, int l) // l: length of seq
	{
		List<List<int[]>> twostates = twoStates(genelist,l);
		List<int[]> thirdState = new ArrayList<int[]>();
		for(int i = 0; i < (genelist.size()); i++ )
		{
			int[] codonpair = {twostates.get(0).get(i) [0], twostates.get(0).get(i) [0] + 2};
			List<int[]> clonepairs = twostates.get(0);
			int[] pair = clonepairs.get(i);
			pair[0] = pair[0] + 3;
			clonepairs.set(i, pair);
			twostates.set(0, clonepairs);
			

			thirdState.add(codonpair);
		}
		twostates.add(thirdState);
		return(twostates);
	}
	/*
	public static List<List<int[]>> fiveStates(List<int[]> genelist, int l) // l: length of seq
	{
		List<List<int[]>> threestates = threeStates(genelist,l);
		List<int[]> thirdState = threestates.get(2);
		for(int i = 0; i < (genelist.size()); i++ )
		{
			int[] codonpair = {twostates.get(0).get(i) [0], twostates.get(0).get(i) [0] + 2};
			List<int[]> clonepairs = twostates.get(0);
			int[] pair = clonepairs.get(i);
			pair[0] = pair[0] + 3;
			clonepairs.set(i, pair);
			twostates.set(0, clonepairs);
			

			thirdState.add(codonpair);
		}
		twostates.add(thirdState);
		return(twostates);
	}
	*/
	public static Hmm<ObservationInteger> buildHmm(Sequence seq, List<int[]> genePositions) throws IOException  
	{
		//Probability of starting in Non-Gene vs Gene. The assumptions is we always start in Non-Gene
		double[] pi = {0.33,0.34,0.33};
		//transistion probabilitys of going from state i to state j
		
		List<List<int[]>> two= twoStates(genePositions, seq.getSequence().length());
		List<List<int[]>> three = threeStates(genePositions, seq.getSequence().length());
		double[][] a = calcTransitionProbs(seq, three);
		//List of Observation distributions
		List<Opdf> opdfs = calcEmissionProbs(seq, three );
		Hmm<ObservationInteger> hmm = new Hmm(pi,a, opdfs);
		//System.out.println(hmm.toString());
		
		
		return hmm;	
		
		
	}
	
	
	public static  double[] countEmissionProbs(Sequence seq, List<int[]> stateList )
	{
		String seqString = seq.getSequence();
		double[] emProbs = new double[4];
		for(int i = 0; i < seqString.length(); i++ )
		{
		if(seqString.charAt(i) == 'a')
		{
			emProbs[0]++;
		}
		if(seqString.charAt(i) == 'c')
		{
			emProbs[1]++;
		}

		if(seqString.charAt(i) == 'g')
		{
			emProbs[2]++;
		}
		else
		{
			emProbs[3]++;
		}
		}
		
		double sum = emProbs[0] + emProbs[1] + emProbs[2] + emProbs[3];
		for(int i = 0; i < emProbs.length ; i++)
		{
			emProbs[i] = emProbs[i] /sum;
		}
		return emProbs;
		
		
	}
	public static List<Opdf> calcEmissionProbs(Sequence seq, List<List<int[]>> featurePositions)
	{
		List<Opdf> opdfs = new ArrayList<Opdf>();
		String seqString = seq.getSequence();
				
		//look at all characters
		for(int i = 0; i < featurePositions.size(); i++ ){
			Opdf state = new OpdfInteger(countEmissionProbs(seq, featurePositions.get(i)));
			opdfs.add(state);
			
		}
	
		return opdfs;

				
	}
	
	
	
	
	public static int isInState(int x,  List<List<int[]>> featurePositions)
	{
		int ftNum = featurePositions.size();
		for(int i = 0; i < ftNum; i++)
		{		
			for(int j = 0; j < featurePositions.get(i).size(); j++ )
			{
			
				int [] tuple = featurePositions.get(i).get(j);
					if(x >= tuple[0] & x <= tuple[1])
						{
						return i;
						}
	
		}
		
	}
		return -1;
	
	}
	
	public static double arraySum (double [] array)
	{
		double sum = 0;
		for(int i = 0; i < array.length; i++)
		{
			sum += array[i];
		}
		return sum;
	}
	
	public static double[][] calcTransitionProbs(Sequence seq, List<List<int[]>> featurePositions)
	{
		int ftNum = featurePositions.size();
		String seqString = seq.getSequence();
		//initialize transition probability matrix 
		double[][] transProbs = new double[ftNum][ftNum];
		//look at all transitions
		int [] states = {isInState(0, featurePositions), isInState(1, featurePositions)};
		for(int i = 0; i < seqString.length()-1; i++ )
		{
			if( states[0] > -1 & states[1] > -1)
				{
				
				transProbs[states[0]][states[1]]++;
				states[0] = states[1];
				states[1] = isInState(i+1, featurePositions);
				
				
				}
			
			
			else
			{
				System.out.println("Index not in Positions!");
			}
		}
		//normalize over all transitions in one state
		
		double sum = seqString.length() -1;
		for(int i = 0; i < transProbs.length; i ++)	
		{
			double arraySum = arraySum(transProbs[i]);
			for(int j = 0; j < transProbs[0].length; j++)
			{

				transProbs[i][j] = transProbs[i][j]/arraySum;
		}
		}
		
		return transProbs;
	}
	
	
	
	
	/*checks if a position is the given list of Gene Positions
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
*/
	
	//execute the Baum-Welch-algorithm for a given hmm
	public static <O extends Observation>  Hmm<O> baumWelch(Hmm<O> hmm, int nbOfIteration)
	{
		//create a MarkovGenerator object 
		MarkovGenerator<O> mg = new MarkovGenerator<O>(hmm);
		//List<List<O>> object for the sequences
		List<List<O>> sequences = new ArrayList<List<O>>();
		
		//build 200 sequences of length 500 based on the given hmm
		for (int i = 0; i < 200; i++)
			sequences.add(mg.observationSequence(500));
			System.out.println(mg.observationSequence(500));;
		
		//initialize a Baum-Welch-object
		BaumWelchLearner bwl = new BaumWelchLearner();
		//set numbers of iterations
		bwl.setNbIterations(nbOfIteration);
		//learn new hmm from old hmm and builded sequences
		System.out.println("execute Baum-Welch algorithm");
		Hmm<O> learntHmm = bwl.learn(hmm, sequences);
		return learntHmm;
	}
	
}
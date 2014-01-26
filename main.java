import be.ac.ulg.montefiore.run.jahmm.Hmm;
import be.ac.ulg.montefiore.run.jahmm.Observation;
import be.ac.ulg.montefiore.run.jahmm.ObservationInteger;
import be.ac.ulg.montefiore.run.jahmm.Opdf;
import be.ac.ulg.montefiore.run.jahmm.OpdfInteger;
import be.ac.ulg.montefiore.run.jahmm.OpdfIntegerFactory;
import be.ac.ulg.montefiore.run.jahmm.ViterbiCalculator;
import be.ac.ulg.montefiore.run.jahmm.io.FileFormatException;
import be.ac.ulg.montefiore.run.jahmm.io.ObservationIntegerReader;
import be.ac.ulg.montefiore.run.jahmm.io.ObservationSequencesReader;
import be.ac.ulg.montefiore.run.jahmm.learn.BaumWelchLearner;
import be.ac.ulg.montefiore.run.jahmm.toolbox.MarkovGenerator;
import br.ufpr.bioinfo.genbank.Entry;
import br.ufpr.bioinfo.genbank.Feature;
import br.ufpr.bioinfo.genbank.Sequence;
import br.ufpr.bioinfo.genbank.io.GenBankReader;
import br.ufpr.bioinfo.genbank.parser.GenBankParser;

import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.Reader;
import java.io.StringReader;
import java.lang.String;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;



public class main {

		
	public static void main(String[] args) throws IOException, FileFormatException {
		
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

		
		
		/*
		List<List<ObservationInteger>> genes = new ArrayList<List<ObservationInteger>>();
		List<List<ObservationInteger>> nongenes = new ArrayList<List<ObservationInteger>>();
		
		List<ObservationInteger> seqObs = getSeqObs(seq.getSequence());
		
		System.out.println("start gene search");
		genes = getGeneSubSequences(seq, features, seqObs);
		System.out.println("start non-gene search");
		nongenes = getNonGeneSubSequences(seq, features, seqObs);
		System.out.println("done");
		*/
		
		
		//if feature is not on the complementary strand and is a gene add its position
		
		
		//success
		System.out.println("Gene positions extracted");
		//Hmm<ObservationInteger> hmm = buildHmm(seq, features);
		//System.out.println(hmm.toString());
		//leaveOneOut(features, 4 , seq);
		Hmm < ObservationInteger > dummyHmm = new Hmm < ObservationInteger > (2 , new OpdfIntegerFactory(4) );
		System.out.println(dummyHmm.toString());
		dummyHmm = baumWelch(dummyHmm ,1, seq);
		
		System.out.println(dummyHmm.toString());
		
		/*
		int[] states = getMostProbableStates(getSeqObs(seq.getSequence()), dummyHmm);
		for(int i = 0; i < states.length; i++){
			System.out.println(states[i]);
		}
*/

	}
	
	//given a list of obeservations and a hmm calculate the most probable sequence of states
	public static int[] getMostProbableStates(ArrayList<ObservationInteger> obs, Hmm<ObservationInteger> hmm)
	{
	ViterbiCalculator viterbi = new ViterbiCalculator(obs, hmm);
	return viterbi.stateSequence();
		
		
	}
	//given to strings returns the percentage of matching letters
	public static double checkPercentageIdentity(int[] a, int[] b)
	{
		double count = 0;
		if(a.length != b.length)
		{
			System.out.println("Lengths of strings do not match!");
			return -1;
		}
		for(int i = 0; i < a.length; i++)
		{
			if(a[i] == b[i])
			{
				count++;
			}
		}
		return (count/a.length);
		
	}
	//return seqObs list from  given Sequence
	public static  ArrayList<ObservationInteger> getSeqObs(String seqString){

		ArrayList<ObservationInteger> seqObs = new ArrayList<ObservationInteger>();
	
	for (int i=0; i<seqString.length();i++){
		if(seqString.charAt(i) == 'a'){
			ObservationInteger help = new ObservationInteger(0);
			seqObs.add(help);
		}
		else if(seqString.charAt(i) == 'c'){
			ObservationInteger help = new ObservationInteger(1);
			seqObs.add(help);
		}
		else if(seqString.charAt(i) == 'g'){
			ObservationInteger help = new ObservationInteger(2);
			seqObs.add(help);
		}
		else {
			ObservationInteger help = new ObservationInteger(3);
			seqObs.add(help);
		}
	}
	return seqObs;	
	}
	
	public static <O extends Observation> ArrayList<O> getSeqObs2(String seqString){

		ArrayList<O> seqObs = new ArrayList<O>();
	
	for (int i=0; i<seqString.length();i++){
		if(seqString.charAt(i) == 'a'){
			ObservationInteger help = new ObservationInteger(0);
			seqObs.add((O) help);
		}
		else if(seqString.charAt(i) == 'c'){
			ObservationInteger help = new ObservationInteger(1);
			seqObs.add((O) help);
		}
		else if(seqString.charAt(i) == 'g'){
			ObservationInteger help = new ObservationInteger(2);
			seqObs.add((O) help);
		}
		else {
			ObservationInteger help = new ObservationInteger(3);
			seqObs.add((O) help);
		}
	}
	return seqObs;	
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

	//given a list of the gene positions also calculate the positions of the non genes
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
	//given a list of gene positions also calculate the positions of non genes and the start codons
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
	
	// for a Sequence and a number of Subsequences train a Hmm with all but the Subsequence and check the Percentage of correctly asserted states
	public static void leaveOneOut(List<Feature> features, int n, Sequence seq ) throws IOException
	
	
	{
		System.out.println("starting leave one out validation");
		
		System.out.println("getting gene positions");
		List<int[]> genePos = getGenePositions(features);
		System.out.println("getting two states positions");
		List<List<int[]>> twoPos= twoStates(genePos, seq.getSequence().length());
		//List<List<int[]>> twoPos= threeStates(genePos, seq.getSequence().length());

		double[] pi = {0.5, 0.5};
		System.out.println("Cutting into subsequences");
		Map<Integer, String> map= cutToSubsequences(seq, n);
		Object[] keys = map.keySet().toArray();
		int subsequenceLength = (int) (Math.floor(seq.getSequence().length()/n));
		System.out.println("Beginning model evaluation");
		for(int i = 0; i < keys.length; i++)
		{
			int key = (Integer) keys[i];
			int[] notTraining = {key, key + subsequenceLength };
			System.out.println("training hmm with subsequence");
			Hmm<ObservationInteger> hmm = buildHmm(seq, twoPos, pi, notTraining );
			String subSequence = map.get(key);
			double identity = checkPercentageIdentity(stringToStates(twoPos, subSequence , key), getMostProbableStates(getSeqObs(subSequence), hmm));
			System.out.println("Percentage of correctly asserted states");
			System.out.println(identity * 100);
		}
		
		
		
		
		
		
	}
	
	//cuts a Sequence in n subsequences and returns them with the absolute positions in the Sequence as keys
	public static Map<Integer, String> cutToSubsequences(Sequence seq, int n)
	{
		
		Map<Integer, String> map = new HashMap<Integer, String>();
		
		int subsequenceLength = (int) (Math.floor(seq.getSequence().length()/n));
		for(int i = 0; i < seq.getSequence().length() - 10;)
		{	
			if(i+subsequenceLength > seq.getSequence().length() )
			{
			map.put(i, (seq.getSubSequence(i, i+subsequenceLength)));
			}
			else
			{
				map.put(i, (seq.getSubSequence(i, seq.getSequence().length() )));
			}
			i += subsequenceLength;
		}
		return map;
			
		
	}	
	
	//given the int[] of all States, subsequence and the absolute position of the subsequence in the whole sequence
	//calculate the int[] of the states
	public static int[] stringToStates(List<List <int[]>> positions, String sequence, int absPos)
	{
		int[] stateSequence = new int[sequence.length()];
		
			
			for(int j = 0; j < sequence.length(); j++)
			{
				
					stateSequence[j] = isInState(j + absPos, positions);
				
			}
			
			
	return stateSequence;		
	}
			
		
		
// get the gene positions from a list of the positions of all features
	public static List<int[]> getGenePositions(List<Feature> features)
	{
		List<int[]> genePositions = new ArrayList<int[]>();
		//for all features
		for(int i = 1; i < features.size(); i++){
					
			Feature f = features.get(i);
		//initialize List for all Gene positions
		
		if(!f.getLocation().isComplement()){ 
			if(f.getKey().equals("gene")){
					//getting positions
					int[] positions = {f.getLocation().getMinor(),f.getLocation().getMajor()};
					genePositions.add(positions);
		
			}
			}
		}	
		return genePositions;
		
	}
	//builds a Hmm given the Seqeunce, the positions of the states, the intial probabilities and the positions where the Hmm should not be trained
	public static Hmm<ObservationInteger> buildHmm(Sequence seq, List<List<int[]>> positions, double[] pi, int[] notIn ) throws IOException  
	{	
		
		
		//transistion probabilitys of going from state i to state j
		
		//List<List<int[]>> three = threeStates(genePositions, seq.getSequence().length());
		double[][] a = calcTransitionProbs(seq, positions, notIn);
		//List of Observation distributions
		List<Opdf> opdfs = calcEmissionProbs(seq, positions, notIn );
		
		Hmm<ObservationInteger> hmm = new Hmm(pi,a, opdfs);
		//System.out.println(hmm.toString());
		
		
		return hmm;	
		
		
	}
	public static  double[] countEmissionProbs(Sequence seq, List<int[]> stateList, int[] notIn )
	{
		String seqString = seq.getSequence();
		double[] emProbs = new double[4];
		for(int i = 0; i < seqString.length(); i++ )
		{
		if(i < notIn[0] & i > notIn[1] )	
		{
			if(isInAState(i, stateList) )
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
		
		}
	}
		double sum = emProbs[0] + emProbs[1] + emProbs[2] + emProbs[3];
		for(int i = 0; i < emProbs.length ; i++)
		{
			emProbs[i] = emProbs[i] /sum;
		}
		return emProbs;
		
		
	}
	
	
	public static List<Opdf> calcEmissionProbs(Sequence seq, List<List<int[]>> featurePositions, int[] notIn)
	{
		List<Opdf> opdfs = new ArrayList<Opdf>();
		String seqString = seq.getSequence();
				
		//look at all characters
		for(int i = 0; i < featurePositions.size(); i++ ){
			Opdf state = new OpdfInteger(countEmissionProbs(seq, featurePositions.get(i), notIn));
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
	
	public static double[][] calcTransitionProbs(Sequence seq, List<List<int[]>> featurePositions, int[] notIn)
	{
		int ftNum = featurePositions.size();
		String seqString = seq.getSequence();
		//initialize transition probability matrix 
		double[][] transProbs = new double[ftNum][ftNum];
		//look at all transitions
		int [] states = {isInState(0, featurePositions), isInState(1, featurePositions)};
		for(int i = 0; i < seqString.length()-1; i++ )
		{
			if(i < notIn[0] & i > notIn[1])
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

	//checks if a position is the given list of Gene Positions
	public static boolean isInAState(int pos, List<int[]> statePositions)
	{
		//interates over all positions pairs
		for(int i = 0; i < statePositions.size(); i++)
		{
			int[] checkPos = statePositions.get(i);
			//if position is in range, then true
			if(pos >= checkPos[0] & pos <= checkPos[1] )
			{
				return true;
			}
		}
		return false;
	}
	
	
	//function return all gene subsequences of a given sequence
	public static <O extends Observation> List<List<ObservationInteger>> getGeneSubSequences(Sequence seq, List<Feature> features, List<ObservationInteger> seqObs) throws IOException, FileFormatException{
		
		List<List<ObservationInteger>> genes = new ArrayList<List<ObservationInteger>>();
				
		for(int i = 1; i < features.size(); i++){
			
			Feature f = features.get(i);
			if(!f.getLocation().isComplement()){
				if(f.getKey().equals("gene")){
					int geneBegin = f.getLocation().getMinor();
					int geneEnd = f.getLocation().getMajor();
					
					genes.add(seqObs.subList(geneBegin, geneEnd));
				}							
			}		
		}
		
		return genes;
	}
	
	
	//function return all non-gene subsequences of a given sequence
	public static <O extends Observation> List<List<ObservationInteger>> getNonGeneSubSequences(Sequence seq, List<Feature> features, List<ObservationInteger> seqObs) throws IOException, FileFormatException{
		
		List<List<ObservationInteger>> non_genes = new ArrayList<List<ObservationInteger>>();
				
		for(int i = 1; i < features.size(); i++){
			
			Feature f = features.get(i);
			if(!f.getLocation().isComplement()){
				if(!f.getKey().equals("gene")){
					int nonGeneBegin = f.getLocation().getMinor();
					int nonGeneEnd = f.getLocation().getMajor();
					non_genes.add(seqObs.subList(nonGeneBegin, nonGeneEnd));
				}			
			}			
		}
		
		return non_genes;
	}
	
	
	//execute the Baum-Welch-algorithm for a given hmm
		public static <O extends Observation>  Hmm<O> baumWelch(Hmm<O> hmm, int nbOfIteration, Sequence seq)
		{
			List<O> seqObs = getSeqObs2(seq.getSequence());
			
			
			//create a MarkovGenerator object 
			MarkovGenerator<O> mg = new MarkovGenerator<O>(hmm);
			//List<List<O>> object for the sequences
			List<List<O>> sequences = new ArrayList<List<O>>();
				
			//build 200 sequences of length 500 based on the given hmm
			for (int i = 0; i < 200; i++)
			{
				List<O> list = mg.observationSequence(500);
				
					System.out.print(list.toString());
				
				
				sequences.add(list);
				System.out.println(list.getClass());

			}
			System.out.println("\n");
			
			List<List<O>> sequences2 = new ArrayList<List<O>>();
			
			//sequences.add((List<O>)seqObs);
			int count = (int)Math.ceil(seqObs.size()/500);
				
			for(int i = 0 ; i < count ; i++)
			{
				if (i == count-1)
				{	
					
					List<O> list = new ArrayList<O>(seqObs.subList(i*500, seqObs.size()-1));
					
						System.out.print(list.toString());
						System.out.println(list.getClass());
					
					sequences2.add( (seqObs.subList(i*500, seqObs.size()-1)));
				}
				else
				{
					List<O> list = new ArrayList<O>(seqObs.subList(i*500, (i+1)*500));
					sequences2.add(list);
				}
				
			}
			
			//initialize a Baum-Welch-object

			BaumWelchLearner bwl = new BaumWelchLearner();
			//set numbers of iterations
			bwl.setNbIterations(nbOfIteration);
			//learn new hmm from old hmm and builded sequences
			System.out.println("execute Baum-Welch algorithm");
			Hmm<O> learntHmm = bwl.learn(hmm, sequences2);
			return learntHmm;
		}
	
}
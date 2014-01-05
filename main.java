import br.ufpr.bioinfo.genbank.Entry;
import br.ufpr.bioinfo.genbank.Feature;
import br.ufpr.bioinfo.genbank.io.GenBankReader;
import br.ufpr.bioinfo.genbank.parser.GenBankParser;

import java.io.IOException;
import java.lang.String;
import java.util.List;



public class main {

		
	public static void main(String[] args) throws IOException {
		readGene();

	}
	
	public static void readGene() throws IOException  
	{
		
		String buffer = GenBankReader.readFile("Data/coli.gb");
		GenBankParser gp = new GenBankParser();
		Entry entry = gp.process(buffer);
		System.out.println(entry);
		List<Feature>features = entry.getFeatures();
		for(int i = 1; i < features.size(); i++){
			
		Feature f = features.get(i);
		String featurename = f.getKey();
		if(!f.getLocation().isComplement()){ 
			if(featurename.equals("gene")){
		
		System.out.println(featurename);

		System.out.println(f.getLocation().getMinor());		
		System.out.println(f.getLocation().getMajor());	
			}
			}
		}

		
	}

}

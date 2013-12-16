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
		for(int i =0; i < features.size(); i++){
		System.out.println(features.get(i).getKey());	
		}

		
	}

}

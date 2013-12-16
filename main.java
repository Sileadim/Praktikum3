import br.ufpr.bioinfo.genbank.Entry;
import br.ufpr.bioinfo.genbank.io.GenBankReader;
import br.ufpr.bioinfo.genbank.parser.GenBankParser;

import java.lang.String;



public class main {

		
	public static void main(String[] args) {
		readGene();

	}
	
	public static void readGene()  
	{
		
		String buffer = GenBankReader.readFile("/Data/coli.gb");
		GenBankParser gp = new GenBankParser();
		Entry entry = gp.process(buffer);
		System.out.println("File read");
	}

}

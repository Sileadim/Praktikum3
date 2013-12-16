import be.ac.ulg.montefiore.run.jahmm.*;
import br.ufpr.bioinfo.genbank.parser.*;
import org.biojava.*;
import java.*;


public void readGenes ()  
{
	
	String buffer = GenBankReader.readFile("/Data/coli.gb");
	GenBankParser gp = new GenBankParser();
	Entry entry = gp.process(buffer);
	System.out.println("File read")
}


public static void main()
{
	readGenes();
}






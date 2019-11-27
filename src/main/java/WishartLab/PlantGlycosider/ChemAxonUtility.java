package WishartLab.PlantGlycosider;

import java.io.IOException;

public class ChemAxonUtility {
	
	
	
//	public static int calculateLogD() {
//		try {
//		    // instantiate the plugin objects
//		    logDPlugin logDPlugin = new logDPlugin(); 
//		    IUPACNamingPlugin iupacNamingPlugin = new IUPACNamingPlugin();
//
//		    // set the parameters for the calculations
//		    // MajorMicrospeciesPlugin parameters
//		    mmsPlugin.setpH(7.4); // major microspecies generation at pH = 7.4
//		    // TPSAPlugin parameters
//		    tpsaPlugin.setpH(7.4); // surface area calculation at pH = 7.4
//		    // logDPlugin parameters
//		    // set the Cl- and Na+/K+ concentration
//		    logDPlugin.setCloridIonConcentration(0.15);
//		    logDPlugin.setNaKIonConcentration(0.15);
//		    // set the pH range and pH step size
//		    logDPlugin.setpHLower(5.4);
//		    logDPlugin.setpHUpper(9.4);
//		    logDPlugin.setpHStep(2.0);
//		    
//		    MolImporter importer = new MolImporter(args[0]);
//		    Molecule mol;
//		    logDPlugin.setMolecule(mol);
//		    logDPlugin.run();
//		    
//		    double[] pHs = logDPlugin.getpHs();
//		    double[] logDs = logDPlugin.getlogDs();
//		    
//		    // store the pH - logD pairs in a string in format
//	        // "[pH1]:logD1 [pH2]:logD2 ..."
//		    String logDresult = "";
//	        for (int i = 0; i < pHs.length; i++) {
//			    	double pH = pHs[i];
//			    	double logD = logDs[i];
//			    	logD = Math.rint(logD * 100)/100; // round
//			    	logDresult += ("[" + pH + "]:" + logD + " ");
//			}
//			importer.close();
//			
//		} catch (IOException e) {
//		    System.err.println("I/O error has occurred.");
//		    e.printStackTrace();
//		} catch (PluginException e) {
//		    System.err.println("Plugin processing or calculation error.");
//		    e.printStackTrace();
//		}
//	    }
//	}
	/**
	 * molecule.toFormat("name") generates the "preferred IUPAC name". 
	 * You can use molecule.toFormat("name:t") to generate the traditional name (same syntax name:t with MolConverter).
	 */
//	public static String getName() {
//		Molecule m = ...;
//		String name = m.toFormat("name");
//		String traditionalName = m.toFormat("name");
// 		String IUPACName = m.toFormat(
//	}


}

/**
 * 
 * @author Djoumbou Feunang, Yannick, PhD
 *
 */

package WishartLab.PlantGlycosider;



/**
 * A class that implements functions to analyze the structure of a molecule
 * (e.g. through structure search). It also implements functions on collections
 * of molecules, such as the removal of duplicates.
 */

/**
 * @author Yannick Djoumbou
 *
 */

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;


import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.graph.ConnectivityChecker;
import org.openscience.cdk.inchi.InChIGenerator;
import org.openscience.cdk.inchi.InChIGeneratorFactory;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.qsar.descriptors.molecular.ALOGPDescriptor;
import org.openscience.cdk.qsar.result.IDescriptorResult;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.cdk.tools.manipulator.MolecularFormulaManipulator;

import ambit2.smarts.query.SMARTSException;
import ambit2.smarts.query.SmartsPatternCDK;


// import ambit2.structure2name.IUPACNameGenerator;


public class ChemStructureExplorer {

	protected static IChemObjectBuilder 	builder 		= SilentChemObjectBuilder.getInstance();
	protected static SmilesParser	    smiParser		= new SmilesParser(builder);
	public static SmilesGenerator 		smiGen			= SmilesGenerator.isomeric();
	public static InChIGeneratorFactory inchiGenFactory;
	
	
	/**
	 * Biodegradation Not Predicted rules based on http://eawag-bbd.ethz.ch/predict/notbepredicted.html
	 * @param molecule
	 * @return
	 * @throws CDKException
	 */
	public static boolean isBtValid(IAtomContainer molecule) throws CDKException{
		return AtomContainerManipulator.getNaturalExactMass(molecule) < 1000.0 && containsCarbon(molecule) && !isMixture(molecule) ;
	}	
	
	/**
	 * contain the validation with isBtValid (i.e. more stricted rules than isBtValid(IAtomContainer molecule)
	 * @param molecule
	 * @return
	 * @throws CDKException
	 */
	public static boolean isPPSValid(IAtomContainer molecule) throws CDKException {
		return isPpsValid(molecule);
	}
	
	
	/**
	 * The pure secondary metabolite shouldn't include any glucoside, sulfate and glycinate compound/attachment
	 * @param molecule
	 * @return
	 * @throws SMARTSException
	 * @throws CDKException
	 * @throws CloneNotSupportedException
	 */
	public static boolean isPureSecondaryMetoblite(IAtomContainer molecule) throws SMARTSException, CDKException, CloneNotSupportedException{
		boolean ppc = false;
		
		if(!(isGlycosylatedCompound(molecule) || !isSulfatedCompound(molecule) || !isGlycinatedCompound(molecule))){
			ppc = true;
		}
		
		return ppc;
		
	}
	
	/**
	 * skip unneccessary compounds
	 * Might want to update based on http://eawag-bbd.ethz.ch/servlets/pageservlet?ptype=smallcompsview
	 * @param molecule
	 * @return
	 * @throws SMARTSException
	 */
	public static boolean isUnneccessaryMetabolite(IAtomContainer molecule) throws SMARTSException{
		
		String patterns = "["
				+ "$([H][H]),"	// Dihydrogen
				+ "$(O=C=O)," 	// Carbon dioxide
				+ "$([#6;A;H2X3]=[O;X1]),"   	// Carbon monoxide
				+ "$([OX2H2])," 	// Water
				+ "$([H])," // Hydrogen
				+ "$([OX1]=[SX1]),"
				+ "$([NX3H3]),"
				+ "$([CX4H3])" // Methane
				+ "$([F,Cl,Br,I;-]),"
				+ "$([F,Cl,Br,I;H1X1]),"
				+ "$([#8;X2H1,X1-][S;X4]([#8;X2H1,X1-])(=[O;X1])=[O;X1])," // Sulfate
				+ "$([#8;X2H1,X1-][P;X4]([#8;X2H1,X1-])([#8;X2H1,X1-])=[O;X1])" // Phosphate
				+ "]";
		
		SmartsPatternCDK cdkPattern =  new SmartsPatternCDK(patterns);
		return (cdkPattern.hasSMARTSPattern(molecule)>0);
	}
	
	
	/**
	 * https://www.intechopen.com/books/terpenes-and-terpenoids/introductory-chapter-terpenes-and-terpenoids
	 * isoprene
	 * Monoterpenes
	 * sequitrepene
	 * diterpene
	 * sesterpene
	 * triterpene
	 * tetrateprene
	 * @param mole
	 * @return
	 * @throws SMARTSException 
	 */
	public static boolean isTerpenoid(IAtomContainer mole) throws SMARTSException {
		boolean is_prene = false;
		String isoprene = "[#6]-[#6](=[#6])-[#6]=[#6]";
		SmartsPatternCDK prenePattern 	=  new SmartsPatternCDK(isoprene);	
		if(prenePattern.hasSMARTSPattern(mole) > 0) {
			is_prene = true; 
		}
		return is_prene; 
		
	}
	
	
	/**
	 * Given a molecule A and a SMARTS expression S, list all occurrences of S
	 * in A. Each occurrence is represented as an ArrayList of atom indexes.
	 * 
	 * @param smartsString
	 * @param molecule
	 * @return A list of occurrences of S in A represented as an ArrayList of
	 *         atom indexes
	 * @throws SMARTSException
	 */
	public static List<List<Integer>> findAllOccurences(String smartsString,
		IAtomContainer molecule) throws SMARTSException {
		// might want to use IsomorphismTester as in SMIRKSManager
		List<List<Integer>> matches = new ArrayList<List<Integer>>();
		SmartsPatternCDK smartsPattern = new SmartsPatternCDK(smartsString);
		if (smartsPattern.hasSMARTSPattern(molecule) > 0) {
			matches = smartsPattern.getUniqueMatchingAtoms(molecule);
		}

		return matches;
	}

	
	public static boolean containsCarbon(IAtomContainer molecule) {
		boolean carbon = false;
		for(IAtom at : molecule.atoms()){
			if(at.getAtomicNumber() == 6){
				carbon = true;
				break;
			}
		}	
		return carbon;
	}
	
	
	public static boolean isPolyphenolOrDerivative(IAtomContainer molecule) throws SMARTSException, CDKException, CloneNotSupportedException{
		boolean polyphenol =  false;
		
		if(isMetabolizablePolyphenolOrDerivative(molecule)){
			// This must be removed and the definition of isPolyphenolOrDerivative should be improved. I noticed that
			// Urolithin B (OC1=CC2=C(OC(=O)C3=CC=CC=C23)C=C1) is not a polyphenolOr derivative but isMetabolizablePolyphenolOrDerivative.
			// Which is not very logical !!!
			polyphenol = true;
		} else {
			String constraintsValid = "["
					+ "$(O=[#6;R0](-[#6;R0]-,=[#6;R0]-[#6;R1]-1=[#6;R1]-[#6;R1]=[#6;R1]-[#6;R1]=[#6;R1]-1)-[#6;R1]-1=[#6;R1]-[#6;R1]=[#6;R1]-[#6;R1]=[#6;R1]-1)," // (dihydro)chalcone
					+ "$(O=[#6](-[#6;R0]-,=[#6;R0]-[c;R1]-,:1[c;R1]-,:[c;R1][c;R1]-,:[c;R1][c;R1]-,:1)-[c;R1]-,:1[c;R1]-,:[c;R1][c;R1]-,:[c;R1][c;R1]-,:1)," // (dihydro)chalcone
					+ "$([H][#8;R0]-[#6;R0](-[#6;R0]-[#6;R1]-1=[#6;R1]-[#6;R1]=[#6;R1]-[#6;R1]=[#6;R1]-1)-[#6;R0]-[#6;R1]-1=[#6;R1]-[#6;R1]=[#6;R1]-[#6;R1]=[#6;R1]-1-[#8][H])," // 
					+ "$([H][#8]-[#6;R1]-1=[#6;R1]-[#6;R1]=[#6;R1]-[#6;R1]=[#6;R1]-1-[#6;R0]-[#6;R0](=[O;R0])-[#6;R0]-[#6;R1]-1=[#6;R1]-[#6;R1]=[#6;R1]-[#6;R1]=[#6;R1]-1)," // 
					+ "$([#6;R0](=[#6;R0]/[#6;R1]-1=[#6;R1]-[#6;R1]=[#6;R1]-[#6;R1]=[#6;R1]1)[#6;R1]-1=[#6;R1]-[#6;R1]=[#6;R1]-[#6;R1]=[#6;R1]1)," // stilbene
					+ "$([#1,$([#8;R0]-[#1,C,$(S(=O)(=O)[OX2H1,OX1-])])]-[#6;R1]=,:1[#6;R1](-[$(C(=O)[O;R0]-[#1,C,$(S(=O)(=O)[OX2H1,OX1-])]),$(C([H])=C([H])C(=O)[O;R0]-[#1,C,$(S(=O)(=O)[OX2H1,OX1-])]),$(C([H])([H])C(=O)[O;R0]-[#1,C,$(S(=O)(=O)[OX2H1,OX1-])]),$(C([H])([H])C([H])([H])C(=O)[O;R0]-[#1,C,$(S(=O)(=O)[OX2H1,OX1-])])])=,:[#6;R1](-[#1,$([#8;R0]-[#1,C,$(S(=O)(=O)[OX2H1,OX1-])])])[#6;R1](-[#1,$([#8;R0]-[#1,C,$(S(=O)(=O)[OX2H1,OX1-])])])=,:[#6;R1](-[#1,$([#8;R0]-[#1,C,$(S(=O)(=O)[OX2H1,OX1-])])])[#6;R1]=,:1-[#1,$([#8;R0]-[#1,C,$(S(=O)(=O)[OX2H1,OX1-])])])," // phenyl- and benzoic acids and conjugates
					+ "$([H][#8;X2]-[#6;R1]1=,:[#6;R1](-[R0;*,#1])[#6]2=,:[#6]([#8+]=,:[#6;R1]1-[#6;R1]=,:1[#6;R1](-[R0;*,#1])=,:[#6;R1](-[#1,#8])[#6;R1](-[#1,#8])=,:[#6;R1](-[#1,#8])[#6;R1]=,:1-[R0;*,#1])[#6](-[R0;*,#1])=,:[#6](-[#8;X2])[#6](-[R0;*,#1])=,:[#6]2-[#8;X2])," // anthocyanidin
					+ "$([H][#6;R1]1=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]([H])=,:[#6;R1]1[C;R1]-,:1([H])-,:[#8][#6]=,:2[#6]=,:[#6](!@-[#8;X2])[#6]=,:[#6](!@-[#8;X2])[#6]=,:2-,:[C;R1]([H])([H])-,:[C;R1]-,:1([H])[#8])," // flavan-3-ol
					//+ "$([H][#6;R1]-1=[#6;R1]-[#6;R1]=[#6;R1]-[#6;R1]([H])=[#6;R1]-1[C;R1]1([H])[#8]-[#6]-2=[#6](-[#6](!@-[#8;X2])=[#6]-[#6](!@-[#8;X2])=[#6]-2)[C;R1]([H])([H])[C;R1]1([H])[#8])," // flavan-3-ol
					+ "$([H][#6;R1]1=,:[#6;R1](-[R0;*,#1])[#6;R1](-[R0;*,#1])=,:[#6;R1](-[R0;*,#1])[#6;R1]([H])=,:[#6;R1]1[C;R1]1([H])[#8;R1]-[#6;R2]=,:2[#6;R1](-[R0;*,#1])=,:[#6;R1](-[R0;*,#1])[#6;R1](-[R0;*,#1])=,:[#6;R1](-[R0;*,#1])[#6;R2]=,:2-[#6;R1](=[O;X1])[C;R1]1([H])[#1,OX2H1,OX1-,$([#8]-[#6])])," // flavanone/flavanonol
					+ "$([H][#6;R1]1=,:[#6;R1]([#8][#6]=,:2[#6]=,:[#6][#6]=,:[#6][#6]=,:2[#6;R1]1=[O;X1])-[c;R1]1[c;R1][c;R1][c;R1][c;R1][c;R1]1)," //flavone
					+ "$([H][#8]-[#6;R1]1=,:[#6;R1]([#8][#6]=,:2[#6]=,:[#6][#6]=,:[#6][#6]=,:2[#6;R1]1=[O;X1])-[c;R1]1[c;R1][c;R1][c;R1][c;R1][c;R1]1)," // flavonol
					+ "$([O;X1]=[#6;R1]1[#6;R2]2=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R2]2[#8;R1][#6;R1]=,:[#6;R1]1-[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:1)," // Isoflavone
					+ "$([H][C;R1]1([H])[#8;R1]-[#6;R2]2=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R2]2-[#6;R1](=[O;X1])[C;R1]1([H])[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:1)," // Isoflavanone
					+ "$([O;R0]=[#6;R1]-1-[#6;R1]-[#6;R1]-[#6;R1](-[#6;R0]-[#6;R1]=,:2[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:2)-[#8;R1]-1)," // phenylvalerolactone
					+ "$([H][#6;R1]1=,:[#6;R1]([H])[#6;R1](-[#8;R0]-[#1,#6,#16])=,:[#6;R1](-[#8;R0]-[#1,#6,#16])[#6;R1](-[#8;R0]-[#1,#6,#16])=,:[#6;R1]1[H])," // pyrrogallol or conjugates
					+ "$([H][#6;R1]=,:1[#6;R1]([H])=,:[#6;R1]([H])[#6;R1](-[#8;R0]-[#1,#6,#16])=,:[#6;R1](-[#8;R0]-[#1,#6,#16])[#6;R1]=,:1[H])," // catechol and conjugates
					+ "$([H][#8;R0]-[#6;R1]1=,:[#6;R1]([H])[#6;R1]([H])=,:[#6;R1]([H])[#6;R1]([H])=,:[#6;R1]1-[#8;R0][H])," // catechol
					+ "$([H][#8;R0]-[#6;R1]1=,:[#6;R1]([H])[#6;R1]([H])=,:[#6;R1]([H])[#6;R1](-[#8;R0][H])=,:[#6;R1]1-[#8;R0][H])," // pyrrogallol
					+ "$([H][#8;R0]-[#6;R1]=,:1[#6;R1]([H])=,:[#6;R1](-[#8;R0][H])[#6;R1]([H])=,:[#6;R1](-[#8;R0][H])[#6;R1]=,:1[H])," // phloroglucinol
					+ "$([#8]-[#6;X3](=[O;R0])-[#6]=,:1[#6]=,:[#6][#6]=,:[#6][#6]=,:1)," // benzoic acid
					+ "$([H][#6](=O)-[#6;R1]=,:1[#6;R1](-[#1,$([#8;R0]-[#1,C,$(S(=O)(=O)[OX2H1,OX1-])])])=,:[#6;R1](-[#1,$([#8;R0]-[#1,C,$(S(=O)(=O)[OX2H1,OX1-])])])[#6;R1](-[#1,$([#8;R0]-[#1,C,$(S(=O)(=O)[OX2H1,OX1-])])])=,:[#6;R1](-[#1,$([#8;R0]-[#1,C,$(S(=O)(=O)[OX2H1,OX1-])])])[#6;R1]=,:1-[#1,$([#8;R0]-[#1,C,$(S(=O)(=O)[OX2H1,OX1-])])])," // benzaldehyde derivative and conjugates
					+ "$([#8;R0]-[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1](=,:[#6;R1][#6;R1]=,:1)-[#6;R0]-,=[#6;R0]-[#6;R0](=[O;R0])-[#6;R0]-[#6;R0](=[O;R0])-[#6;R0]-,=[#6;R0]-[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:1)," // curcominoid
					+ "$([H][#8]-[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R2]-2[#6;R2]=,:1-[#6;R2]1=,:[#6;R1]([#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R2]1-[#6](=O)-[#8]C1([H])C([H])([#6;R1]-[#8]-[#6]-2=O)[#8]C([H])([#8][H])C2([H])[#8]-[#6](=O)-[#6;R2]3=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1](-[#8][H])=,:[#6;R2]3-[#6;R2]3=,:[#6;R1]([#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R2]3-[#6](=O)-[#8]C12[H])-[#8][H])-[#8][H])," // ellagitannin pattern1
					+ "$([H][#8]-[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R2]-2[#6;R2]=,:1-[#6;R2]1=,:[#6;R1]([#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R2]1-[#6](=O)-[#8]C1([H])C([H])([#8])C([H])([#8]C([H])([#8][H])C1([H])[#8]-[#6]-2=O)C([H])([H])[#8])-[#8][H])," // ellagitannin pattern2
					+ "$([H][#8]-[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R2]-2[#6;R2]=,:1-[#6;R2]1=,:[#6;R1]([#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R2]1-[#6](=O)-[#8]C1([H])C([H])([#6;R1]-[#8]-[#6]-2=O)[#8]C([H])([#8][H])C([H])([#8])C1([H])[#8])-[#8][H])," // ellagitannin pattern3
					+ "$([H][#8]!@-[#6]1=,:[#6][#6]2=,:[#6]3[#6]([#8][#6](=O)[#6]4=,:[#6][#6](!@-[#8][H])=,:[#6](!@-[#8][H])[#6]([#8][#6]2=O)=,:[#6]34)=,:[#6]1!@-[#8][H])," // ellagic acid
					+ "$([H][#8]-[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R2]2[#6;R2]=,:3[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R2]=,:3[#8;R1][#6;R1](=[O;R0])[#6;R2]=,:12)," // Urolithin backbone
					+ "$([H][#8]-[#6;R1]1=,:[#6;R1][#6;R1]=,:[#6;R2]2[#6;R2]=,:3[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R2]=,:3[#8;R1][#6;R1](=[O;R0])[#6;R2]2=,:[#6;R1]1)," // Urolithin backbone
					+ "$([H][#8]-[#6;R1]1=,:[#6;R1][#6;R1]=,:[#6;R2]2[#6;R2](=,:[#6;R1]1)[#6;R2]=,:1[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R2]=,:1[#8;R1][#6;R1]2=[O;R0])," // Urolithin backbone
					+ "$([H][#8]-[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R2]2[#6;R2]=,:1[#6;R2]=,:1[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R2]=,:1[#8;R1][#6;R1]2=[O;R0])," // Urolithin backbone
					+ "$([H][#8]-[#6;R1]1=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R2]=,:2[#8;R1][#6;R1](=[O;R0])[#6;R2]3=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R2]3[#6;R2]1=,:2)," // Urolithin backbone
					+ "$([H][#8]-[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R2]=,:2[#8;R1][#6;R1](=[O;R0])[#6;R2]3=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R2]3[#6;R2]=,:2[#6;R1]=,:1)," // Urolithin backbone
					+ "$([H][#8]-[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R2]2=,:[#6;R2]([#6;R1]=,:1)[#8;R1][#6;R1](=[O;R0])[#6;R2]1=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R2]21)," // Urolithin backbone
					+ "$([H][#8]-[#6;R1]1=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R2]=,:2[#6;R2]3=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R2]3[#6;R1](=[O;R0])[#8;R1][#6;R2]1=,:2)," // Urolithin backbone
					+ "$([H][#8]-[#6;R0](=[O;R0])-[#6;R1]=,:1[#6;R2]=,:[#6;R2][#6;R1]=,:[#6;R1][#6;R1]=,:1)," // fused carboxylated benzene
					+ "$([H][#8]-[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R2][#6;R2]=,:1)," // fused hydroxylated benzene
					+ "$([H][#8]-[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R2]=,:[#6;R2][#6;R1]=,:1),"  // fused hydroxylated benzene
					+ "$([#8;A;X2H1,X1-][#6](=O)-[#6;R0]-[#6;R0]-[#6;R0]-[#6;R0]-[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:1)," // phenylvaleric acid
					+ "$([#8;A;X2H1,X1-][#6](=O)-[#6;R0]-[#6;R0]-[#6;R0]-[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:1)," // phenylbutyric acid
									
					// benzo[c]chromen‐6‐one (Urolithins must be hydroxylated at one or more positions)
					// Episin, J.C. (2013); Biological Significance of Urolithins, the Gut Microbial Ellagic Acid-Derived Metabolites: The Evidence So Far; 
					// Evidence-Based Complementary and Alternative Medicine; Volume 2013, Article ID 270418; 15 pages; http://dx.doi.org/10.1155/2013/270418
					+ "$([H][#8]-[#6;R1]1=,:[#6;R1]([H])[#6;R2]=,:2[#8;R1][#6;R1](=[O;R0])[#6;R2]=,:3[#6;R2]([#6;R2]=,:2[#6;R1]([H])=,:[#6;R1]1[H])=,:[#6;R1]([H])[#6;R1]([H])=,:[#6;R1]([H])[#6;R1]=,:3[H])"
	
					+ "]";
					// add depside
			
			SmartsPatternCDK smartsPatternValid = new SmartsPatternCDK(constraintsValid);	
			String constraintsInvalid_3 ="[$([H][#8][C;R1]1([H])[C;R1]([H])([H])[#6]=,:2[#6]=,:[#6][#6]=,:[#6](-[CX3,c])[#6]=,:2-[#8][C;R1]1([H])[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:1),$([H][#8][C;R1]1([H])[C;R1]([H])([#6;CX3,c])[#6]=,:2[#6]=,:[#6][#6]=,:[#6][#6]=,:2-[#8][C;R1]1([H])[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:1)]";
					
			String athocyanidin = "[$([#8]-[#6;R1]-1=[#6;R1](-[#8;R1]-[#6]-2=[#6]-[#6](=O)-[#6]=[#6]-[#6]-2=[#6;R1]-1)-[c;R1]1[c;R1][c;R1]c([#8;A;H1X2])[c;R1][c;R1]1),"
					+ "$([#8]-[#6;R1]1=,:[#6;R1]c2ccc([#8;A;H1X2])cc2[#8;R1+]=,:[#6;R1]1-[c;R1]1[c;R1][c;R1]c([#8;A;H1X2])[c;R1][c;R1]1),"
					+ "$([#8]-[#6;R1]1=,:[#6;R1]c2ccc([#8;A;H1X2])cc2[#8;R1+]=,:[#6;R1]1-[c;R1]1[c;R1][c;R1]c([#8;A;H1X2])[c;R1][c;R1]1),"
					+ "$([#8;A;H1X2][#6]-1=[#6;R1]-[#6;R1]=[#6;R1](-[#6;R1]=[#6;R1]-1)[C;R1]1([#8;A;H1X2])[#8;R1]-[#6]-2=[#6]-[#6]([#8;A;H1X2])=[#6]-[#6]=[#6]-2-[#6;R1]=[#6;R1]1-[#8]),"
					+ "$([H][#8;X2]-[#6;R1]1=,:[#6;R1](-[R0;*,#1])[#6]2=,:[#6]([#8+]=,:[#6;R1]1-[#6;R1]=,:1[#6;R1](-[R0;*,#1])=,:[#6;R1](-[R0;#1,$([O][H]),$([O]-[C]([H])[H])])[#6;R1](-"
					+ "[R0;#1,$([O][H]),$([O]-[C]([H])[H])])=,:[#6;R1](-[R0;#1,$([O][H]),$([O]-[C]([H])[H])])[#6;R1]=,:1-[R0;*,#1])[#6](-[R0;*,#1])=,:[#6](-[#8;X2]-[R0;#1,$([C](-[H])(-[H])[H])])[#6](-[R0;*,#1])=,:[#6]2-[#8;X2]-[R0;#1,$([C](-[H])(-[H])[H])])]";
	
			SmartsPatternCDK flavonoidPattern 		=  new SmartsPatternCDK("[$([#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1](=,:[#6;R1][#6;R1]=,:1)!@-[#6]-1-,=[#6;R1]-[#6;R1]-[#6]=,:2[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6]=,:2-[#8;R1]-1),"
					+ "$(O=[#6]1[#6]=,:[#6]([#8][#6]2=,:[#6][#6]=,:[#6][#6]=,:[#6]12)-[#6]=,:1[#6]=,:[#6][#6]=,:[#6][#6]=,:1)"
					+ "]");	                      
			
	
			SmartsPatternCDK otherFlavonoidsPattern = new SmartsPatternCDK("["
			                         + "$([H][#6]=,:1[#6]=,:[#6][#6]([H])=,:[#6]([#6]=,:1[H])-[#6]1=,:[#6][#6][#6]2=,:[#6]([#8]1)[#6]([H])=,:[#6](-[#8])[#6]([H])=,:[#6]2-[#8]),"
			                         + "$([H][#6]=,:1[#6]=,:[#6][#6]([H])=,:[#6]([#6]=,:1[H])C1([H])[#8]-[#6]=,:2[#6]([H])=,:[#6](-[#8])[#6]([H])=,:[#6](-[#8])[#6]=,:2-[#6]-[#6]1[H]),"
			                         + "$([H][#6]=,:1[#6]=,:[#6][#6]([H])=,:[#6]([#6]=,:1[H])-[#6]=,:1[#6][#6]2=,:[#6]([#8][#6]=,:1[H])[#6]([H])=,:[#6](-[#8])[#6]([H])=,:[#6]2-[#8]),"
			                         + "$([H][#6]=,:1[#6]=,:[#6][#6]([H])=,:[#6]([#6]=,:1[H])C1([H])[#6]-[#6]=,:2[#6](-[#8])=,:[#6]([H])[#6](-[#8])=,:[#6]([H])[#6]=,:2-[#8]C1([H])[H])"
			                         + "]");
			SmartsPatternCDK anthocyanidinPattern = new SmartsPatternCDK(athocyanidin);
			SmartsPatternCDK isoflavonoidPattern 	=  new SmartsPatternCDK("[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1](=,:[#6;R1][#6;R1]=,:1)!@-[#6;R1]-1-,=[#6]-[#8;R1]-[#6]=,:2[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6]=,:2-[#6;R1]-1");
			SmartsPatternCDK isoflavonePattern = new SmartsPatternCDK("[R0;*,#1]-[#6;R1]=,:1[#8]c2c([#6;R1](=[O;X1])[#6;R1]=,:1-[c;R1]1[c;R1](-[R0;*,#1])[c;R1](-[R0;*,#1])[c;R1](-[R0;*,#1])[c;R1](-[R0;*,#1])[c;R1]1-[R0;*,#1])c(-[R0;*,#1])c(-[R0;*,#1])c(-[R0;*,#1])c2-[R0;*,#1]");
			
			if( smartsPatternValid.hasSMARTSPattern(molecule)>0 || anthocyanidinPattern.hasSMARTSPattern(molecule)>0 || 
				 flavonoidPattern.hasSMARTSPattern(molecule)>0 || otherFlavonoidsPattern.hasSMARTSPattern(molecule)>0 ||
				 isoflavonoidPattern.hasSMARTSPattern(molecule)>0 || isoflavonePattern.hasSMARTSPattern(molecule)>0 || 			
					ChemStructureExplorer.findAllOccurences(constraintsInvalid_3,molecule).size()>0) {
					polyphenol = true;
			}
			
			
		}
			
		return polyphenol;
	}
	
	
	public static boolean isMetabolizablePolyphenolOrDerivative(IAtomContainer mol) throws SMARTSException, CDKException, CloneNotSupportedException{
		IAtomContainer molecule = mol.clone();
		AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(molecule);
		AtomContainerManipulator.convertImplicitToExplicitHydrogens(molecule);
		
		boolean polyphenol =  false;
		String constraintsValid = "["
				+ "$(O=[#6;R0](-[#6;R0]-,=[#6;R0]-[#6;R1]-1=[#6;R1]-[#6;R1]=[#6;R1]-[#6;R1]=[#6;R1]-1)-[#6;R1]-1=[#6;R1]-[#6;R1]=[#6;R1]-[#6;R1]=[#6;R1]-1)," // (dihydro)chalcone
				+ "$(O=[#6](-[#6;R0]-,=[#6;R0]-[c;R1]-,:1[c;R1]-,:[c;R1][c;R1]-,:[c;R1][c;R1]-,:1)-[c;R1]-,:1[c;R1]-,:[c;R1][c;R1]-,:[c;R1][c;R1]-,:1)," // (dihydro)chalcone
				+ "$([H][#8;R0]-[#6;R0](-[#6;R0]-[#6;R1]-1=[#6;R1]-[#6;R1]=[#6;R1]-[#6;R1]=[#6;R1]-1)-[#6;R0]-[#6;R1]-1=[#6;R1]-[#6;R1]=[#6;R1]-[#6;R1]=[#6;R1]-1-[#8][H])," // 
				+ "$([H][#8]-[#6;R1]-1=[#6;R1]-[#6;R1]=[#6;R1]-[#6;R1]=[#6;R1]-1-[#6;R0]-[#6;R0](=[O;R0])-[#6;R0]-[#6;R1]-1=[#6;R1]-[#6;R1]=[#6;R1]-[#6;R1]=[#6;R1]-1)," // 
				+ "$([#6;R0](=[#6;R0]/[#6;R1]-1=[#6;R1]-[#6;R1]=[#6;R1]-[#6;R1]=[#6;R1]1)[#6;R1]-1=[#6;R1]-[#6;R1]=[#6;R1]-[#6;R1]=[#6;R1]1)," // stilbene
				+ "$([#1,$([#8;R0]-[#1,C,$(S(=O)(=O)[OX2H1,OX1-])])]-[#6;R1]=,:1[#6;R1](-[$(C(=O)[O;R0]-[#1,C,$(S(=O)(=O)[OX2H1,OX1-])]),$(C([H])=C([H])C(=O)[O;R0]-[#1,C,$(S(=O)(=O)[OX2H1,OX1-])]),$(C([H])([H])C(=O)[O;R0]-[#1,C,$(S(=O)(=O)[OX2H1,OX1-])]),$(C([H])([H])C([H])([H])C(=O)[O;R0]-[#1,C,$(S(=O)(=O)[OX2H1,OX1-])])])=,:[#6;R1](-[#1,$([#8;R0]-[#1,C,$(S(=O)(=O)[OX2H1,OX1-])])])[#6;R1](-[#1,$([#8;R0]-[#1,C,$(S(=O)(=O)[OX2H1,OX1-])])])=,:[#6;R1](-[#1,$([#8;R0]-[#1,C,$(S(=O)(=O)[OX2H1,OX1-])])])[#6;R1]=,:1-[#1,$([#8;R0]-[#1,C,$(S(=O)(=O)[OX2H1,OX1-])])])," // phenyl- and benzoic acids and conjugates
				+ "$([H][#8;X2]-[#6;R1]1=,:[#6;R1](-[R0;*,#1])[#6]2=,:[#6]([#8+]=,:[#6;R1]1-[#6;R1]=,:1[#6;R1](-[R0;*,#1])=,:[#6;R1](-[#1,#8])[#6;R1](-[#1,#8])=,:[#6;R1](-[#1,#8])[#6;R1]=,:1-[R0;*,#1])[#6](-[R0;*,#1])=,:[#6](-[#8;X2])[#6](-[R0;*,#1])=,:[#6]2-[#8;X2])," // anthocyanidin
				+ "$([H][#6;R1]1=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]([H])=,:[#6;R1]1[C;R1]-,:1([H])-,:[#8][#6]=,:2[#6]=,:[#6](!@-[#8;X2])[#6]=,:[#6](!@-[#8;X2])[#6]=,:2-,:[C;R1]([H])([H])-,:[C;R1]-,:1([H])[#8])," // flavan-3-ol
				//+ "$([H][#6;R1]-1=[#6;R1]-[#6;R1]=[#6;R1]-[#6;R1]([H])=[#6;R1]-1[C;R1]1([H])[#8]-[#6]-2=[#6](-[#6](!@-[#8;X2])=[#6]-[#6](!@-[#8;X2])=[#6]-2)[C;R1]([H])([H])[C;R1]1([H])[#8])," // flavan-3-ol
				+ "$([H][#6;R1]1=,:[#6;R1](-[R0;*,#1])[#6;R1](-[R0;*,#1])=,:[#6;R1](-[R0;*,#1])[#6;R1]([H])=,:[#6;R1]1[C;R1]1([H])[#8;R1]-[#6;R2]=,:2[#6;R1](-[R0;*,#1])=,:[#6;R1](-[R0;*,#1])[#6;R1](-[R0;*,#1])=,:[#6;R1](-[R0;*,#1])[#6;R2]=,:2-[#6;R1](=[O;X1])[C;R1]1([H])[#1,OX2H1,OX1-,$([#8]-[#6])])," // flavanone/flavanonol
				+ "$([H][#6;R1]1=,:[#6;R1]([#8][#6]=,:2[#6]=,:[#6][#6]=,:[#6][#6]=,:2[#6;R1]1=[O;X1])-[c;R1]1[c;R1][c;R1][c;R1][c;R1][c;R1]1)," //flavone
				+ "$([H][#8]-[#6;R1]1=,:[#6;R1]([#8][#6]=,:2[#6]=,:[#6][#6]=,:[#6][#6]=,:2[#6;R1]1=[O;X1])-[c;R1]1[c;R1][c;R1][c;R1][c;R1][c;R1]1)," // flavonol
				+ "$([O;X1]=[#6;R1]1[#6;R2]2=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R2]2[#8;R1][#6;R1]=,:[#6;R1]1-[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:1)," // Isoflavone
				+ "$([H][C;R1]1([H])[#8;R1]-[#6;R2]2=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R2]2-[#6;R1](=[O;X1])[C;R1]1([H])[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:1)," // Isoflavanone
				+ "$([O;R0]=[#6;R1]-1-[#6;R1]-[#6;R1]-[#6;R1](-[#6;R0]-[#6;R1]=,:2[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:2)-[#8;R1]-1)," // phenylvalerolactone
				+ "$([H][#6;R1]1=,:[#6;R1]([H])[#6;R1](-[#8;R0]-[#1,#6,#16])=,:[#6;R1](-[#8;R0]-[#1,#6,#16])[#6;R1](-[#8;R0]-[#1,#6,#16])=,:[#6;R1]1[H])," // pyrrogallol or conjugates
				+ "$([H][#6;R1]=,:1[#6;R1]([H])=,:[#6;R1]([H])[#6;R1](-[#8;R0]-[#1,#6,#16])=,:[#6;R1](-[#8;R0]-[#1,#6,#16])[#6;R1]=,:1[H])," // catechol and conjugates
				+ "$([H][#8;R0]-[#6;R1]1=,:[#6;R1]([H])[#6;R1]([H])=,:[#6;R1]([H])[#6;R1]([H])=,:[#6;R1]1-[#8;R0][H])," // catechol
				+ "$([H][#8;R0]-[#6;R1]1=,:[#6;R1]([H])[#6;R1]([H])=,:[#6;R1]([H])[#6;R1](-[#8;R0][H])=,:[#6;R1]1-[#8;R0][H])," // pyrrogallol
				+ "$([H][#8;R0]-[#6;R1]=,:1[#6;R1]([H])=,:[#6;R1](-[#8;R0][H])[#6;R1]([H])=,:[#6;R1](-[#8;R0][H])[#6;R1]=,:1[H])," // phloroglucinol
				+ "$([#8]-[#6;X3](=[O;R0])-[#6]=,:1[#6]=,:[#6][#6]=,:[#6][#6]=,:1)," // benzoic acid
				+ "$([H][#6](=O)-[#6;R1]=,:1[#6;R1](-[#1,$([#8;R0]-[#1,C,$(S(=O)(=O)[OX2H1,OX1-])])])=,:[#6;R1](-[#1,$([#8;R0]-[#1,C,$(S(=O)(=O)[OX2H1,OX1-])])])[#6;R1](-[#1,$([#8;R0]-[#1,C,$(S(=O)(=O)[OX2H1,OX1-])])])=,:[#6;R1](-[#1,$([#8;R0]-[#1,C,$(S(=O)(=O)[OX2H1,OX1-])])])[#6;R1]=,:1-[#1,$([#8;R0]-[#1,C,$(S(=O)(=O)[OX2H1,OX1-])])])," // benzaldehyde derivative and conjugates
				+ "$([#8;R0]-[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1](=,:[#6;R1][#6;R1]=,:1)-[#6;R0]-,=[#6;R0]-[#6;R0](=[O;R0])-[#6;R0]-[#6;R0](=[O;R0])-[#6;R0]-,=[#6;R0]-[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:1)," // curcominoid
				+ "$([H][#8]-[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R2]-2[#6;R2]=,:1-[#6;R2]1=,:[#6;R1]([#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R2]1-[#6](=O)-[#8]C1([H])C([H])([#6;R1]-[#8]-[#6]-2=O)[#8]C([H])([#8][H])C2([H])[#8]-[#6](=O)-[#6;R2]3=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1](-[#8][H])=,:[#6;R2]3-[#6;R2]3=,:[#6;R1]([#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R2]3-[#6](=O)-[#8]C12[H])-[#8][H])-[#8][H])," // ellagitannin pattern1
				+ "$([H][#8]-[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R2]-2[#6;R2]=,:1-[#6;R2]1=,:[#6;R1]([#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R2]1-[#6](=O)-[#8]C1([H])C([H])([#8])C([H])([#8]C([H])([#8][H])C1([H])[#8]-[#6]-2=O)C([H])([H])[#8])-[#8][H])," // ellagitannin pattern2
				+ "$([H][#8]-[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R2]-2[#6;R2]=,:1-[#6;R2]1=,:[#6;R1]([#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R2]1-[#6](=O)-[#8]C1([H])C([H])([#6;R1]-[#8]-[#6]-2=O)[#8]C([H])([#8][H])C([H])([#8])C1([H])[#8])-[#8][H])," // ellagitannin pattern3
				+ "$([H][#8]!@-[#6]1=,:[#6][#6]2=,:[#6]3[#6]([#8][#6](=O)[#6]4=,:[#6][#6](!@-[#8][H])=,:[#6](!@-[#8][H])[#6]([#8][#6]2=O)=,:[#6]34)=,:[#6]1!@-[#8][H])," // ellagic acid
				+ "$([H][#8]-[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R2]2[#6;R2]=,:3[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R2]=,:3[#8;R1][#6;R1](=[O;R0])[#6;R2]=,:12)," // Urolithin backbone
				+ "$([H][#8]-[#6;R1]1=,:[#6;R1][#6;R1]=,:[#6;R2]2[#6;R2]=,:3[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R2]=,:3[#8;R1][#6;R1](=[O;R0])[#6;R2]2=,:[#6;R1]1)," // Urolithin backbone
				+ "$([H][#8]-[#6;R1]1=,:[#6;R1][#6;R1]=,:[#6;R2]2[#6;R2](=,:[#6;R1]1)[#6;R2]=,:1[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R2]=,:1[#8;R1][#6;R1]2=[O;R0])," // Urolithin backbone
				+ "$([H][#8]-[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R2]2[#6;R2]=,:1[#6;R2]=,:1[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R2]=,:1[#8;R1][#6;R1]2=[O;R0])," // Urolithin backbone
				+ "$([H][#8]-[#6;R1]1=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R2]=,:2[#8;R1][#6;R1](=[O;R0])[#6;R2]3=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R2]3[#6;R2]1=,:2)," // Urolithin backbone
				+ "$([H][#8]-[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R2]=,:2[#8;R1][#6;R1](=[O;R0])[#6;R2]3=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R2]3[#6;R2]=,:2[#6;R1]=,:1)," // Urolithin backbone
				+ "$([H][#8]-[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R2]2=,:[#6;R2]([#6;R1]=,:1)[#8;R1][#6;R1](=[O;R0])[#6;R2]1=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R2]21)," // Urolithin backbone
				+ "$([H][#8]-[#6;R1]1=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R2]=,:2[#6;R2]3=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R2]3[#6;R1](=[O;R0])[#8;R1][#6;R2]1=,:2)," // Urolithin backbone
				+ "$([H][#8]-[#6;R0](=[O;R0])-[#6;R1]=,:1[#6;R2]=,:[#6;R2][#6;R1]=,:[#6;R1][#6;R1]=,:1)," // fused carboxylated benzene
				+ "$([H][#8]-[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R2][#6;R2]=,:1)," // fused hydroxylated benzene
				+ "$([H][#8]-[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R2]=,:[#6;R2][#6;R1]=,:1),"  // fused hydroxylated benzene
				+ "$([#8;A;X2H1,X1-][#6](=O)-[#6;R0]-[#6;R0]-[#6;R0]-[#6;R0]-[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:1)," // phenylvaleric acid
				+ "$([#8;A;X2H1,X1-][#6](=O)-[#6;R0]-[#6;R0]-[#6;R0]-[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:1)," // phenylbutyric acid
				
				// benzo[c]chromen‐6‐one (Urolithins must be hydroxylated at one or more positions)
				// Episin, J.C. (2013); Biological Significance of Urolithins, the Gut Microbial Ellagic Acid-Derived Metabolites: The Evidence So Far; 
				// Evidence-Based Complementary and Alternative Medicine; Volume 2013, Article ID 270418; 15 pages; http://dx.doi.org/10.1155/2013/270418
				+ "$([H][#8]-[#6;R1]1=,:[#6;R1]([H])[#6;R2]=,:2[#8;R1][#6;R1](=[O;R0])[#6;R2]=,:3[#6;R2]([#6;R2]=,:2[#6;R1]([H])=,:[#6;R1]1[H])=,:[#6;R1]([H])[#6;R1]([H])=,:[#6;R1]([H])[#6;R1]=,:3[H])"
				+ "]";
				// add depside
		
		SmartsPatternCDK smartsPatternValid = new SmartsPatternCDK(constraintsValid);
		

		
		
		
		/**
		 * (R1) Marín, L. et al. (2015); Bioavailability of Dietary Polyphenols and Gut Microbiota Metabolism: Antimicrobial Properties; Biomed Res Int. 2015; 2015: 905215.; doi:  10.1155/2015/905215
		 * (R2) Deprez, S. et al. (2001); “Transport of proanthocyanidin dimer, trimer, and polymer across monolayers of human intestinal epithelial Caco-2 cells,” Antioxidants and Redox Signaling, vol. 3, no. 6, pp. 957–967
		 * (R3) Monagas, M. et al. (2010); “Insights into the metabolism and microbial biotransformation of dietary  avan-3-ols and the bioactivity of their metabolites,” Food and Function, vol. 1, no. 3, pp. 233–253.
		 * 
		 * Oligomers with a degree of polymerization >3 are not absorbed in the small intestine, and therefore they are metabolized in the colon (R1-R3).
		 */
		String constraintsInvalid_3 ="[$([H][#8][C;R1]1([H])[C;R1]([H])([H])[#6]=,:2[#6]=,:[#6][#6]=,:[#6](-[CX3,c])[#6]=,:2-[#8][C;R1]1([H])[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:1),$([H][#8][C;R1]1([H])[C;R1]([H])([#6;CX3,c])[#6]=,:2[#6]=,:[#6][#6]=,:[#6][#6]=,:2-[#8][C;R1]1([H])[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:1)]";
				
		String athocyanidin = "[$([#8]-[#6;R1]-1=[#6;R1](-[#8;R1]-[#6]-2=[#6]-[#6](=O)-[#6]=[#6]-[#6]-2=[#6;R1]-1)-[c;R1]1[c;R1][c;R1]c([#8;A;H1X2])[c;R1][c;R1]1),"
				+ "$([#8]-[#6;R1]1=,:[#6;R1]c2ccc([#8;A;H1X2])cc2[#8;R1+]=,:[#6;R1]1-[c;R1]1[c;R1][c;R1]c([#8;A;H1X2])[c;R1][c;R1]1),"
				+ "$([#8]-[#6;R1]1=,:[#6;R1]c2ccc([#8;A;H1X2])cc2[#8;R1+]=,:[#6;R1]1-[c;R1]1[c;R1][c;R1]c([#8;A;H1X2])[c;R1][c;R1]1),"
				+ "$([#8;A;H1X2][#6]-1=[#6;R1]-[#6;R1]=[#6;R1](-[#6;R1]=[#6;R1]-1)[C;R1]1([#8;A;H1X2])[#8;R1]-[#6]-2=[#6]-[#6]([#8;A;H1X2])=[#6]-[#6]=[#6]-2-[#6;R1]=[#6;R1]1-[#8]),"
				+ "$([H][#8;X2]-[#6;R1]1=,:[#6;R1](-[R0;*,#1])[#6]2=,:[#6]([#8+]=,:[#6;R1]1-[#6;R1]=,:1[#6;R1](-[R0;*,#1])=,:[#6;R1](-[R0;#1,$([O][H]),$([O]-[C]([H])[H])])[#6;R1](-"
				+ "[R0;#1,$([O][H]),$([O]-[C]([H])[H])])=,:[#6;R1](-[R0;#1,$([O][H]),$([O]-[C]([H])[H])])[#6;R1]=,:1-[R0;*,#1])[#6](-[R0;*,#1])=,:[#6](-[#8;X2]-[R0;#1,$([C](-[H])(-[H])[H])])[#6](-[R0;*,#1])=,:[#6]2-[#8;X2]-[R0;#1,$([C](-[H])(-[H])[H])])]";
		
		SmartsPatternCDK otherFlavonoidsPattern = new SmartsPatternCDK("["
		                         + "$([H][#6]=,:1[#6]=,:[#6][#6]([H])=,:[#6]([#6]=,:1[H])-[#6]1=,:[#6][#6][#6]2=,:[#6]([#8]1)[#6]([H])=,:[#6](-[#8])[#6]([H])=,:[#6]2-[#8]),"
		                         + "$([H][#6]=,:1[#6]=,:[#6][#6]([H])=,:[#6]([#6]=,:1[H])C1([H])[#8]-[#6]=,:2[#6]([H])=,:[#6](-[#8])[#6]([H])=,:[#6](-[#8])[#6]=,:2-[#6]-[#6]1[H]),"
		                         + "$([H][#6]=,:1[#6]=,:[#6][#6]([H])=,:[#6]([#6]=,:1[H])-[#6]=,:1[#6][#6]2=,:[#6]([#8][#6]=,:1[H])[#6]([H])=,:[#6](-[#8])[#6]([H])=,:[#6]2-[#8]),"
		                         + "$([H][#6]=,:1[#6]=,:[#6][#6]([H])=,:[#6]([#6]=,:1[H])C1([H])[#6]-[#6]=,:2[#6](-[#8])=,:[#6]([H])[#6](-[#8])=,:[#6]([H])[#6]=,:2-[#8]C1([H])[H])"
		                         + "]");
		                      
		SmartsPatternCDK anthocyanidinPattern = new SmartsPatternCDK(athocyanidin);
		SmartsPatternCDK isoflavonePattern = new SmartsPatternCDK("[R0;*,#1]-[#6;R1]=,:1[#8]c2c([#6;R1](=[O;X1])[#6;R1]=,:1-[c;R1]1[c;R1](-[R0;*,#1])[c;R1](-[R0;*,#1])[c;R1](-[R0;*,#1])[c;R1](-[R0;*,#1])[c;R1]1-[R0;*,#1])c(-[R0;*,#1])c(-[R0;*,#1])c(-[R0;*,#1])c2-[R0;*,#1]");
		SmartsPatternCDK sulfatedRadicalPattern = new SmartsPatternCDK("[#6]-[#8;X2]S([#8;A;X2H1,X1-])(=[O;X1])=[O;X1]");
		SmartsPatternCDK glucuronidePattern = new SmartsPatternCDK("[#8;R0]-[#6;R1]-1-[#8;R1]-[#6;R1](-[#6;R1](-[#8;R0])-[#6;R1](-[#8;R0])-[#6;R1]-1-[#8;R0])-[#6](-[#8])=O");

		SmartsPatternCDK glycosylMoietyPattern = new SmartsPatternCDK("["
				+ "$(CC1OC(O)C(O)C(O)C1O),"
				+ "$(CC1OC(O)C(O)CC1O),"
				+ "$([#6]!@-[#6]-1-[#8]-[#6](!@-[#8])-[#6](!@-[*,#1;OX2H1,$(NC(=O)C)])-[#6](!@-[#8])-[#6]-1!@-[#8]),"
				+ "$([#6]!@-[#6]-1-[#8]-[#6](!@-[#8])-[#6](!@-[OX2H1,$(NC(=O)C)])-[#6]-[#6]-1!@-[#8])"
				+ "]"
				);
		
		SmartsPatternCDK oMethylPattern 		=  new SmartsPatternCDK("[#6;A;H3X4][#8;X2R0]-[#6;R1]");		
		SmartsPatternCDK flavonoidPattern 	=  new SmartsPatternCDK("[$([#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1](=,:[#6;R1][#6;R1]=,:1)!@-[#6]-1-,=[#6;R1]-[#6;R1]-[#6]=,:2[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6]=,:2-[#8;R1]-1),"
				+ "$(O=[#6]1[#6]=,:[#6]([#8][#6]2=,:[#6][#6]=,:[#6][#6]=,:[#6]12)-[#6]=,:1[#6]=,:[#6][#6]=,:[#6][#6]=,:1)"
				+ "]");		
		
		SmartsPatternCDK isoflavonoidPattern 	=  new SmartsPatternCDK("[#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1](=,:[#6;R1][#6;R1]=,:1)!@-[#6;R1]-1-,=[#6]-[#8;R1]-[#6]=,:2[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6]=,:2-[#6;R1]-1");
		
		glycosylMoietyPattern.match(molecule);
		oMethylPattern.match(molecule);
		sulfatedRadicalPattern.match(molecule);
		
		if( (flavonoidPattern.hasSMARTSPattern(molecule)>0 || isoflavonoidPattern.hasSMARTSPattern(molecule)>0 || anthocyanidinPattern.hasSMARTSPattern(molecule)>0 || otherFlavonoidsPattern.hasSMARTSPattern(molecule)>0)  
				&& (glycosylMoietyPattern.getUniqueMatchingAtoms().size() +
				sulfatedRadicalPattern.getUniqueMatchingAtoms().size() >= 2)){
			
			// I removed oMethylPattern.getUniqueMatchingAtoms().size() from the constraint  because some compounds have 2 or more methyl groups when
			// undergoing reduction.
			
			polyphenol = false;
//			System.err.println("(flavonoidPattern.hasSMARTSPattern(molecule)>0 || isoflavonoidPattern.hasSMARTSPattern(molecule)>0 || anthocyanidinPattern.hasSMARTSPattern(molecule)>0 || otherFlavonoidsPattern.hasSMARTSPattern(molecule)>0)   && (glycosylMoietyPattern.getUniqueMatchingAtoms().size() + oMethylPattern.getUniqueMatchingAtoms().size() + sulfatedRadicalPattern.getUniqueMatchingAtoms().size() >= 2)");
		}
		
		else if ( (ChemStructureExplorer.findAllOccurences(constraintsInvalid_3,molecule).size() <= 3) && (smartsPatternValid.hasSMARTSPattern(molecule) > 0)) {
			polyphenol = true;
		}
		
		return polyphenol;
	}


	public static boolean isGlycosylatedCompound(IAtomContainer molecule) throws SMARTSException {		
		SmartsPatternCDK glucuronidePattern = new SmartsPatternCDK("[#8;R0]-[#6;R1]-1-[#8;R1]-[#6;R1](-[#6;R1](-[#8;R0])-[#6;R1](-[#8;R0])-[#6;R1]-1-[#8;R0])-[#6](-[#8])=O");
		SmartsPatternCDK glycosylMoietyPattern = new SmartsPatternCDK("["
				+ "$(CC1OC(O)C(O)C(O)C1O),"
				+ "$(CC1OC(O)C(O)CC1O),"
				+ "$([#6]!@-[#6]-1-[#8]-[#6](!@-[#8])-[#6](!@-[*,#1;OX2H1,$(NC(=O)C)])-[#6](!@-[#8])-[#6]-1!@-[#8]),"
				+ "$([#6]!@-[#6]-1-[#8]-[#6](!@-[#8])-[#6](!@-[OX2H1,$(NC(=O)C)])-[#6]-[#6]-1!@-[#8])"
				+ "]"
				);
		
		glucuronidePattern.match(molecule);
		glycosylMoietyPattern.match(molecule);		
		boolean isGlycosylated = glucuronidePattern.hasSMARTSPattern(molecule)>0 || glycosylMoietyPattern.hasSMARTSPattern(molecule)>0 ;
			
		return isGlycosylated;
	}
	
	public static boolean isSulfatedCompound(IAtomContainer molecule) throws SMARTSException {
		SmartsPatternCDK sulfatedRadicalPattern = new SmartsPatternCDK("[#6]-[#8;X2]S([#8;A;X2H1,X1-])(=[O;X1])=[O;X1]");	
		return sulfatedRadicalPattern.hasSMARTSPattern(molecule)>0;	
	}
		
	public static boolean isGlutathioneConjugate(IAtomContainer molecule) throws SMARTSException{
		SmartsPatternCDK glutathioneConjugatePattern = new SmartsPatternCDK("[H][#7]([#6;A;H2X4][#6](-[#8])=O)-[#6](=O)[#6;A;H1X4]([#6;A;H2X4][#16][#6,#8,#16;A])[#7]([H])-[#6](=O)[#6;A;H2X4][#6;A;H2X4][#6;A;H1X4]([#7])[#6](-[#8])=O");	
		return glutathioneConjugatePattern.hasSMARTSPattern(molecule)>0;	
	}
	
	public static boolean isGlycinatedCompound(IAtomContainer molecule) throws SMARTSException {
		SmartsPatternCDK glycinatedRadicalPattern = new SmartsPatternCDK("[#6][#7;A;H1X3][#6;A;H2X4][#6;X3]([#8;A;X2H1,X1-])=[O;X1]");	
		return glycinatedRadicalPattern.hasSMARTSPattern(molecule)>0;	
	}	

	public static boolean isPhaseIPolyphenolCandidateOrDerivative(IAtomContainer molecule) throws SMARTSException, CDKException, CloneNotSupportedException{
		boolean ppc = false;
		
		if(!(isGlycosylatedCompound(molecule) || isSulfatedCompound(molecule) || isGlycinatedCompound(molecule))
				&& ChemStructureExplorer.isMetabolizablePolyphenolOrDerivative(molecule)) {
			ppc = true;
		}
		
		return ppc;
		
	}
	
	public static boolean isMixture(IAtomContainer molecule) throws CDKException{
		// compound is not a mixture (checkConnectivity returns 2 or more atomContainers)
		boolean mixture = ConnectivityChecker.partitionIntoMolecules(molecule).getAtomContainerCount()>1;
		return mixture;	
	}
	

	public static boolean isPpsValid(IAtomContainer molecule) throws CDKException{
		// http://eawag-bbd.ethz.ch/predict/notbepredicted.html
		
		InChIGeneratorFactory factory = InChIGeneratorFactory.getInstance();
		InChIGenerator gen1 = factory.getInChIGenerator(molecule);
		String inchikey = gen1.getInchiKey();
		
		boolean valid = (!(isPpsCofactor(inchikey) || isPpsDeadEndCompound(inchikey)))  &&  AtomContainerManipulator.getNaturalExactMass(molecule)<1000.0 && 
				containsCarbon(molecule) && !isMixture(molecule) && 
				( (numberOfAtomWithAtomicNumber(molecule,1) + numberOfAtomWithAtomicNumber(molecule, 6) + 
						numberOfAtomWithAtomicNumber(molecule, 7) + numberOfAtomWithAtomicNumber(molecule, 8) + 
						numberOfAtomWithAtomicNumber(molecule, 15) + numberOfAtomWithAtomicNumber(molecule, 16)
						+ numberOfAtomWithAtomicNumber(molecule, 9) + numberOfAtomWithAtomicNumber(molecule, 17) 
						+ numberOfAtomWithAtomicNumber(molecule, 35) + numberOfAtomWithAtomicNumber(molecule, 53) 
						== molecule.getAtomCount() ) );

		return valid;
	}


	/**
	 * 
	 * @param molecule
	 * @param atomicNumber
	 * @return the number of atoms with the given atomic number 
	 */
	public static int numberOfAtomWithAtomicNumber(IAtomContainer molecule,int atomicNumber){
		int atCount = 0;
		
		for(int l = 0; l <molecule.getAtomCount(); l++){
			if(molecule.getAtom(l).getAtomicNumber() == atomicNumber)
				atCount++;
		}
		return atCount;
	}

	public static boolean isPpsCofactor(String inchikey){
		return (ppsCofactors.containsKey(inchikey));	
	}
	
	public static boolean isPpsCofactor(IAtomContainer molecule) throws CDKException{
		String inchikey = molecule.getProperty("InChIKey");
		if( inchikey != null){
			return (ppsCofactors.containsKey(inchikey));
		} else{
			InChIGeneratorFactory factory = InChIGeneratorFactory.getInstance();
			InChIGenerator gen1 = factory.getInChIGenerator(molecule);
			inchikey = gen1.getInchiKey();
			molecule.setProperty("InChIKey", inchikey);
			molecule.setProperty("InChI", gen1.getInchi());
			return (ppsCofactors.containsKey(inchikey));
		}
			
	}
	
	public static boolean isPpsDeadEndCompound(String inchikey){
		return (ppsDeadEndCompounds.containsKey(inchikey));	
	}
	
	
	public static boolean isPpsDeadEndCompound(IAtomContainer molecule) throws CDKException{
		String inchikey = molecule.getProperty("InChIKey");
		if( inchikey != null){
			return (ppsDeadEndCompounds.containsKey(inchikey));	
		} else{
			InChIGeneratorFactory factory = InChIGeneratorFactory.getInstance();
			InChIGenerator gen1 = factory.getInChIGenerator(molecule);
			inchikey = gen1.getInchiKey();
			molecule.setProperty("InChIKey", inchikey);
			molecule.setProperty("InChI", gen1.getInchi());
			return (ppsDeadEndCompounds.containsKey(inchikey));	
		}
		
		
	}
	


	/**
	 * get isotopemass and alogp
	 * @param molecule
	 * @return
	 * @throws CDKException
	 */
	public static LinkedHashMap<String, String> computePhysicoChemicalProperties(IAtomContainer molecule) throws CDKException {
		LinkedHashMap<String, String> properties = new LinkedHashMap<String, String> ();
		
		ALOGPDescriptor ALogp = new ALOGPDescriptor();
		IDescriptorResult alogp  = ALogp.calculate(molecule).getValue();
		Double majorIsotopeMass = MolecularFormulaManipulator.getMajorIsotopeMass(MolecularFormulaManipulator.getMolecularFormula(molecule));
		properties.put("Major Isotope Mass" , majorIsotopeMass.toString());
		properties.put("ALogP", alogp.toString().split(",")[0]);
		
		return properties;
		
	}
	
	
	
	private static LinkedHashMap<String, String[] > ppsCofactors;
	
	static{
		ppsCofactors = new LinkedHashMap<String, String[]>();
		ppsCofactors.put("ZSLZBFCDCINBPY-NGQYVWNKSA-N", new String[]{"acetyl-CoASH","CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OCC1OC(C(O)C1OP(O)(O)=O)N1C=NC2=C(N)N=CN=C12"});
		ppsCofactors.put("QGZKDVFQNNGYKY-UHFFFAOYSA-N", new String[]{"ammonia","N"});
		ppsCofactors.put("BVKZGUZCCUSVTD-UHFFFAOYSA-M", new String[]{"bicarbonate","OC([O-])=O"});
		ppsCofactors.put("CURLTUGMZLYLDI-UHFFFAOYSA-N", new String[]{"carbon dioxide","O=C=O"});
		ppsCofactors.put("UGFAIRIUMAVXCW-UHFFFAOYSA-N", new String[]{"carbon monoxide","[C-]#[O+]"});
		ppsCofactors.put("BVKZGUZCCUSVTD-UHFFFAOYSA-N", new String[]{"carbonate","OC(O)=O"});
		ppsCofactors.put("RGJOEKWQDUBAIZ-UHFFFAOYSA-N", new String[]{"CoenzymeASH","CC(C)(COP(O)(=O)OP(O)(=O)OCC1OC(C(O)C1OP(O)(O)=O)n1cnc2c(N)ncnc12)C(O)C(=O)NCCC(=O)NCCS"});
		ppsCofactors.put("OTMSDBZUPAUEDD-UHFFFAOYSA-N", new String[]{"ethane","CC"});
		ppsCofactors.put("RWSXRVCMGQZWBV-UHFFFAOYSA-N", new String[]{"glutathione","NC(CCC(=O)NC(CS)C(=O)NCC(O)=O)C(O)=O"});
		ppsCofactors.put("RWSOTUBLDIXVET-UHFFFAOYSA-N", new String[]{"hydrogen sulfide","S"});
		ppsCofactors.put("IOVCWXUNBOPUCH-UHFFFAOYSA-M", new String[]{"nitrite","[O-]N=O"});
		ppsCofactors.put("NBIIXXVUZAFLBC-UHFFFAOYSA-N", new String[]{"phosphate","OP(O)(O)=O"});
		ppsCofactors.put("NBIIXXVUZAFLBC-UHFFFAOYSA-L", new String[]{"phosphate dianion","OP([O-])([O-])=O"});
		ppsCofactors.put("NBIIXXVUZAFLBC-UHFFFAOYSA-K", new String[]{"phosphate trianion","[O-]P([O-])([O-])=O"});
		ppsCofactors.put("QAOWNCQODCNURD-UHFFFAOYSA-L", new String[]{"sulfate","[O-]S([O-])(=O)=O"});
		ppsCofactors.put("LSNNMFCWUKXFEE-UHFFFAOYSA-L", new String[]{"sulfite","[O-]S([O-])=O"});
		ppsCofactors.put("LSNNMFCWUKXFEE-UHFFFAOYSA-N", new String[]{"sulfurous acid","OS(O)=O"});
		ppsCofactors.put("XPRXJFAFJFZWTC-UHFFFAOYSA-M", new String[]{"trioxidosulfate","[O]S([O-])=O"});
	
	}
	
	
	/**
	 * A dictionary with dead-end compounds. Remember that AMBIT does not return stereo-specific configurations
	 */
	private static LinkedHashMap<String, String[] > ppsDeadEndCompounds;

	static{
		ppsDeadEndCompounds = new LinkedHashMap<String, String[]>();
		ppsDeadEndCompounds.put("GTZCVFVGUGFEME-IWQZZHSRSA-M", new String[]{"cis-aconitate","OC(=O)\\C=C(\\CC([O-])=O)C(O)=O"});
		ppsDeadEndCompounds.put("IKHGUXGNUITLKF-UHFFFAOYSA-N", new String[]{"acetaldehyde","CC=O"});
		ppsDeadEndCompounds.put("QTBSBXVTEAMEQO-UHFFFAOYSA-N", new String[]{"acetate","CC(O)=O"});
		ppsDeadEndCompounds.put("WDJHALXBUFZDSR-UHFFFAOYSA-N", new String[]{"acetoacetate","CC(=O)CC(O)=O"});
		ppsDeadEndCompounds.put("CSCPPACGZOOCGX-UHFFFAOYSA-N", new String[]{"acetone","CC(C)=O"});
		ppsDeadEndCompounds.put("HSFWRNGVRCDJHI-UHFFFAOYSA-N", new String[]{"acetylene","C#C"});
		ppsDeadEndCompounds.put("GFFGJBXGBJISGV-UHFFFAOYSA-N", new String[]{"adenine","NC1=C2N=CNC2=NC=N1"});
		ppsDeadEndCompounds.put("WNLRTRBMVRJNCN-UHFFFAOYSA-N", new String[]{"adipate","OC(=O)CCCCC(O)=O"});
		ppsDeadEndCompounds.put("QNAYBMKLOCPYGJ-UHFFFAOYSA-N", new String[]{"alanine","CC(N)C(O)=O"});
		ppsDeadEndCompounds.put("WQZGKKKJIJFFOK-DVKNGEFBSA-N", new String[]{"alpha-D-glucose","OC[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@@H]1O"});
		ppsDeadEndCompounds.put("VUTBELPREDJDDH-UHFFFAOYSA-N", new String[]{"4-amino-5-hydroxymethyl-2-methylpyrimidine","CC1=NC(N)=C(CO)C=N1"});
		ppsDeadEndCompounds.put("TYQCGQRIZGCHNB-JLAZNSOCSA-N", new String[]{"ascorbate","OC[C@H](O)[C@H]1OC(O)=C(O)C1=O"});
		ppsDeadEndCompounds.put("OHJMTUPIZMNBFR-UHFFFAOYSA-N", new String[]{"biuret","NC(=O)NC(N)=O"});
		ppsDeadEndCompounds.put("LRHPLDYGYMQRHN-UHFFFAOYSA-N", new String[]{"1-butanol","CCCCO"});
		ppsDeadEndCompounds.put("FERIUCNNQQJTOY-UHFFFAOYSA-N", new String[]{"butyrate","CCCC(O)=O"});
		ppsDeadEndCompounds.put("KXDHJXZQYSOELW-UHFFFAOYSA-N", new String[]{"carbamate","NC(O)=O"});
		ppsDeadEndCompounds.put("XLJMAIOERFSOGZ-UHFFFAOYSA-N", new String[]{"cyanate","OC#N"});
		ppsDeadEndCompounds.put("XFXPMWWXUTWYJX-UHFFFAOYSA-N", new String[]{"cyanide","[C-]#N"});
		ppsDeadEndCompounds.put("UFULAYFCSOUIOV-UHFFFAOYSA-N", new String[]{"cysteamine","NCCS"});
		ppsDeadEndCompounds.put("XUJNEKJLAYXESH-UHFFFAOYSA-N", new String[]{"cysteine","NC(CS)C(O)=O"});
		ppsDeadEndCompounds.put("OPTASPLRGRRNAP-UHFFFAOYSA-N", new String[]{"cytosine","NC1=NC(=O)NC=C1"});
		ppsDeadEndCompounds.put("WQZGKKKJIJFFOK-GASJEMHNSA-N", new String[]{"D-glucose","OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O"});
		ppsDeadEndCompounds.put("GHVNFZFCNZKVNT-UHFFFAOYSA-N", new String[]{"decanoate","CCCCCCCCCC(O)=O"});
		ppsDeadEndCompounds.put("QSJXEFYPDANLFS-UHFFFAOYSA-N", new String[]{"diacetyl","CC(=O)C(C)=O"});
		ppsDeadEndCompounds.put("LFQSCWFLJHTTHZ-UHFFFAOYSA-N", new String[]{"ethanol","CCO"});
		ppsDeadEndCompounds.put("HZAXFHJVJLSVMW-UHFFFAOYSA-N", new String[]{"ethanolamine","NCCO"});
		ppsDeadEndCompounds.put("LYCAIKOWRPUZTN-UHFFFAOYSA-N", new String[]{"ethylene glycol","OCCO"});
		ppsDeadEndCompounds.put("WSFSSNUMVMOOMR-UHFFFAOYSA-N", new String[]{"formaldehyde","C=O"});
		ppsDeadEndCompounds.put("ZHNUHDYFZUAESO-UHFFFAOYSA-N", new String[]{"formamide","NC=O"});
		ppsDeadEndCompounds.put("BDAGIHXWWSANSR-UHFFFAOYSA-N", new String[]{"formate","OC=O"});
		ppsDeadEndCompounds.put("VZCYOOQTPOCHFL-OWOJBTEDSA-N", new String[]{"fumarate","OC(=O)\\C=C\\C(O)=O"});
		ppsDeadEndCompounds.put("JFCQEDHGNNZCLN-UHFFFAOYSA-N", new String[]{"glutarate","OC(=O)CCCC(O)=O"});
		ppsDeadEndCompounds.put("PEDCQBHIVMGVHV-UHFFFAOYSA-N", new String[]{"glycerol","OCC(O)CO"});
		ppsDeadEndCompounds.put("DHMQDGOQFOQNFH-UHFFFAOYSA-N", new String[]{"glycine","NCC(O)=O"});
		ppsDeadEndCompounds.put("WGCNASOHLSPBMP-UHFFFAOYSA-N", new String[]{"glycoladehyde","O=CCO"});
		ppsDeadEndCompounds.put("AEMRFAOFKBGASW-UHFFFAOYSA-N", new String[]{"glycolate","OCC(O)=O"});
		ppsDeadEndCompounds.put("HHLFWLYXYJOTON-UHFFFAOYSA-N", new String[]{"glyoxylate","OC(=O)C=O"});
		ppsDeadEndCompounds.put("UYTPUPDQBNUYGX-UHFFFAOYSA-N", new String[]{"guanine","NC1=NC2=C(N=CN2)C(=O)N1"});
		ppsDeadEndCompounds.put("AFENDNXGAFYKQO-UHFFFAOYSA-N", new String[]{"2-hydroxybutyrate","CCC(O)C(O)=O"});
		ppsDeadEndCompounds.put("WHBMMWSBFZVSSR-UHFFFAOYSA-N", new String[]{"3-hydroxybutyrate","CC(O)CC(O)=O"});
		ppsDeadEndCompounds.put("SJZRECIVHVDYJC-UHFFFAOYSA-N", new String[]{"4-hydroxybutyrate","OCCCC(O)=O"});
		ppsDeadEndCompounds.put("ALRHLSYJTWAHJZ-UHFFFAOYSA-N", new String[]{"3-hydroxypropanoate","OCCC(O)=O"});
		ppsDeadEndCompounds.put("HHDDCCUIIUWNGJ-UHFFFAOYSA-N", new String[]{"hydroxypyruvate","O=C(O)C(=O)CO"});
		ppsDeadEndCompounds.put("KIAHVTZFUFPMSB-UHFFFAOYSA-N", new String[]{"imidazole pyruvate","CC(=O)C(=O)ON1C=CN=C1"});
		ppsDeadEndCompounds.put("KQNPFQTWMSNSAP-UHFFFAOYSA-N", new String[]{"isobutyrate","CC(C)C(O)=O"});
		ppsDeadEndCompounds.put("GWYFCOCPABKNJV-UHFFFAOYSA-N", new String[]{"isovalerate","CC(C)CC(O)=O"});
		ppsDeadEndCompounds.put("TYEYBOSBBBHJIV-UHFFFAOYSA-N", new String[]{"2-ketobutyrate","CCC(=O)C(=O)O"});
		ppsDeadEndCompounds.put("TYEYBOSBBBHJIV-UHFFFAOYSA-N", new String[]{"L-2-oxobutyrate","CCC(=O)C(O)=O"});
		ppsDeadEndCompounds.put("JVTAAEKCZFNVCJ-UHFFFAOYSA-N", new String[]{"lactate","CC(O)C(O)=O"});
		ppsDeadEndCompounds.put("POULHZVOKOAJMA-UHFFFAOYSA-N", new String[]{"lauric acid","CCCCCCCCCCCC(O)=O"});
		ppsDeadEndCompounds.put("BJEPYKJPYRNKOW-UHFFFAOYSA-N", new String[]{"malate","OC(CC(O)=O)C(O)=O"});
		ppsDeadEndCompounds.put("FSQQTNAZHBEJLS-UPHRSURJSA-N", new String[]{"maleamate","NC(=O)\\C=C/C(O)=O"});
		ppsDeadEndCompounds.put("VZCYOOQTPOCHFL-UPHRSURJSA-N", new String[]{"maleate","OC(=O)\\C=C/C(O)=O"});
		ppsDeadEndCompounds.put("OFOBLEOULBTSOW-UHFFFAOYSA-N", new String[]{"malonate","C(C(=O)O)C(=O)O"});
		ppsDeadEndCompounds.put("OFOBLEOULBTSOW-UHFFFAOYSA-N", new String[]{"malonate semialdehyde","OC(=O)CC(O)=O"});
		ppsDeadEndCompounds.put("VNWKTOKETHGBQD-UHFFFAOYSA-N", new String[]{"methane","C"});
		ppsDeadEndCompounds.put("OKKJLVBELUTLKV-UHFFFAOYSA-N", new String[]{"methanol","CO"});
		ppsDeadEndCompounds.put("BAVYZALUXZFZLV-UHFFFAOYSA-N", new String[]{"methylamine","CN"});
		ppsDeadEndCompounds.put("YYPNJNDODFVZLE-UHFFFAOYSA-N", new String[]{"3-methylcrotonate","CC(C)=CC(O)=O"});
		ppsDeadEndCompounds.put("HNEGQIOMVPPMNR-IHWYPQMZSA-N", new String[]{"2-methylmaleate","C\\C(=C\\C(O)=O)C(O)=O"});
		ppsDeadEndCompounds.put("ZIYVHBGGAOATLY-UHFFFAOYSA-N", new String[]{"methylmalonate","CC(C(O)=O)C(O)=O"});
		ppsDeadEndCompounds.put("PVNIIMVLHYAWGP-UHFFFAOYSA-N", new String[]{"niacin","OC(=O)C1=CC=CN=C1"});
		ppsDeadEndCompounds.put("WWZKQHOCKIZLMA-UHFFFAOYSA-N", new String[]{"ocatanoate","CCCCCCCC(=O)O"});
		ppsDeadEndCompounds.put("MUBZPKHOEPUJKR-UHFFFAOYSA-N", new String[]{"oxalate","OC(=O)C(O)=O"});
		ppsDeadEndCompounds.put("KHPXUQMNIQBQEV-UHFFFAOYSA-N", new String[]{"oxaloacetate","OC(=O)CC(=O)C(O)=O"});
		ppsDeadEndCompounds.put("KPGXRSRHYNQIFN-UHFFFAOYSA-N", new String[]{"2-oxoglutarate","OC(=O)CCC(=O)C(O)=O"});
		ppsDeadEndCompounds.put("NOXRYJAWRSNUJD-UHFFFAOYSA-N", new String[]{"2-oxopent-4-enoate","OC(=O)C(=O)CC=C"});
		ppsDeadEndCompounds.put("ATUOYWHBWRKTHZ-UHFFFAOYSA-N", new String[]{"propane","CCC"});
		ppsDeadEndCompounds.put("KFZMGEQAYNKOFK-UHFFFAOYSA-N", new String[]{"2-propanol","CC(C)O"});
		ppsDeadEndCompounds.put("XBDQKXXYIPTUBI-UHFFFAOYSA-N", new String[]{"propionate","CCC(O)=O"});
		ppsDeadEndCompounds.put("DNIAPMSPPWPWGF-UHFFFAOYSA-N", new String[]{"propylene glycol","CC(O)CO"});
		ppsDeadEndCompounds.put("LCTONWCANYUPML-UHFFFAOYSA-N", new String[]{"pyruvate","CC(=O)C(O)=O"});
		ppsDeadEndCompounds.put("FSYKKLYZXJSNPZ-UHFFFAOYSA-N", new String[]{"sarcosine","CNCC(O)=O"});
		ppsDeadEndCompounds.put("KDYFGRWQOYBRFD-UHFFFAOYSA-N", new String[]{"succinate","OC(=O)CCC(O)=O"});
		ppsDeadEndCompounds.put("AYFVYJQAPQTCCC-GBXIJSLDSA-N", new String[]{"threonine","C[C@@H](O)[C@H](N)C(O)=O"});
		ppsDeadEndCompounds.put("RWQNBRDOKXIBIV-UHFFFAOYSA-N", new String[]{"thymine","CC1=CNC(=O)NC1=O"});
		ppsDeadEndCompounds.put("ISAKRJDGNUQOIC-UHFFFAOYSA-N", new String[]{"uracil","O=C1NC=CC(=O)N1"});
		ppsDeadEndCompounds.put("XSQUKJJJFZCRTK-UHFFFAOYSA-N", new String[]{"urea","NC(N)=O"});
		ppsDeadEndCompounds.put("LRFVTYWOQMYALW-UHFFFAOYSA-N", new String[]{"xanthine","O=C1NC2=C(NC=N2)C(=O)N1"});

	}
	
		
	

	
	
	
	
}

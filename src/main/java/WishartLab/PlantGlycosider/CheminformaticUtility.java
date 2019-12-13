package WishartLab.PlantGlycosider;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.Stack;

import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.graph.ConnectivityChecker;
import org.openscience.cdk.inchi.InChIGenerator;
import org.openscience.cdk.inchi.InChIGeneratorFactory;
import org.openscience.cdk.inchi.InChIToStructure;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.io.SDFWriter;
import org.openscience.cdk.layout.StructureDiagramGenerator;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import ambit2.smarts.query.SMARTSException;
import ambit2.smarts.query.SmartsPatternCDK;
import net.sf.jniinchi.INCHI_RET;

public class CheminformaticUtility {
	
	public ArrayList<Integer[]> combination_list = new ArrayList<Integer[]>(); // this list constantly changes
	

	/**
	 * getter function
	 * @return
	 */
	public ArrayList<Integer[]> get_combination_list() {
		return this.combination_list;
	}
	
	
	/**
	 * get all possible set of sugar 
	 * you can only call it after permutations() is called and store everything into combination_list;
	 * @param input
	 * @param factor
	 * @return
	 */
	public ArrayList<Integer[]> GetAllPossibleSet(Integer[] input, int factor) {
		ArrayList<Integer[]> output = new ArrayList<Integer[]>();
		ArrayList<String> output_for_check_dup = new ArrayList<String>();
		Set<Integer> set = new HashSet<Integer>();
		for(int i = 0; i < input.length; i++) {
			set.add(input[i]);
		}
		
		permutations(set, new Stack<Integer>(), set.size());
		
		for(int f =  1; f <= factor; f++) {
			for (int i = 0; i < combination_list.size(); i++) {
				ArrayList<Integer[]> temp = GetPermutation(combination_list.get(i), f);
				for (int n = 0; n < temp.size(); n++) {
					if (output_for_check_dup.contains(Arrays.toString(temp.get(n)))){
						continue;
					}else {
						output.add((Integer[])temp.get(n));
						output_for_check_dup.add(Arrays.toString(temp.get(n)));
					}
					
				}
			}
		}
		
		return output;
	}

	/**
	 * given list and factor (number of item in set), generate all combination)
	 * 
	 * @param input
	 * @param factor
	 * @return
	 */
	public ArrayList<Integer[]> GetPermutation(Integer[] input, int factor) {
		ArrayList<Integer[]> output = new ArrayList<Integer[]>();

		int[] s = new int[factor]; // here we'll keep indices
		// pointing to elements in input array

		if (factor <= input.length) {
			// first index sequence: 0, 1, 2, ...
			for (int i = 0; (s[i] = i) < factor - 1; i++)
				;
			output.add(getSubset(input, s));
			for (;;) {
				int i;
				// find position of item that can be incremented
				for (i = factor - 1; i >= 0 && s[i] == input.length - factor + i; i--)
					;
				if (i < 0) {
					break;
				}
				s[i]++; // increment this item
				for (++i; i < factor; i++) { // fill up remaining items
					s[i] = s[i - 1] + 1;
				}
				output.add(getSubset(input, s));
			}
		}

		return output;

	}

	
	/**
	 * Do permutation
	 * @param items
	 * @param permutation
	 * @param size
	 */
	public void permutations(Set<Integer> items, Stack<Integer> permutation, int size) {
		
		/* permutation stack has become equal to size that we require */
		if (permutation.size() == size) {
			/* print the permutation */
			combination_list.add(permutation.toArray(new Integer[0]));
		}

		/* items available for permutation */
		Integer[] availableItems = items.toArray(new Integer[0]);
		for (Integer i : availableItems) {
			
			/* add current item */
			permutation.push(i);

			/* remove item from available item set */
			items.remove(i);

			/* pass it on for next permutation */
			permutations(items, permutation, size);

			/* pop and put the removed item back */
			items.add(permutation.pop());
		}
		
	}
	
	

	
	/**
	 * Different than TransformSmilesToContainer because this function don't consider the sugar attach position
	 * @param smiles
	 * @return
	 */
	public static IAtomContainer parseSmilesToContainer(String smiles) {
		
		
		IChemObjectBuilder builder = SilentChemObjectBuilder.getInstance();
		SmilesParser smiParser = new SmilesParser(builder);
		IChemObjectBuilder containerBuilder = SilentChemObjectBuilder.getInstance();
		IAtomContainer mole = containerBuilder.newAtomContainer();
		
		try {
			
			mole = smiParser.parseSmiles(smiles);
			
		} catch (InvalidSmilesException e) {
			
			mole = null;
		}
		
		return mole;
		
	}
	
	
	/**
	 * write the IAtomContainer to SDF
	 * @param uniqueContainer
	 * @param compound_inchikey
	 */
	public static void saveIAtomContainerSetToSDF(IAtomContainerSet uniqueContainer, String compound_inchikey) {
		String current_dir = System.getProperty("user.dir");
		try {
			FileWriter fw = new FileWriter(String.format("%s/generatedSDF/%s.sdf", current_dir, compound_inchikey), true);
			SDFWriter sdfwriter = new SDFWriter(fw);			        
	        
	        for(int i = 0; i < uniqueContainer.getAtomContainerCount(); i++) {
				IAtomContainer mole = CheminformaticUtility.Generate2DCoordinate(uniqueContainer.getAtomContainer(i));
				if (mole != null) {
					sdfwriter.write(AtomContainerManipulator.removeHydrogens(mole));
				}	
			}
	        
			sdfwriter.close();
	        fw.close();
		} catch (IOException e) {
			e.printStackTrace();
		} catch (CDKException e) {
			e.printStackTrace();
		}
        
	}
	
	
	
	
	/**
	 * transform smiles to IAtomContainer based on sugar smiles and sugar_index
	 * @param sugar_name
	 * @param suger_index
	 * @return
	 * @throws InvalidSmilesException
	 */
	public static IAtomContainer TransformSugarSmilesToContainer(String sugar_smiles, int sugar_index) throws InvalidSmilesException {
		
		SmilesParser sp = new SmilesParser(SilentChemObjectBuilder.getInstance());
		IAtomContainer sugar_mole = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainer.class);
		sugar_mole = sp.parseSmiles(sugar_smiles);
		sugar_mole.getAtom(sugar_index).setProperty("index", "sugar");
		
		return sugar_mole;
	}
	
	
	
	/**
	 * get list of OH by matching the smarts 
	 * @param mole
	 * @return
	 * @throws SMARTSException
	 */
	public static Set<Integer> GetOHList(IAtomContainer mole) throws SMARTSException {
		Set<Integer> output = new HashSet<Integer>();
		SmartsPatternCDK smarts = new SmartsPatternCDK();
		smarts.setSmarts("[OH]");
		
		if(smarts.hasSMARTSPattern(mole) != 0) {
			List<List<Integer>> matach_atom_ind = smarts.getUniqueMatchingAtoms(mole); // two matching group will be list.size() = 2
			for(int i = 0; i< matach_atom_ind.size(); i++) {
				List<Integer> tmp = matach_atom_ind.get(i);			// usually the first index is the som;
				int som = tmp.get(0);
				output.add(som);
			}
		}
		
		
		return output;
	}
	
	/**
	 * 
	 * @param mole
	 * @return String[0] inchi, String[1] inchikey
	 * @throws CDKException
	 */
	public static String[] GenerateInChi(IAtomContainer mole) throws CDKException {
		InChIGeneratorFactory factory = InChIGeneratorFactory.getInstance();
		InChIGenerator gen = factory.getInChIGenerator(mole);

		INCHI_RET ret = gen.getReturnStatus();
		if (ret == INCHI_RET.WARNING) {
			// InChI generated, but with warning message
			System.out.println("InChI warning: " + gen.getMessage());
			return null;
		} else if (ret != INCHI_RET.OKAY) {
			// InChI generation failed
			return null;
			// throw new CDKException("InChI failed: " + ret.toString() + " [" +
			// gen.getMessage() + "]");

		}

		String[] return_value = new String[] { gen.getInchi(), gen.getInchiKey() };
		return return_value;
	}
	
	
	
	/**
	 * split OH out of container and return the longest one 
	 * 
	 * @param mole
	 * @return
	 * @throws CDKException
	 */
	public static String SplitContainer(String mole_smiles) throws CDKException {		
		
		String[] splited_mol = mole_smiles.split("\\.");
		String longest = "";
		
		if(splited_mol.length == 1) {
			longest = mole_smiles;
		}else if (splited_mol.length > 1) {
			for(int i = 0; i < splited_mol.length; i++) {
				if(splited_mol[i].length() > longest.length()) {
					longest = splited_mol[i];
				}
			}
		}
		
		return longest;
	}
	
	
	/**
	 * split OH out of container and return the longest one 
	 * based on compare number of atoms
	 * @param mole
	 * @return
	 * @throws CDKException
	 */
	public static IAtomContainer SplitContainer(IAtomContainer mole) throws CDKException {		
		
		IAtomContainerSet moleset = ConnectivityChecker.partitionIntoMolecules(mole);
		int num_container = moleset.getAtomContainerCount();
		if(num_container == 1) {
			return moleset.getAtomContainer(0);
		}
		
		int highest_count = 0;
		int index = 0;
		for(int i = 0; i < num_container; i++) {
			int num_atom = moleset.getAtomContainer(i).getAtomCount();
			if(moleset.getAtomContainer(i).getAtomCount() > highest_count) {
				highest_count = num_atom;
				index = i;
			}
		}
		
		return moleset.getAtomContainer(index);
	}
	
	/**
	 * generate 2d coordinate by giving molecules
	 * @param mole
	 * @return
	 */
	public static IAtomContainer Generate2DCoordinate(IAtomContainer mole) {

		AtomContainerManipulator.suppressHydrogens( mole);
		AtomContainerManipulator.convertImplicitToExplicitHydrogens(mole);
		StructureDiagramGenerator sdg = new StructureDiagramGenerator();
		sdg.setMolecule(mole);

		try {
			sdg.generateCoordinates();
			IAtomContainer new_mole = sdg.getMolecule();
			return new_mole;
		} catch (CDKException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			return null;
		}
		
		
	}
	
	/**
	 * get iatomcontainer from inchi
	 * @param inchi
	 * @return
	 * @throws Exception
	 */
	public static IAtomContainer getAtomContainerFromInChI(String inchi) throws Exception {
		InChIGeneratorFactory inchiFactory = InChIGeneratorFactory.getInstance();
		InChIToStructure its = inchiFactory.getInChIToStructure(inchi, DefaultChemObjectBuilder.getInstance());
		if(its == null) {
			throw new Exception("InChI problem: " + inchi);
		}
		INCHI_RET ret = its.getReturnStatus();
		if (ret == INCHI_RET.WARNING) {
		} else if (ret != INCHI_RET.OKAY) {
			throw new Exception("InChI problem: " + inchi);
		}
		return its.getAtomContainer();
	}
	
	
	/**
	 * *helper function 
	 * generate actual subset by index sequence
	 * @param input
	 * @param subset
	 * @return
	 */
	public Integer[] getSubset(Integer[] input, int[] subset) {
		Integer[] result = new Integer[subset.length];
		for (int i = 0; i < subset.length; i++)
			result[i] = input[subset[i]];
		return result;
	}
	
	
	
	
	

}

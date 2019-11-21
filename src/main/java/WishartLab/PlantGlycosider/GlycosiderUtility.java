package WishartLab.PlantGlycosider;

import java.util.ArrayList;
import java.util.List;

import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.inchi.InChIGenerator;
import org.openscience.cdk.inchi.InChIGeneratorFactory;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

public class GlycosiderUtility {
	
	
//	public static void saveIAtomContainerToSDF(IAtomContainer uniqueContainer) {
//		String current_dir = System.getProperty("user.dir");
//		try {
//			FileWriter fw = new FileWriter(String.format("%s/generatedfolder/test.sdf", current_dir), true);
//			SDFWriter sdfwriter = new SDFWriter(fw);			        
//	        
//			sdfwriter.write(AtomContainerManipulator.removeHydrogens(uniqueContainer));
//			
//			sdfwriter.close();
//	        fw.close();
//		} catch (IOException e) {
//			e.printStackTrace();
//		} catch (CDKException e) {
//			e.printStackTrace();
//		}
//        
//	}
	
	/**
	 * combine two IAtomContainer, they must have charge on the atom of interests
	 * atom is not mapped when comes in. Reference: https://forum.knime.com/t/cdk-atom-numbering/4431/4
	 * order: add sugar to mole
	 * mole_index.length should be same as sugar_set.size()
	 * @param mole don't contain chargesm but the the charge info is in charge_table 
	 * @param mole_index the index for mole to add sugar
	 * @param sugar_set Already contain the charges 
	 * @return
	 * @throws CDKException 
	 */
	public static IAtomContainer CombineContainers(IAtomContainer mole, int[] mole_index, 
			IAtomContainerSet sugar_set) throws CDKException {
		for(int c = 0; c < mole_index.length; c++) {
			// for each iteration; add property to compound's O(-H) atom
			// since the sugar_set.getAtomContainer(n) already contain the property sugar
			// remove the mole and sugar property, so property will be reset in each iteration
			mole.getAtom(mole_index[c]).setProperty("index", "mole");
			
			IAtomContainer sugar = sugar_set.getAtomContainer(c);
			AtomContainerManipulator.convertImplicitToExplicitHydrogens(mole);
			AtomContainerManipulator.convertImplicitToExplicitHydrogens(sugar);
			
			mole.add(sugar);
			
			int original_atom_index = new Integer(0);
			int sugar_atom_index = new Integer(0);
			int sugar_atom_remove_index = new Integer(0);
			for(int k = 0; k < mole.getAtomCount();k++) {

				if (mole.getAtom(k).getProperty("index") != null) {
					if(mole.getAtom(k).getProperty("index").equals("mole")) {
						original_atom_index = k;
					}else if(mole.getAtom(k).getProperty("index").equals("sugar")) {
						// find the nearest carbon for sugar; since mole is uncertain
						sugar_atom_index = k;
						List<IAtom> nearestAtom = mole.getConnectedAtomsList(mole.getAtom(k));
						for (int i = 0; i < nearestAtom.size(); i++) {
							if(nearestAtom.get(i).getSymbol().equals("C")) {
								sugar_atom_remove_index = mole.indexOf(nearestAtom.get(i));
							}
						}

					}else {
						continue;
					}

				}
			}
			
			

			mole.addBond(original_atom_index, sugar_atom_remove_index, IBond.Order.SINGLE);
			mole.removeBond(mole.getAtom(sugar_atom_index), mole.getAtom(sugar_atom_remove_index));
			
			// remove extra valence on atom O
			List<IAtom> nearestAtom = mole.getConnectedAtomsList(mole.getAtom(original_atom_index));
			for (int i = 0; i < nearestAtom.size(); i++) {
				//System.out.println(nearestAtom.get(i).getSymbol());
				if(nearestAtom.get(i).getSymbol().equals("H")) {
					mole.removeBond(nearestAtom.get(i), mole.getAtom(original_atom_index));
				}
			}
			
			// remove all properties
			for(int k = 0; k < mole.getAtomCount();k++) {
				if (mole.getAtom(k).getProperty("index") != null) {
					mole.getAtom(k).removeProperty("index");
				}
			}
			
		}
		return mole;
	}
	
	/**
	 * remove duplication based on the inchikey
	 * @param transformedSmiles
	 * @return
	 * @throws CDKException 
	 */
	public static IAtomContainerSet RemoveDupliation(IAtomContainerSet transformedSmiles) throws CDKException {
		ArrayList<String> non_dup = new ArrayList<String>();
		IChemObjectBuilder containerSetBuilder = SilentChemObjectBuilder.getInstance();
		IAtomContainerSet non_dup_container = containerSetBuilder.newInstance(IAtomContainerSet.class);
		InChIGeneratorFactory factory = InChIGeneratorFactory.getInstance();
		
		for(int i = 0; i < transformedSmiles.getAtomContainerCount(); i++) {
			try {
				IAtomContainer mole = transformedSmiles.getAtomContainer(i);
				
				InChIGenerator gen = factory.getInChIGenerator(mole);
				String inchikey = gen.getInchiKey();
				if (non_dup.contains(inchikey)) {
					continue;
				}else {
					non_dup.add(inchikey);
					non_dup_container.addAtomContainer(mole);
				}
			} catch (InvalidSmilesException e) {
				
				e.printStackTrace();
				continue;
			} catch (CDKException e) {
				
				e.printStackTrace();
				continue;
			}
		}
		
		return non_dup_container;
	}
	
	
	/**
	 * remove duplication based on the inchikey
	 * @param transformedSmiles
	 * @return
	 * @throws CDKException 
	 */
	public static IAtomContainerSet RemoveDupliation(ArrayList<String> transformedSmiles) throws CDKException {
		ArrayList<String> non_dup = new ArrayList<String>();
		IChemObjectBuilder containerSetBuilder = SilentChemObjectBuilder.getInstance();
		IAtomContainerSet non_dup_container = containerSetBuilder.newInstance(IAtomContainerSet.class);
		InChIGeneratorFactory factory = InChIGeneratorFactory.getInstance();
		IChemObjectBuilder builder = SilentChemObjectBuilder.getInstance();
		SmilesParser smiParser = new SmilesParser(builder);
		
		for(int i = 0; i < transformedSmiles.size(); i++) {
			try {
				IAtomContainer mole = smiParser.parseSmiles(transformedSmiles.get(i));
				
				InChIGenerator gen = factory.getInChIGenerator(mole);
				String inchikey = gen.getInchiKey();
				if (non_dup.contains(inchikey)) {
					continue;
				}else {
					non_dup.add(inchikey);
					non_dup_container.addAtomContainer(mole);
				}
			} catch (InvalidSmilesException e) {
				
				e.printStackTrace();
				continue;
			} catch (CDKException e) {
				
				e.printStackTrace();
				continue;
			}
		}
		
		return non_dup_container;
	}
	
	
}

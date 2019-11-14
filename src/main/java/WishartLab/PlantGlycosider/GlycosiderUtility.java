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

public class GlycosiderUtility {
	
	
	/**
	 * combine two IAtomContainer, they must have charge on the atom of interests
	 * 
	 * order: add sugar to mole
	 * @param mole
	 * @param sugar
	 * @return
	 * @throws CDKException 
	 */
	public static IAtomContainer CombineContainers(IAtomContainer mole, double[] charge_table, IAtomContainerSet sugar) throws CDKException {

		for(int c = 0; c < charge_table.length; c++) {
			mole.add(sugar.getAtomContainer(c));
			int index_1 = new Integer(0);
			int index_2 = new Integer(0);
			double charge = charge_table[c];
			for (int k = 0; k < mole.getAtomCount(); k++) {
				
				IAtom atoms = mole.getAtom(k);
				if (atoms.getCharge() != null) {					
					// find sugar molecule mark
					// with garantee of number of sugar can't exceed the number of OH, add sugar one by one
					// need to charge for order
					if (atoms.getCharge() == 1.0) {
						index_1 = k;
						atoms.setCharge(null);
					} 
					
					// find original molecule mark
					else if (atoms.getCharge() == charge) {
						List<IAtom> nearestAtom = mole.getConnectedAtomsList(atoms);
						for (int i = 0; i < nearestAtom.size(); i++) {
							if(nearestAtom.get(i).getSymbol().equals("C")) {
								int real_connect = mole.indexOf(nearestAtom.get(i));
								index_2 = real_connect;
								mole.removeBond(atoms, nearestAtom.get(i));
							}
						}
						
						atoms.setCharge(null);
					}
				}
				
				
			}
			mole.addBond(index_1, index_2, IBond.Order.SINGLE);
			for (int k = 0; k < mole.getAtomCount(); k++) {
				IAtom atoms = mole.getAtom(k);
				if (atoms.getCharge() != null) {
					if(atoms.getCharge() == 1.0) {
						atoms.setCharge(null);
					}
					
				}
			}
			
		}
		
		// remove the charges
		for (int k = 0; k < mole.getAtomCount(); k++) {
			IAtom atoms = mole.getAtom(k);
			
			if (atoms.getCharge() != null) {
				atoms.setCharge(null);
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

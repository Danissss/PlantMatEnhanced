package WishartLab.PlantGlycosider;

import org.openscience.cdk.graph.ConnectivityChecker;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;

/**
 * Pre validation only validate two things: 
 * Is the compound good for glycosidation? 
 * Is there any constraint for this particular compound?
 * 
 * @author xuan
 *
 */
public class PreValidation {
	
	
	public boolean isCompoundValidate(IAtomContainer mole) {
		boolean isValidate = true;
		
		// check for mixture
		if(isMixture(mole)) {
			return false;
		}
		
		// check for organic
		if(!containsCarbon(mole)) {
			return false;
		}
		
		
		
		
		
		return isValidate;
		
	}
	
	
	/**
	 * return true if mole is mixture
	 * @param mole
	 * @return
	 */
	private boolean isMixture(IAtomContainer mole) {
		boolean is_mixture = ConnectivityChecker.partitionIntoMolecules(mole).getAtomContainerCount()>1;
		return is_mixture;
		
	}
	
	/**
	 * 
	 * @param mole
	 * @return
	 */
	private boolean containsCarbon(IAtomContainer mole) {
		boolean carbon = false;
		for(IAtom at : mole.atoms()){
			if(at.getAtomicNumber() == 6){
				carbon = true;
				break;
			}
		}	
		return carbon;
	}
}

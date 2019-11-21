package WishartLab.PlantGlycosider;

import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IChemObjectBuilder;

public class PostValidation {
	
	public IAtomContainerSet validateCompound(IAtomContainerSet moleset) {
		IChemObjectBuilder builder = DefaultChemObjectBuilder.getInstance();
		IAtomContainerSet result_mole = builder.newInstance(IAtomContainerSet.class);
		
		for(int i = 0; i < moleset.getAtomContainerCount(); i++) {
			IAtomContainer mole = moleset.getAtomContainer(i);
			if(isValidate(mole)) {
				result_mole.addAtomContainer(mole);
			}
		}
		
		return result_mole;
	}
	
	
	/**
	 * validation logic
	 * @param mole
	 * @return
	 */
	private boolean isValidate(IAtomContainer mole) {
		boolean is_validate = true;
		
		
		
		return is_validate;
	}
}

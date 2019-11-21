package WishartLab.PlantGlycosider;

import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;

import ambit2.smarts.query.SMARTSException;

/**
 * Pre validation only validate two things: 
 * Is the compound good for glycosidation? 
 * Is there any constraint for this particular compound?
 * 
 * @author xuan
 *
 */
public class PreValidation {
	
	/**
	 * reject mixture
	 * reject inorganic
	 * reject more than 1000.0 mass
	 * reject any compounds that include glucoside, sulfate and glycinate compound/attachments
	 * 
	 * @param mole
	 * @return
	 */
	public boolean isCompoundValidate(IAtomContainer mole) {
		boolean isValidate = false;

		try {
			boolean is_ppsValid = ChemStructureExplorer.isPPSValid(mole); // return true if is pps valid; false otherwise
			boolean is_pure_sec_metabolites = ChemStructureExplorer.isPureSecondaryMetoblite(mole);
			boolean is_polyphenol = ChemStructureExplorer.isPolyphenolOrDerivative(mole);
			System.out.println(String.format("is_ppsValid = %b; "
					+ "\nis_pure_sec_metabolites = %b; "
					+ "\nis_polyphenol = %b; ", is_ppsValid,is_pure_sec_metabolites,is_polyphenol));
			if(is_ppsValid && is_pure_sec_metabolites && is_polyphenol) {
				isValidate = true;
			}
		} catch (CDKException e) {
			e.printStackTrace();
		} catch (SMARTSException e) {
			e.printStackTrace();
		} catch (CloneNotSupportedException e) {
			e.printStackTrace();
		}
		
		return isValidate;
		
	}
}

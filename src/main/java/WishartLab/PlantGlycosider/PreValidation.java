package WishartLab.PlantGlycosider;

import java.util.Set;

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
	 * only accept polyphenol
	 * phenol class
	 * Flavonoids stand as the first class of polyphenols; They can also be divided into three groups: anthocyanins, flavones and flavonols.
	 * Tannin is the name derived from French “Tanin” (tanning substance) and used for a range of natural polyphenols.; very high molecule weight
	 * secondary metabolites extracted from plants are subdivided in three major classes: terpenoids, alkaloids and phenolics.
	 * is_alkaloid validation will be a fingerprint based machine learning model
	 * https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-019-2963-6
	 * @param mole
	 * @return
	 */
	public boolean isCompoundValidate(IAtomContainer mole) {
		boolean isValidate = false;

		try {
			boolean is_ppsValid = ChemStructureExplorer.isPPSValid(mole); // return true if is pps valid; false otherwise
			boolean is_pure_sec_metabolites = ChemStructureExplorer.isPureSecondaryMetoblite(mole);
			boolean is_polyphenol = ChemStructureExplorer.isPolyphenolOrDerivative(mole);
			boolean is_reasonable_oh = isReasonableOH(mole);
			boolean is_terpene = ChemStructureExplorer.isTerpenoid(mole);
			System.out.println(String.format("is_ppsValid = %b; "
					+ "\nis_pure_sec_metabolites = %b; "
					+ "\nis_polyphenol = %b; ", is_ppsValid,is_pure_sec_metabolites,is_polyphenol));
			if(is_ppsValid && is_pure_sec_metabolites && is_reasonable_oh) {
				if(is_polyphenol || is_terpene) {
					isValidate = true;
				}
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
	
	private boolean isReasonableOH(IAtomContainer mole) throws SMARTSException {
		boolean is_reasonable = true;
		Set<Integer> OH_list = CheminformaticUtility.GetOHList(mole); // return list of OH
		
		if (OH_list.size() > 4) {
			is_reasonable = false;
		}
		
		return is_reasonable;
	}
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
}

package WishartLab.PlantGlycosider;

import java.io.FileNotFoundException;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IChemObjectBuilder;

public class PostValidation {
	
	
	/**
	 * Glycosylation allows the conversion of water-insoluble and unstable organic compounds into the corresponding 
	 * water-soluble and stable ones to improve their bio- and pharmacological properties, Enzymatic glycosylation of terpenoids Francisco Rivas
	 * https://link.springer.com/article/10.1007/s11101-013-9301-9
	 * 
	 * check for water soluable https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3965570/
	 * logD logP and Lipophilicity: https://www.cambridgemedchemconsulting.com/resources/physiochem/logD.html
	 * hydrophilic ("water-loving") substances tend to dissolve in water and other hydrophilic substances.
	 * partition coefficient (P) or distribution coefficient (D): https://en.wikipedia.org/wiki/Partition_coefficient
	 *  hydrophilic ("water-loving") or hydrophobic ("water-fearing") 
	 * @param moleset
	 * @return
	 * @throws Exception 
	 * @throws FileNotFoundException 
	 */
	public IAtomContainerSet validateCompound(IAtomContainerSet moleset) throws FileNotFoundException, Exception {
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
	 * @throws Exception 
	 * @throws FileNotFoundException 
	 */
	private boolean isValidate(IAtomContainer mole) throws FileNotFoundException, Exception {
		// Note The code assumes that aromaticity has been detected before evaluating this descriptor. 
		// The code also expects that the molecule will have hydrogens explicitly set. For SD files, 
		// this is usually not a problem since hydrogens are explicit. But for the case of molecules 
		// obtained from SMILES, hydrogens must be made explicit.
		// return three values;
		// 1.ALogP - Ghose-Crippen LogKow
		// 2.ALogP2
		// 3.amr - molar refractivity
		CreateSolubilityModel logSmodel = new CreateSolubilityModel();
		double logs = logSmodel.predictLogS(mole);
		boolean is_LogSvalidate = logSmodel.isReasonableLogS(logs);
		
		CreateLogDModel logDmodel = new CreateLogDModel();
		double logd = logDmodel.predictLogD(mole);
		boolean is_LogDvalidate = logDmodel.isReasonableLogD(logd);
		
		return (is_LogSvalidate&&is_LogDvalidate);
	}
	
	/**
	 * Solubility information: displays intrinsic solubility value, solubility at pH 7.4 and a qualitative category for the predicted solubility. These categories are:
	 * low: if solubility is < 0.01 mg/ml
	 * moderate: if solubility is between 0.01 and 0.06 mg/ml
	 * high: if solubility is > 0.06 mg/ml
	 * 
	 * 
	 * float vs double java
	 * float is represented in 32 bits, with 1 sign bit, 8 bits of exponent, 
	 * double is represented in 64 bits, with 1 sign bit, 11 bits of exponent, and 52 bits of significant.
	 * 
	 * via machine learning  https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3965570/#notes-1
	 * or via chemaxon (purely based on structure info) https://apidocs.chemaxon.com/jchem/doc/dev/java/api/
	 * @return
	 */
//	private double getSoluabilityByFragment() {
	// for chemaxon api https://pro.chemicalize.com/app/calculations/api/calculation
	// use wishartlab chemaxon api key.
//		double soluable = 0.0;
//		SolubilityCalculator calculator = new SolubilityCalculator();
//		// calculate intrinsic solubility
//		SolubilityResult result1 = calculator.calculateIntrinsicSolubility(mol);
//		double sol1 = result.getSolubility();                            // intrinsic solubility
//		String category1 = result.getSolubilityCategory().shortName();   // intrinsic solubility category
//
//		// calculate pH-dependent solubility
//		SolubilityResult result2 = calculator.calculatePhDependentSolubility(mol, 7.4);
//		double sol2 = result2.getSolubility();                            // solubility at pH 7.4
//		String category2 = result2.getSolubilityCategory().shortName();   // solubility category at pH 7.4
//		
//		return soluable;
//	}
//	
	
}














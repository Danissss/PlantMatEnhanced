/**
 * 
 * @author Djoumbou Feunang, Yannick, PhD
 *
 */

package WishartLab.PlantGlycosider;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.Map;

import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import ambit2.smarts.query.SMARTSException;
import ambit2.smarts.query.SmartsPatternCDK;

public class ChemicalClassFinder {
	
	public ChemicalClassFinder() {
		// TODO Auto-generated constructor stub
	}

	public enum ChemicalClassName {
		ALPHA_HYDROXYLATED_FATTY_ACID, BETA_HYDROXYLATED_FATTY_ACID,		
		BILE_ACID, ETHER_LIPID, 
		FATTY_ACID,
		GLYCEROLIPID, GLYCEROPHOSPHOLIPID, 
		OMEGA_HYDROXYLATED_FATTY_ACID,
		SPHINGOLIPID, 
		UNSATURATED_FATTY_ACID,
		UNFUNCTIONALIZED_UNSATURATED_FATTY_ACID,
		SATURATED_FATTY_ACID,
		UNFUNCTIONALIZED_SATURATED_FATTY_ACID,
		GLYCEROL_3_PHOSPHATE_INOSITOL,
		C23_BILE_ACID, C24_BILE_ACID,
		SULFATE_ESTER
		
	}
	
	public static LinkedHashMap<ChemicalClassName, LinkedHashMap<String, String[]>>			chemicalClassDefinitions;
	static {
		chemicalClassDefinitions = new LinkedHashMap<ChemicalClassName, LinkedHashMap<String, String[]>>();
		
		chemicalClassDefinitions.put(ChemicalClassName.ALPHA_HYDROXYLATED_FATTY_ACID, new LinkedHashMap<String, String[]>());		
		chemicalClassDefinitions.get(ChemicalClassName.ALPHA_HYDROXYLATED_FATTY_ACID).put("smarts", new String[]{

		});
		chemicalClassDefinitions.get(ChemicalClassName.ALPHA_HYDROXYLATED_FATTY_ACID).put("negativeSmarts", new String[]{

		});
		chemicalClassDefinitions.put(ChemicalClassName.BETA_HYDROXYLATED_FATTY_ACID, new LinkedHashMap<String, String[]>());		
		chemicalClassDefinitions.get(ChemicalClassName.BETA_HYDROXYLATED_FATTY_ACID).put("smarts", new String[]{

		});
		chemicalClassDefinitions.get(ChemicalClassName.BETA_HYDROXYLATED_FATTY_ACID).put("negativeSmarts", new String[]{

		});
		chemicalClassDefinitions.put(ChemicalClassName.BILE_ACID, new LinkedHashMap<String, String[]>());		
		chemicalClassDefinitions.get(ChemicalClassName.BILE_ACID).put("smarts", new String[]{

		});
		chemicalClassDefinitions.get(ChemicalClassName.BILE_ACID).put("negativeSmarts", new String[]{

		});
		chemicalClassDefinitions.put(ChemicalClassName.ETHER_LIPID, new LinkedHashMap<String, String[]>());		
		chemicalClassDefinitions.get(ChemicalClassName.ETHER_LIPID).put("smarts", new String[]{

		});
		chemicalClassDefinitions.get(ChemicalClassName.ETHER_LIPID).put("negativeSmarts", new String[]{

		});
		chemicalClassDefinitions.put(ChemicalClassName.FATTY_ACID, new LinkedHashMap<String, String[]>());		
		chemicalClassDefinitions.get(ChemicalClassName.FATTY_ACID).put("smarts", new String[]{
				"[#6]!@-,!@=[#6]!@-,!@=[#6]-[#6;X3]([#8;A;X1-,X2H1])=[O;X1]",
				"C=C"

		});
		chemicalClassDefinitions.get(ChemicalClassName.FATTY_ACID).put("negativeSmarts", new String[]{

		});
		
		chemicalClassDefinitions.put(ChemicalClassName.SATURATED_FATTY_ACID, new LinkedHashMap<String, String[]>());		
		chemicalClassDefinitions.get(ChemicalClassName.SATURATED_FATTY_ACID).put("smarts", new String[]{
				"[#6;A;X4;H2,H3][#6;A;H2X4][#6;A;H2X4][#6;X3]([#8;A;X2H1,X1-])=[O;X1]"
		});
		chemicalClassDefinitions.get(ChemicalClassName.SATURATED_FATTY_ACID).put("negativeSmarts", new String[]{
				"C=C",
				"C#C"
		});

		chemicalClassDefinitions.put(ChemicalClassName.UNFUNCTIONALIZED_SATURATED_FATTY_ACID, new LinkedHashMap<String, String[]>());		
		chemicalClassDefinitions.get(ChemicalClassName.UNFUNCTIONALIZED_SATURATED_FATTY_ACID).put("smarts", new String[]{
				"[#6;A;X4;H2,H3][#6;A;H2X4][#6;A;H2X4][#6;X3]([#8;A;X2H1,X1-])=[O;X1]"
		});
		chemicalClassDefinitions.get(ChemicalClassName.UNFUNCTIONALIZED_SATURATED_FATTY_ACID).put("negativeSmarts", new String[]{
				"C=C",
				"C#C",
				"[#6;X4][#8;A;H1X2]"
		});
		
		chemicalClassDefinitions.put(ChemicalClassName.GLYCEROLIPID, new LinkedHashMap<String, String[]>());		
		chemicalClassDefinitions.get(ChemicalClassName.GLYCEROLIPID).put("smarts", new String[]{
				"[$([OX1-,OX2H1,$([OX2](-[CX4])[CX4;R0;H2]),$([OX2]-[CX3]=[CX3]),$([OX2]-[CX3](=O)-[#1,#6])][#6;A;H2X4R0][#6;A;H1X4R0]([OX1-,OX2H1,$([OX2](-[CX4])[CX4;R0;H2]),$([OX2]-[CX3]=[CX3]),$([OX2]-[CX3](=O)-[#1,#6])])[#6;A;H2X4R0][#8]-[CX4,$([CX3]=[CX3]),$([CX3](=O)-[#1,#6])]),$([OX1-,OX2H1,$([OX2](-[CX4])[CX4;R0;H2]),$([OX2]-[CX3]=[CX3]),$([OX2]-[CX3](=O)-[#1,#6])][#6;A;H2X4R0][#6;A;H1X4R0]([#6;A;H2X4R0][OX1-,OX2H1,$([OX2](-[CX4])[CX4;R0;H2]),$([OX2]-[CX3]=[CX3]),$([OX2]-[CX3](=O)-[#1,#6])])[#8;X2]-[CX4,$([CX3]=[CX3]),$([CX3](=O)-[#1,#6])])]"	
		});
		chemicalClassDefinitions.get(ChemicalClassName.GLYCEROLIPID).put("negativeSmarts", new String[]{
		});
		
		chemicalClassDefinitions.put(ChemicalClassName.GLYCEROPHOSPHOLIPID, new LinkedHashMap<String, String[]>());		
		chemicalClassDefinitions.get(ChemicalClassName.GLYCEROPHOSPHOLIPID).put("smarts", new String[]{
				"["
				+ "$([#8]P([!#1!#6;OX1-,OX2H1,$([O]-[#6])])(=O)[#8;R0][#6;A;H2X4R0][#6;A;H1X4R0]([R0;OX1-,OX2H1,$([OX2]-[#1,CX4,$(C(=O)-[#1,#6])])])[#6;A;H2X4R0][R0;OX1-,OX2H1,$([OX2]-[#1,CX4,$(C(=O)-[#1,#6])])]),"
				+ "$([#8]P([!#1!#6;OX1-,OX2H1,$([O]-[#6])])(=O)[#8;R0][#6;A;H2X4R0][#6;A;X3R0](=O)[#6;A;H2X4R0][R0;OX1-,OX2H1,$([OX2]-[#1,CX4,$(C(=O)-[#1,#6])])])]"
		});
		chemicalClassDefinitions.get(ChemicalClassName.GLYCEROPHOSPHOLIPID).put("negativeSmarts", new String[]{

		});
		
		chemicalClassDefinitions.put(ChemicalClassName.SPHINGOLIPID, new LinkedHashMap<String, String[]>());		
		chemicalClassDefinitions.get(ChemicalClassName.SPHINGOLIPID).put("smarts", new String[]{
				"["
//				+ "$([#6;A;H3X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H1X4]([#8;A;X2H1,X1-])[#6;A;H1X4]([!#1!#6;NX4H3+,NX3H2,$([#7X3]C(=O)C)])[#6;A;H2X4][#8]),"
//				+ "$([#6;A;H3X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;H1X3]=[#6;H1X3][#6;A;H1X4]([#8;A;X2H1,X1-])[#6;A;H1X4]([!#1!#6;NX4H3+,NX3H2,$([#7X3]C(=O)C)])[#6;A;H2X4][#8])"

				+ "$([#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H1X4]([#8;A;X2H1,X1-])[#6;A;H1X4]([!#1!#6;NX4H3+,NX3H2,$([#7X3]C(=O)C)])[#6;A;H2X4][#8]),"
//				+ "$([#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;H1X3]=[#6;H1X3][#6;A;H1X4]([#8;A;X2H1,X1-])[#6;A;H1X4]([!#1!#6;NX4H3+,NX3H2,$([#7X3]C(=O)C)])[#6;A;H2X4][#8])"
				+ "$([H][#6]([#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4])-,=[#6]([H])[#6;A;H2X4]C([H])([#1,OX2H1])[#6;H1X3]=[#6;H1X3][#6;A;H1X4]([#8;A;X2H1,X1-])[#6;A;H1X4]([!#1!#6;NX4H3+,NX3H2,$([#7X3]C(=O)C)])[#6;A;H2X4][#8]),"
				// 	5-hydroxy,3E-sphingosine (LMSP01080004)
				+ "$([H][#8]C([H])([#6;A;H2X4][#6]([H])-,=[#6]([H])[#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4])[#6;H1X3]=[#6;H1X3][#6;A;H1X4]([!#1!#6;NX4H3+,NX3H2,$([#7X3]C(=O)C)])[#6;A;H2X4][#8])"
				+ "]"
		});
		chemicalClassDefinitions.get(ChemicalClassName.SPHINGOLIPID).put("negativeSmarts", new String[]{
		});		
		
		chemicalClassDefinitions.put(ChemicalClassName.ETHER_LIPID, new LinkedHashMap<String, String[]>());		
		chemicalClassDefinitions.get(ChemicalClassName.ETHER_LIPID).put("smarts", new String[]{
				"[$([#8;X2][#6;A;H2X4]!@-[#6;A;X4](!@-[!#1!#6;$([OX2H1]),$([OX2]-[CX4H2]-[#6;A])])!@-[#6;A;H2X4][!#1!#6;$([OX2H1]),$([OX2]-[CX4H2]-[#6;A])]),"
				+ "$([#8]!@-[#6;A;X4](!@-[#6;A;H2X4][!#1!#6;$([OX2H1]),$([OX2]-[CX4H2]-[#6;A])])!@-[#6;A;H2X4][!#1!#6;$([OX2H1]),$([OX2]-[CX4H2]-[#6;A])]),"
				+ "$([!#1!#6;$([OX2H1]),$([OX2]-[CX4H2]-[#6;A])][#6;A;H2X4]!@-[#6;A;X3](!@=[O;X1])!@-[#6;A;H2X4][!#1!#6;$([OX2H1]),$([OX2]-[CX4H2]-[#6;A])]),"
				+ "$([#8;X2][#6;A;H2X4]!@-[#6;A;X3](!@=[O;X1])!@-[#6;A;H2X4][!#1!#6;$([OX2H1]),$([OX2]-[CX4H2]-[#6;A])])"
				+ "]"
		});
		chemicalClassDefinitions.get(ChemicalClassName.ETHER_LIPID).put("negativeSmarts", new String[]{
		});
		
			
		chemicalClassDefinitions.put(ChemicalClassName.GLYCEROL_3_PHOSPHATE_INOSITOL, new LinkedHashMap<String, String[]>());		
		chemicalClassDefinitions.get(ChemicalClassName.GLYCEROL_3_PHOSPHATE_INOSITOL).put("smarts", new String[]{
				"[#8][#6;A;H1X4R1]1[#6;A;H1X4R1]([#8])[#6;A;H1X4R1]([#8])[#6;A;H1X4R1]([#8]P([#8;X2H1,X1-])(=O)[#8]-[#6;H2X4]-[#6;H1X4](-[#6;H2X4]-[#8]-[#6]([#6,#1;A])=O)-[#8]-[#6]([#6,#1;A])=O)[#6;A;H1X4R1]([#8])[#6;A;H1X4R1]1[#8]"
		});
		chemicalClassDefinitions.get(ChemicalClassName.GLYCEROL_3_PHOSPHATE_INOSITOL).put("negativeSmarts", new String[]{
		});
		
		chemicalClassDefinitions.put(ChemicalClassName.C24_BILE_ACID, new LinkedHashMap<String, String[]>());		
		chemicalClassDefinitions.get(ChemicalClassName.C24_BILE_ACID).put("smarts", new String[]{
				"["
				+ "$([#6;A;H3X4][#6;A;H1X4]([#6;A;H2X4][#6;A;H2X4][#6]([!#1!#6;OX1,OX2H1,NX3H1])=O)[#6]1-[#6]-,=[#6]-[#6]2-[#6]-,=3-[#6;CX4H2,$([CX4H1]-[OX2H1])]-[#6]-,=[#6]4-,=[#6]-,=[#6]([#8;A;X2H1,X1-])-[#6]-[#6]C4([#6;A;H3X4])[#6]-,=3-[#6]-[#6;CX4H2,$([CX4H1]-[OX2H1])]C12[#6;A;H3X4])"
				+ "]"
		});
		chemicalClassDefinitions.get(ChemicalClassName.C24_BILE_ACID).put("negativeSmarts", new String[]{
		});		
		
		chemicalClassDefinitions.put(ChemicalClassName.C23_BILE_ACID, new LinkedHashMap<String, String[]>());		
		chemicalClassDefinitions.get(ChemicalClassName.C23_BILE_ACID).put("smarts", new String[]{
				"["
				+ "$([#6;A;H3X4][#6;A;H1X4]([#6;A;H2X4][#6]([!#1!#6;OX1,OX2H1,NX3H1])=O)[#6]1-[#6]-,=[#6]-[#6]2-[#6]-,=3-[#6;CX4H2,$([CX4H1]-[OX2H1])]-[#6]-,=[#6]4-,=[#6]-,=[#6]([#8;A;X2H1,X1-])-[#6]-[#6]C4([#6;A;H3X4])[#6]-,=3-[#6]-[#6;CX4H2,$([CX4H1]-[OX2H1])]C12[#6;A;H3X4])"
				+ "]"
		});
		chemicalClassDefinitions.get(ChemicalClassName.C23_BILE_ACID).put("negativeSmarts", new String[]{
		});		

		chemicalClassDefinitions.put(ChemicalClassName.SULFATE_ESTER, new LinkedHashMap<String, String[]>());		
		chemicalClassDefinitions.get(ChemicalClassName.SULFATE_ESTER).put("smarts", new String[]{
				"["
				+ "$([#6]-[#8;X2][S;X4]([#8;A;X1-,X2H1])(=[O;X1])=[O;X1])"
				+ "]"
		});
		chemicalClassDefinitions.get(ChemicalClassName.SULFATE_ESTER).put("negativeSmarts", new String[]{
		});
		
//		chemicalClassDefinitions.put(ChemicalClassName.GLYCEROPHOSPHOLIPID, "");
//		chemicalClassDefinitions.put(ChemicalClassName.OMEGA_HYDROXYLATED_FATTY_ACID, "");
//		chemicalClassDefinitions.put(ChemicalClassName.SPHINGOLIPID, "");
//		chemicalClassDefinitions.put(ChemicalClassName.UNSATURATED_FATTY_ACID, "");
	
	}
	
	
	public static ChemicalClassName findChemicalClass(IAtomContainer molecule) throws SMARTSException, CloneNotSupportedException, CDKException{
		ChemicalClassName chemClass = null;
//		IAtomContainer molecule = mol.clone();
//		AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(molecule);
//		AtomContainerManipulator.convertImplicitToExplicitHydrogens(molecule);
		
		if(isAlphaHydroxyFattyAcid(molecule)){
			chemClass = ChemicalClassName.ALPHA_HYDROXYLATED_FATTY_ACID;
		} else if(isBetaHydroxyFattyAcid(molecule)){
			chemClass = ChemicalClassName.BETA_HYDROXYLATED_FATTY_ACID;
		} else if(isOmegaHydroxyFattyAcid(molecule)){
			chemClass = ChemicalClassName.OMEGA_HYDROXYLATED_FATTY_ACID;
		} else if(isUnsaturatedFattyAcid(molecule)){
			chemClass = ChemicalClassName.UNSATURATED_FATTY_ACID;
		} else if(isFattyAcid(molecule)){
			chemClass = ChemicalClassName.FATTY_ACID;
		}  		
		else {
			for(Map.Entry<ChemicalClassName, LinkedHashMap<String, String[]>> cc : chemicalClassDefinitions.entrySet()){
				if( 
						cc.getKey() == ChemicalClassName.GLYCEROL_3_PHOSPHATE_INOSITOL || 
						cc.getKey() == ChemicalClassName.GLYCEROLIPID || 
						cc.getKey() == ChemicalClassName.GLYCEROPHOSPHOLIPID ||
						cc.getKey() == ChemicalClassName.ETHER_LIPID || 
						cc.getKey() == ChemicalClassName.SPHINGOLIPID ||
						cc.getKey() == ChemicalClassName.C24_BILE_ACID ||
						cc.getKey() == ChemicalClassName.C23_BILE_ACID){
					
					boolean b = true;
					
					for(String smart : chemicalClassDefinitions.get(cc.getKey()).get("smarts")){
						SmartsPatternCDK pattern = new SmartsPatternCDK(smart);	
						if(!(pattern.hasSMARTSPattern(molecule)>0)){
							b = false;
							break;
						}
					}					
					for(String negativeSmart : chemicalClassDefinitions.get(cc.getKey()).get("negativeSmarts")){
						SmartsPatternCDK npattern = new SmartsPatternCDK(negativeSmart);	
						if(npattern.hasSMARTSPattern(molecule)>0){
							b = false;
							break;
						}
					}
					
					if(b == true){
						chemClass = cc.getKey();
					}					
				}
			}
		}
		

		
		
		
		return chemClass;
	}
	
	
	public static boolean isUnsubstitutedSatudatedFattyAcid(IAtomContainer molecule){	
		
		return false;
	}
	
	
	public static boolean isFattyAcid(IAtomContainer molecule){	
		
		return false;
	}
	
	
	public static boolean isUnsaturatedFattyAcid(IAtomContainer molecule){	
		
		return false;
	}
	
	public static boolean isSaturatedFattyAcid(IAtomContainer molecule) throws SMARTSException{	
		boolean sfa = true;

//		for(String smart : chemicalClassDefinitions.get(ChemicalClassName.SATURATED_FATTY_ACID).get("smarts")){
//			SmartsPatternCDK pattern = new SmartsPatternCDK(smart);	
//			if(!(pattern.hasSMARTSPattern(molecule)>0)){
//				sfa = false;
//				break;
//			}
//		}	
//		
//		for(String negativeSmart : chemicalClassDefinitions.get(ChemicalClassName.ETHER_LIPID).get("negativeSmarts")){
//			SmartsPatternCDK npattern = new SmartsPatternCDK(negativeSmart);	
//			if(npattern.hasSMARTSPattern(molecule)>0){
//				sfa = false;
//				break;
//			}
//		}
//		
//		// As per observation, 
		
		
		
		return sfa;
	}
	
	public static boolean isAlphaHydroxyFattyAcid(IAtomContainer molecule){	
		
		return false;
	}
	public static boolean isBetaHydroxyFattyAcid(IAtomContainer molecule){	
		
		return false;
	}
	public static boolean isOmegaHydroxyFattyAcid(IAtomContainer molecule){	
		
		return false;
	}
	public static boolean isEtherLipid(IAtomContainer molecule) throws SMARTSException{	
		
		boolean b = true;
		
		for(String smart : chemicalClassDefinitions.get(ChemicalClassName.ETHER_LIPID).get("smarts")){
			SmartsPatternCDK pattern = new SmartsPatternCDK(smart);	
			if(!(pattern.hasSMARTSPattern(molecule)>0)){
				b = false;
				break;
			}
		}					
		for(String negativeSmart : chemicalClassDefinitions.get(ChemicalClassName.ETHER_LIPID).get("negativeSmarts")){
			SmartsPatternCDK npattern = new SmartsPatternCDK(negativeSmart);	
			if(npattern.hasSMARTSPattern(molecule)>0){
				b = false;
				break;
			}
		}	
		return b;
	}	

	public static boolean isGlyceroLipid(IAtomContainer molecule) throws SMARTSException{	
		
		boolean b = true;
		
		for(String smart : chemicalClassDefinitions.get(ChemicalClassName.GLYCEROLIPID).get("smarts")){
			SmartsPatternCDK pattern = new SmartsPatternCDK(smart);	
			if(!(pattern.hasSMARTSPattern(molecule)>0)){
				b = false;
				break;
			}
		}					
		for(String negativeSmart : chemicalClassDefinitions.get(ChemicalClassName.GLYCEROLIPID).get("negativeSmarts")){
			SmartsPatternCDK npattern = new SmartsPatternCDK(negativeSmart);	
			if(npattern.hasSMARTSPattern(molecule)>0){
				b = false;
				break;
			}
		}
		
		return b;

	}
	
	public static boolean isGlycerophosphoLipid(IAtomContainer molecule) throws SMARTSException{	
		
		boolean b = true;
		
		for(String smart : chemicalClassDefinitions.get(ChemicalClassName.GLYCEROPHOSPHOLIPID).get("smarts")){
			SmartsPatternCDK pattern = new SmartsPatternCDK(smart);	
			if(!(pattern.hasSMARTSPattern(molecule)>0)){
				b = false;
				break;
			}
		}					
		for(String negativeSmart : chemicalClassDefinitions.get(ChemicalClassName.GLYCEROPHOSPHOLIPID).get("negativeSmarts")){
			SmartsPatternCDK npattern = new SmartsPatternCDK(negativeSmart);	
			if(npattern.hasSMARTSPattern(molecule)>0){
				b = false;
				break;
			}
		}
		
		return b;

	}	
	
	public static boolean isSphingoLipid(IAtomContainer molecule) throws SMARTSException{	
		
		boolean b = true;
		
		for(String smart : chemicalClassDefinitions.get(ChemicalClassName.SPHINGOLIPID).get("smarts")){
			SmartsPatternCDK pattern = new SmartsPatternCDK(smart);	
			if(!(pattern.hasSMARTSPattern(molecule)>0)){
				b = false;
				break;
			}
		}					
		for(String negativeSmart : chemicalClassDefinitions.get(ChemicalClassName.SPHINGOLIPID).get("negativeSmarts")){
			SmartsPatternCDK npattern = new SmartsPatternCDK(negativeSmart);	
			if(npattern.hasSMARTSPattern(molecule)>0){
				b = false;
				break;
			}
		}
		
		return b;

	}	

	public static boolean isGlycerol_3_PhosphateInositol(IAtomContainer molecule) throws SMARTSException{	
		
		boolean b = true;
		
		for(String smart : chemicalClassDefinitions.get(ChemicalClassName.GLYCEROL_3_PHOSPHATE_INOSITOL).get("smarts")){
			SmartsPatternCDK pattern = new SmartsPatternCDK(smart);	
			if(!(pattern.hasSMARTSPattern(molecule)>0)){
				b = false;
				break;
			}
		}					
		for(String negativeSmart : chemicalClassDefinitions.get(ChemicalClassName.GLYCEROL_3_PHOSPHATE_INOSITOL).get("negativeSmarts")){
			SmartsPatternCDK npattern = new SmartsPatternCDK(negativeSmart);	
			if(npattern.hasSMARTSPattern(molecule)>0){
				b = false;
				break;
			}
		}
		
		return b;

	}	

	public static boolean isC23BileAcid(IAtomContainer molecule) throws SMARTSException{	
		
		boolean b = true;
		
		for(String smart : chemicalClassDefinitions.get(ChemicalClassName.C23_BILE_ACID).get("smarts")){
			SmartsPatternCDK pattern = new SmartsPatternCDK(smart);	
			if(!(pattern.hasSMARTSPattern(molecule)>0)){
				b = false;
				break;
			}
		}					
		for(String negativeSmart : chemicalClassDefinitions.get(ChemicalClassName.C23_BILE_ACID).get("negativeSmarts")){
			SmartsPatternCDK npattern = new SmartsPatternCDK(negativeSmart);	
			if(npattern.hasSMARTSPattern(molecule)>0){
				b = false;
				break;
			}
		}
		
		return b;

	}
	
	public static boolean isC24BileAcid(IAtomContainer molecule) throws SMARTSException{	
		
		boolean b = true;
		
		for(String smart : chemicalClassDefinitions.get(ChemicalClassName.C24_BILE_ACID).get("smarts")){
			SmartsPatternCDK pattern = new SmartsPatternCDK(smart);	
			if(!(pattern.hasSMARTSPattern(molecule)>0)){
				b = false;
				break;
			}
		}					
		for(String negativeSmart : chemicalClassDefinitions.get(ChemicalClassName.C24_BILE_ACID).get("negativeSmarts")){
			SmartsPatternCDK npattern = new SmartsPatternCDK(negativeSmart);	
			if(npattern.hasSMARTSPattern(molecule)>0){
				b = false;
				break;
			}
		}
		
		return b;

	}	
	
	
	public static boolean isSulfateEster(IAtomContainer molecule) throws SMARTSException{	
		
		boolean b = true;
		
		for(String smart : chemicalClassDefinitions.get(ChemicalClassName.SULFATE_ESTER).get("smarts")){
			SmartsPatternCDK pattern = new SmartsPatternCDK(smart);	
			if(!(pattern.hasSMARTSPattern(molecule)>0)){
				b = false;
				break;
			}
		}					
		for(String negativeSmart : chemicalClassDefinitions.get(ChemicalClassName.SULFATE_ESTER).get("negativeSmarts")){
			SmartsPatternCDK npattern = new SmartsPatternCDK(negativeSmart);	
			if(npattern.hasSMARTSPattern(molecule)>0){
				b = false;
				break;
			}
		}	
		return b;
	}
	
	
	
	
	public static void main(String[] args) throws SMARTSException, CloneNotSupportedException, CDKException{
		ChemicalClassFinder ccf = new ChemicalClassFinder();
		IChemObjectBuilder 	builder = SilentChemObjectBuilder.getInstance();
		SmilesParser	smiParser		= new SmilesParser(builder);
		IAtomContainer molecule = smiParser.parseSmiles("CCCCCCCCCCCCCCCC(=O)OCC(COP(=O)(O)OC1C(C(C(C(C1O)O)O)O)O)OC(=O)CCCCCCCCCCCCCCC");
		System.out.println(ccf.findChemicalClass(molecule));
	}
	
}

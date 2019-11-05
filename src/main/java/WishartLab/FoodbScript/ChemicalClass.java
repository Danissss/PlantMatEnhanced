package WishartLab.FoodbScript;

import org.openscience.cdk.interfaces.IAtomContainer;

import ambit2.smarts.query.SMARTSException;
import ambit2.smarts.query.SmartsPatternCDK;

public class ChemicalClass {
	
	/**
	 * check if the compound is in polyphenol class
	 */
	public static boolean is_polyphenol(IAtomContainer mol) {
		
		String polyphenol_smarts_2 = "[#8;A;H1X2][#6;R1]=,:1[#6;R1]=,:[#6;R1][#6;R1]=,:[#6;R1][#6;R1]=,:1";
		String polyphenol_smarts_3 = "[$([#8;A;H1X2][#6;R1]-1=[#6;R1]-[#6;R1]=[#6;R1]-[#6;R1]=[#6;R1]-1),$([#8;A;H1X2][c;R1]1[c;R1][c;R1][c;R1][c;R1][c;R1]1)]";
		
		// SmartsPatternCDK smarts = new SmartsPatternCDK();
		SmartsPatternCDK smarts2 = new SmartsPatternCDK();
		SmartsPatternCDK smarts3 = new SmartsPatternCDK();
		try {
			// smarts.setSmarts(polyphenol_smarts);
			smarts2.setSmarts(polyphenol_smarts_2);
			smarts3.setSmarts(polyphenol_smarts_3);
			// smarts.hasSMARTSPattern(mole) return number of encountered smarts
			if (smarts2.hasSMARTSPattern(mol) != 0 || smarts3.hasSMARTSPattern(mol) != 0) {
				return true;
			} else {
				return false;
			}
		} catch (SMARTSException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			return false;
		}

	}

	/**
	 * check if the compound is in phenol class
	 */
	public static boolean is_phenol(IAtomContainer mol) {
		String phenol_smarts = "[$([#8;A;H1X2][#6;R1]-1=[#6;R1]-[#6;R1]=[#6;R1]-[#6;R1]=[#6;R1]-1),$([#8;A;H1X2][c;R1]1[c;R1][c;R1][c;R1][c;R1][c;R1]1)]";
		SmartsPatternCDK smarts = new SmartsPatternCDK();
		try {
			// smarts.setSmarts(polyphenol_smarts);
			smarts.setSmarts(phenol_smarts);
			// smarts.hasSMARTSPattern(mole) return number of encountered smarts
			if (smarts.hasSMARTSPattern(mol) != 0) {
				return true;
			} else {
				return false;
			}
		} catch (SMARTSException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			return false;
		}
		
	}
	
	
	/**
	 * check if the compound is in terpenoids class
	 */
	public static boolean is_terpenoid(IAtomContainer mol) {
		
		String terpenoids_smarts = "[#6]-[#6](=[#6])-[#6]=[#6]";
		
		SmartsPatternCDK smarts = new SmartsPatternCDK();
		try {
			smarts.setSmarts(terpenoids_smarts);
			// smarts.hasSMARTSPattern(mole) return number of encountered smarts
			if (smarts.hasSMARTSPattern(mol) != 0) {
				return true;
			} else {
				return false;
			}
		} catch (SMARTSException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			return false;
		}

	}
}

package WishartLab.PlantGlycosider;

import java.util.HashMap;

public class Sugars {
	
	
	/**
	 * get int[] for all sugars
	 * @return
	 */
	public Integer[] GetSugarIndexSet() {
		Integer[] output = new Integer[SugarIndex().size()];
		for (int i = 1 ; i <= SugarIndex().size(); i++) {
			output[i-1] = i;
		}
		return output;
	}
	
	
	/**
	 * map for sugar structure
	 * should add parse file ability, then no need to change code every time
	 * should strictly match HashMap<Integer, String> SugarIndex()
	 * @return
	 */
	public HashMap<Integer, String> SugarGroup(){
		HashMap<Integer, String> output = new HashMap<Integer, String>();
		
		// index are based on ChemAxon
		// pentoses -> Aldopentoses
		output.put(1, "OC1COC(O)C(O)C1O");
		// pentoses -> Ketopentoses
		output.put(2, "OCC(O)C(O)C(=O)CO");
		// pentoses -> Deoxy Sugars
		output.put(3, "OCC(O)C(O)CC=O");
		// hexoses -> Aldohexoses
		output.put(4, "OCC1OC(O)C(O)C(O)C1O");
		// hexoses -> Ketohexoses
		output.put(5, "OCC1(O)OCC(O)C(O)C1O");
		
//		// Heptose -> aldoheptose 
//		output.put(6,"OCC(O)C1OC(O)C(O)C(O)C1O");
//		// Heptose -> Ketoheptose (Sedoheptulose)
//		output.put(7, "OCC(O)C(O)C(O)C(O)C(=O)CO");
		
		
		
		return output;
	}
	
	
	/**
	 * map for index
	 * should add parse file ability, then no need to change code every time
	 * should strictly match HashMap<Integer, String> SugarGroup()
	 * @return
	 */
	public HashMap<Integer,Integer> SugarIndex(){
		HashMap<Integer,Integer> output = new HashMap<Integer,Integer>();
		output.put(1,8);
		output.put(2,null);
		output.put(3,null);
		output.put(4,null);
		output.put(5,null);
		
		
		return output;
	}
}

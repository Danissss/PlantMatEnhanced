package WishartLab.PlantGlycosider;

import java.util.HashMap;

public class Sugars {
	
	
	/**
	 * get int[] for all sugars
	 * @return
	 */
	public Integer[] GetSugarIndexSet(HashMap<Integer,Integer> sugar_index) {
		
		Integer[] output = new Integer[sugar_index.size()];
		int i = 0;
		for(Integer index : sugar_index.keySet()) {
			output[i] = index;
			i++;
		}
		
		return output;
	}
	
	
	/**
	 * map for sugar structure
	 * should add parse file ability, then no need to change code every time
	 * should strictly match HashMap<Integer, String> SugarIndex()
	 * Either is 5 membrane or 6 membrane?
	 * Fructose exists to the extent of about 80% in the pyranose form and about 20% as the five-membered furanose form resulting from addition of the -OH group at C5 to the C2 carbonyl group.
	 * http://www.chem.latech.edu/~deddy/chem121/Carbohydrates.htm
	 * 
	 * @return
	 */
	public HashMap<Integer, String> SugarGroup(){
		HashMap<Integer, String> output = new HashMap<Integer, String>();
		
		// index are based on ChemAxon
		// pentoses -> Aldopentoses
		output.put(1, "OC1COC(O)C(O)C1O");
		// pentoses -> Ketopentoses
		output.put(2, "OCC1OC(O)C(O)C1O");
		// pentoses -> Deoxy Sugars // N.V. BHAGAVAN, in Medical Biochemistry (Fourth Edition), 2002 (Deoxy Sugars)
		output.put(3, "OCC1OC(O)CC1O");
//		// hexoses -> 
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
	 * Most ring structure sugar will be attached by the -OH near the -o-;
	 * Most attachment should happen around -o-
	 * @return
	 */
	public HashMap<Integer,Integer> SugarAttachMapIndex(){
		HashMap<Integer,Integer> output = new HashMap<Integer,Integer>();
		output.put(1,5);
		output.put(2,5);
		output.put(3,5);
		output.put(4,5);
		output.put(5,0);
		
		
		return output;
	}
}

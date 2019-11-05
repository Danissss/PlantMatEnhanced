package WishartLab.FoodbScript;

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
//		output.put(1,"C1[C@H]([C@H]([C@@H](C(O1)(CO)O)O)O)O"); // d_fructose
//		output.put(2,"C([C@@H]1[C@H]([C@@H]([C@H](C(O1)O)O)O)O)O"); // d_glucose
//		output.put(3,"C([C@@H]1[C@@H]([C@@H]([C@H](C(O1)O)O)O)O)O"); // d_galactose
//		output.put(4,"C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O[C@]2([C@H]([C@@H]([C@H](O2)CO)O)O)CO)O)O)O)O"); // sucrose
//		output.put(5,"C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O[C@@H]2[C@H](O[C@@H]([C@@H]([C@H]2O)O)O)CO)O)O)O)O"); // maltose
//		output.put(6,"C([C@@H]1[C@@H]([C@@H]([C@H]([C@@H](O1)O[C@@H]2[C@H](OC([C@@H]([C@H]2O)O)O)CO)O)O)O)O"); // lactose
//		output.put(7, "C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O[C@@H]2[C@@H]([C@H]([C@@H]([C@H](O2)CO)O)O)O)O)O)O)O"); // Trehalose
//		output.put(8, "C1[C@H]([C@H]([C@H](C(O1)O)O)O)O"); //d_ribose
//		output.put(9, "C1[C@H]([C@@H]([C@H](C(O1)O)O)O)O"); //d_xylose
//		output.put(10, "C1[C@@H]([C@@H]([C@H](C(O1)O)O)O)O"); // l_Arabinose
//		output.put(11, "C([C@H]([C@H](C(=O)CO)O)O)O"); // d_ribulose
//		output.put(12, "C([C@@H]1[C@H]([C@@H]([C@@H](C(O1)O)O)O)O)O"); // d_mannose
//		output.put(13, "C([C@@H]1[C@@H]([C@@H]([C@H]([C@H](O1)OC[C@@H]2[C@H]([C@@H]([C@H](C(O2)O)O)O)O)O)O)O)O"); //d_melibiose
		
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
//		output.put(6,null);
//		output.put(7,null);
//		output.put(8, 0);
//		output.put(9, 0);
//		output.put(10, 0);
//		output.put(11, 9);
//		output.put(12, 11);
//		output.put(13, 22);
		
		
		return output;
	}
}

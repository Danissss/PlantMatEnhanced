/** 
 * In such case, MS is insufficient to provide confident positional identification, 
 * and thus, positional isomers are not included in PlantMAT identifications.
 * PlantMAT only considers the instances where the aglycone is attached with the maximum of two sugar chains 
 * in linear fashion
 * create permutation of sugar groups
 * add right charges 
 * add charges to OH group in original group
 * for every attachment, remove the charge
 */

package WishartLab.PlantGlycosider;

import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Set;
import java.util.Stack;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.TimeoutException;

import org.apache.commons.lang3.ArrayUtils;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import com.opencsv.CSVParser;
import com.opencsv.CSVParserBuilder;
import com.opencsv.CSVReader;
import com.opencsv.CSVReaderBuilder;
import com.opencsv.CSVWriter;





// run as java -jar PlantMat.jar <core_number> <file_name> <running time>

public class PlantGlycosider {
	
	private PreValidation preval = new PreValidation();
	private PostValidation postval = new PostValidation();
	
	/**
	 * return empty ArrayList<String> if something goes wrong...
	 * @param mole
	 * @param sugar_structure
	 * @param sugerindex
	 * @param sugar_list
	 * @param s_list
	 * @return
	 */
	public IAtomContainerSet permutatedCompoundInSmiles(IAtomContainer mole, HashMap<Integer,String> sugar_structure, 
			HashMap<Integer,Integer> sugerindex, Integer[] sugar_list, Integer[] s_list){
		
		IChemObjectBuilder builder = DefaultChemObjectBuilder.getInstance();
		IAtomContainerSet result_mole = builder.newInstance(IAtomContainerSet.class);
		CheminformaticUtility permuated_chem_util = new CheminformaticUtility();
		
		ArrayList<Integer[]> temp_sugar_list = permuated_chem_util.GetAllPossibleSet(sugar_list,s_list.length); // temp set will be list of all possible of order of sugars
		int[] oh_list = new int[s_list.length];
		for(int i=0;i<s_list.length;i++) {oh_list[i]=s_list[i];}
		for(Integer[] sugar_comb : temp_sugar_list) {
			// sugar_comb contain the index of selected sugar
			String[] smiles = new String[sugar_comb.length];
			int[] smiles_index = new int[sugar_comb.length];
			
			for (int i=0;i<smiles.length;i++) {smiles[i] = sugar_structure.get(sugar_comb[i]);}
			for (int i=0;i<smiles_index.length;i++) {smiles_index[i] = sugerindex.get(sugar_comb[i]);}
			
			try {
				
				IAtomContainer copied_mole = mole.clone();
				IAtomContainer mole_man = AddSugarGroupPermu(copied_mole,oh_list,smiles,smiles_index);
				result_mole.addAtomContainer(CheminformaticUtility.SplitContainer(mole_man));
				
			} catch (CloneNotSupportedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			} catch (InvalidSmilesException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			} catch (CDKException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			
		}
		
		return result_mole;
		
		
		
	}
	/**
	 * Sugar's stereochemistry won't lost
	 * if OH number is larger than number of sugar, then permutated all possible add-on for missing sugars
	 * will be called by AddSugarGroupPermu; AddSugarGroupPermu will be called by PreformTransformation
	 * 1. get list of OH of given compound
	 * 2. do permutation of OH with combination of sugars 
	 * @param mole
	 * @throws Exception
	 */
	public IAtomContainerSet PreformTransformation(IAtomContainer mole, int restrictOH) throws Exception {
		
		
		IChemObjectBuilder builder = DefaultChemObjectBuilder.getInstance();
		IAtomContainerSet result_mole = builder.newInstance(IAtomContainerSet.class);
		CheminformaticUtility new_uti = new CheminformaticUtility();
		
		
		Set<Integer> OH_list = CheminformaticUtility.GetOHList(mole); // return list of OH
		
		if (OH_list.size() > restrictOH) {
			return result_mole;
		}
		
		
		
		ArrayList<Integer[]> combined_list = new ArrayList<Integer[]>();
		for(int i =  1; i <= OH_list.size(); i++) {
			new_uti.permutations(OH_list, new Stack<Integer>(), i);
		}
		combined_list.addAll(new_uti.get_combination_list());
		
		
		// for each combination of OH poisition, get all combination of sugar
		Sugars sugar = new Sugars();
		
		HashMap<Integer,String> sugar_structure = sugar.SugarGroup();
		HashMap<Integer,Integer> sugerindex = sugar.SugarAttachMapIndex();
		Integer[] sugar_list = sugar.GetSugarIndexSet(sugerindex); // index of all available sugars
		
		if (sugar_list.length < OH_list.size()) {
			// if number of sugar is less than number of o-glycoside point, then permuate the possible sugar to add to the list
			// if #OH = 5, #sugar = 4; then we need another sugar, after get that sugar, we run permutatedCompoundInSmiles();
			// if #OH = 5, #sugar = 4; then we need two sugars, perumate the existing sugar for two combination;
			int extra_space = OH_list.size() - sugar_list.length;
			CheminformaticUtility extrace_space_util = new CheminformaticUtility();
			ArrayList<Integer[]> extract_combination = extrace_space_util.GetAllPossibleSet(sugar_list,extra_space); // combination of sugars 
			
			for(int i = 0; i < extract_combination.size(); i++) {
				for(int j = 0; j < combined_list.size(); j++) {
					
					Integer[] new_sugar_list = ArrayUtils.addAll(extract_combination.get(i), combined_list.get(j));
					result_mole.add(permutatedCompoundInSmiles(mole,sugar_structure, sugerindex, sugar_list,new_sugar_list));	
				
				}
			}
			
			
		}else{
			// sugar_list.length > OH_list.size() and sugar_list.length = OH_list.size()
			// if number of sugar is larger than position of OHs; 
			// for each combined_list (list of array of OH position); add each sugar one at time	
			for(Integer[] s_list : combined_list) {
				result_mole.add(permutatedCompoundInSmiles(mole,sugar_structure, sugerindex, sugar_list,s_list));		
			}
			
		}
		
		
		
		return result_mole;
		
		
	}
	
	
	/**
	 * Add sugar group with multiple an selective OHs
	 * Should only generate one container
	 * @param mole
	 * @param mole_index
	 * @param sugar_smiles
	 * @param sugar_index
	 * @return
	 * @throws InvalidSmilesException 
	 */
	public IAtomContainer AddSugarGroupPermu(IAtomContainer mole, int[] mole_index, String[] sugar_smiles, int[] sugar_index) throws InvalidSmilesException {
		
		IChemObjectBuilder builder = DefaultChemObjectBuilder.getInstance();
		IAtomContainer output = builder.newAtomContainer();
		
		if (sugar_smiles.length < mole_index.length) {
			// number of sugar smiles can't exceed the number of OH
			// can be equal
			
		}else if (mole_index.length > sugar_smiles.length) {
			// when there are more OH than smiles
			// randomly add repeated sugar
		}
		
		if(sugar_smiles.length == 1) {
			// all OH have to have same stuff
			String[] new_sugar_smiles = new String[mole_index.length];
			int[] new_sugar_index = new int[mole_index.length];
			for(int i=0; i < mole_index.length; i++) {
				new_sugar_smiles[i] = sugar_smiles[0];
				new_sugar_index[i] = sugar_index[0];
				
			}

			
			// add one sugar at time
			// increment charges start from two
			// keep all sugar charges one
			// connect the sugar and compound with order
			
			IAtomContainerSet sugar_mole_set = builder.newInstance(IAtomContainerSet.class);			
			for(int i=0; i < mole_index.length;i++) {
				sugar_mole_set.addAtomContainer(CheminformaticUtility.TransformSugarSmilesToContainer(new_sugar_smiles[i],new_sugar_index[i]));
			}
			
			try {
				output = GlycosiderUtility.CombineContainers(mole,mole_index,sugar_mole_set);
			} catch (InvalidSmilesException e) {
				e.printStackTrace();
			} catch (CDKException e) {
				e.printStackTrace();
			}
			
		}else if(sugar_smiles.length > 1) {
			// following the order
			// add whatever how many OH contains.
			// when sugar length is less than OH index

			// prepare arrays for CombineContainers
			// mole_index is the array of OH; sugar_smiles is array of sugar smiles;
			String[] new_sugar_smiles = new String[mole_index.length];
			int[] new_sugar_index = new int[mole_index.length];
			IAtomContainerSet sugar_mole_set = builder.newInstance(IAtomContainerSet.class);
			
			int last_sugar_index = sugar_index[sugar_index.length - 1];
			String last_smiles = sugar_smiles[sugar_smiles.length - 1];
			
			for(int i=0; i < sugar_smiles.length; i++) {
				new_sugar_smiles[i] = sugar_smiles[i];
				new_sugar_index[i] = sugar_index[i];
				sugar_mole_set.addAtomContainer(CheminformaticUtility.TransformSugarSmilesToContainer(new_sugar_smiles[i],new_sugar_index[i]));
			}
			for(int i=sugar_smiles.length-1; i < mole_index.length; i++) {
				new_sugar_smiles[i] = last_smiles;
				new_sugar_index[i] = last_sugar_index;
				sugar_mole_set.addAtomContainer(CheminformaticUtility.TransformSugarSmilesToContainer(last_smiles,last_sugar_index));
			}
			
			try {
				output = GlycosiderUtility.CombineContainers(mole,mole_index,sugar_mole_set);
			} catch (InvalidSmilesException e) {
				e.printStackTrace();
			} catch (CDKException e) {
				e.printStackTrace();
			}
		}
		
		return output;
	}
	
	
	
	/**
	 * parse and transformer single compounds from foodb
	 * get all transformerable compounds and save them to directory
	 * then get these file add to new food
	 * make sure that all duplications are removed
	 * @param data_file
	 * @throws IOException 
	 * @throws InterruptedException 
	 */
	public int TransformerFoodbSingleCompound(String smiles,String inchikey) throws IOException, InterruptedException {
		
		int total_num = 0;
		PlantGlycosider fdb = new PlantGlycosider();
		try {
			IAtomContainer mole = CheminformaticUtility.parseSmilesToContainer(smiles);
			if (preval.isCompoundValidate(mole)) {
				IAtomContainerSet transformedSmiles = fdb.PreformTransformation(mole,4);
				if (transformedSmiles.getAtomContainerCount() != 0) {
					IAtomContainerSet postValidateCompound = postval.validateCompound(transformedSmiles);
					CheminformaticUtility.saveIAtomContainerSetToSDF(postValidateCompound, inchikey);
					// may not need RemoveDuplication
					// CheminformaticUtility.saveIAtomContainerSetToSDF(GlycosiderUtility.RemoveDupliation(transformedSmiles),inchikey);
					total_num = postValidateCompound.getAtomContainerCount();
				}
			}
			
			
			
			
			
		} catch (InvalidSmilesException e) {
			 e.printStackTrace();
		} catch (Exception e) {
			 e.printStackTrace();
		}
		
		return total_num;
		
		
	}
	
	
	
	/**
	 * 
	 * @param core Number of CPU core to utilize
	 * @param filename The file path. The file should contain smiles
	 * @param exception_time The time limit for each job
	 * @throws IOException 
	 */
	public void runPlantMatMultiThreading(int core, String file_name, int exception_time) throws IOException {
		HashMap<String,Future<Integer>> maps = new HashMap<String,Future<Integer>>();
		
		int numberCore = Integer.valueOf(core);
		
		if(file_name == null) {
			System.out.println("Set the file name (superbio;allHuman)");
			System.exit(0);
		}
		
		
		String current_dir = System.getProperty("user.dir");
		
		System.out.println("number of core: " + Integer.toString(numberCore));
		ExecutorService executor = Executors.newFixedThreadPool(numberCore);
		
		CSVWriter StatusWriter = new CSVWriter(new FileWriter(String.format("%s/generatedfolder/%s", current_dir,"PlantMatException.csv")));
		final CSVParser parser = new CSVParserBuilder().withSeparator('\t').withIgnoreQuotations(true).build();
		final CSVReaderBuilder builder = new CSVReaderBuilder(new FileReader(String.format("%s/generatedfolder/%s", current_dir,file_name)));
		final CSVReader csvReader = builder.withCSVParser(parser).build();
		
		
		String[] tmp;
		while ((tmp = csvReader.readNext()) != null) {
			String compound_smiles = tmp[1];
			String compound_inchikey = tmp[2];

			// check if the file already exist:
			
			File existing_file = new File(String.format("%s/generatedSDF/%s.sdf", current_dir,compound_inchikey));
			if(existing_file.exists()) {
				// continue is kind of like goto `[is] tmp = csvReader.readNext()) != null`
				continue;
			}
			
			Future<Integer> future = executor.submit(new CallablePlantMat(compound_smiles,compound_inchikey));
			maps.put(compound_inchikey, future);
		}
		
		
		for(String key : maps.keySet()) {
			try {
				// https://docs.oracle.com/javase/6/docs/api/java/util/concurrent/Future.html#get%28long,%20java.util.concurrent.TimeUnit%29
				Integer result = maps.get(key).get(exception_time, TimeUnit.SECONDS);
				StatusWriter.writeNext(new String[] {key,Integer.toString(result)});
				
			}catch (TimeoutException e) {
				maps.get(key).cancel(true);
				StatusWriter.writeNext(new String[] {key,"TimeoutException"});
				
			}catch (InterruptedException e) {                                               
				StatusWriter.writeNext(new String[] {key,"InterruptedException"});
			} catch (ExecutionException e) {
				e.printStackTrace();
				StatusWriter.writeNext(new String[] {key,"ExecutionException"});
			}catch (Exception e) {
				e.printStackTrace();
				StatusWriter.writeNext(new String[] {key,"Exception"});
			}
		}
		

		executor.shutdownNow();
		csvReader.close();
		StatusWriter.close();
	}
	
	
	/**
	 * 
	 * @param core Number of CPU core to utilize
	 * @param filename The file path. The file should contain smiles
	 * @param exception_time The time limit for each job
	 * @throws IOException 
	 */
	public void runPlantMat(String file_name) throws IOException {
		
		PlantGlycosider plantglycosider = new PlantGlycosider();
		// molecule validation
		// TODO: implement pre validation by moleculer selection based on ml model, chemical property
		
		
		if(file_name == null) {
			System.out.println("Set the file name (superbio;allHuman)");
			System.exit(0);
		}
		
		
		String current_dir = System.getProperty("user.dir");
		
		CSVWriter StatusWriter = new CSVWriter(new FileWriter(String.format("%s/generatedfolder/%s", current_dir,"PlantMatException.csv")));
		final CSVParser parser = new CSVParserBuilder().withSeparator('\t').withIgnoreQuotations(true).build();
		final CSVReaderBuilder builder = new CSVReaderBuilder(new FileReader(String.format("%s/generatedfolder/%s", current_dir,file_name)));
		final CSVReader csvReader = builder.withCSVParser(parser).build();
		
		
		String[] tmp;
		while ((tmp = csvReader.readNext()) != null) {
			
			String compound_smiles = tmp[1];
			String compound_inchikey = tmp[2];
			
			File existing_file = new File(String.format("%s/generatedSDF/%s.sdf", current_dir,compound_inchikey));
			if(existing_file.exists()) {
				// continue is kind of like goto `[is] tmp = csvReader.readNext()) != null`
				continue;
			}
			
			try {
				int num = plantglycosider.TransformerFoodbSingleCompound(compound_smiles, compound_inchikey);
				StatusWriter.writeNext(new String[] {compound_inchikey,Integer.toString(num)});
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
				StatusWriter.writeNext(new String[] {compound_inchikey,"ProgramException"});
			}
			// check if the file already exist:
		}
		
		csvReader.close();
		StatusWriter.close();
	}
	
	

	
	/**
	 * IAtomContainer maintains a list of Atoms and ElectronContainers (i.e. a
	 * structure) parse smirks get the smirkReaction use smirkManager and
	 * smirkReaction and reactant to produce set of products smirks tutorial
	 * http://daylight.com/dayhtml/doc/theory/theory.smirks.html // looking for the
	 * alocohol functional group add glucose to the OH group input csv should
	 * contain source_id, food_id, original_food_id, food_name (add to corresponding
	 * food is the most important) output csv should contain same field
	 * 
	 * @param args
	 * @throws Exception
	 */
	public static void main(String[] args) throws Exception {
		
		PlantGlycosider fdb = new PlantGlycosider();
									// int core, String file_name, int exception_time;
		fdb.runPlantMatMultiThreading(Integer.valueOf(args[0]), args[1], Integer.valueOf(args[2]));
//		fdb.runPlantMat(args[0]);
		
	}
	
	
	
	
	/**
	 * concurrent class
	 * @author xuan
	 *
	 */
	private static class CallablePlantMat implements Callable<Integer> {
		private String Smiles = new String();
		private String inchikey = new String();

		
		public CallablePlantMat(String Smiles, String inchikey)
		{
			this.Smiles = Smiles;
			this.inchikey = inchikey;             
		
		}          
		
		public Integer call() throws Exception{
			try {
				PlantGlycosider plantglycosider = new PlantGlycosider();
				
				// molecule validation
				// TODO: implement pre validation by moleculer selection based on ml model, chemical property
				int num = plantglycosider.TransformerFoodbSingleCompound(this.Smiles, this.inchikey);
				// molecule selections from transformered metabolites
				// TODO: implement post validation based on selecting the proper one (e.g. no more than five sugar etc.
				
				return num;
			}
			catch(Exception e) {
				return 0;
			}
			 
	              
		}
	}


}

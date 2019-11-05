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

package WishartLab.FoodbScript;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Reader;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.Stack;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.TimeoutException;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.inchi.InChIGenerator;
import org.openscience.cdk.inchi.InChIGeneratorFactory;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.io.SDFWriter;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;

import com.opencsv.CSVReader;
import com.opencsv.CSVWriter;





// run as java -jar PlantMat.jar <core_number> <file_name> <running time>

public class PlantMat {
	
	
	private String current_dir = System.getProperty("user.dir");
	private SmilesGenerator smiGen = new SmilesGenerator(SmiFlavor.Isomeric);
	/**
	 * Sugar's stereochemistry won't lost
	 * if OH number is larger than number of sugar, then permutated all possible add-on for missing sugars
	 * will be called by AddSugarGroupPermu; AddSugarGroupPermu will be called by PreformTransformation
	 * 1. get list of OH
	 * 2. do permutation of OH
	 * @param mole
	 * @throws Exception
	 */
	public ArrayList<String> PreformTransformation(IAtomContainer mole) throws Exception {
		ArrayList<String> result_smiles = new ArrayList<String>();
		SmilesGenerator smiGen = new SmilesGenerator(SmiFlavor.Isomeric);
		Utility combine_util = new Utility();
		Utility new_uti = new Utility();
		
		
		Set<Integer> OH_list = Utility.GetOHList(mole); // return list of OH
		ArrayList<Integer[]> combined_list = new ArrayList<Integer[]>();
		for(int i =  1; i <= OH_list.size(); i++) {
			new_uti.permutations(OH_list, new Stack<Integer>(), i);
			combined_list.addAll(new_uti.get_combination_list());
		}
		
		// remove duplicates
		Set<Integer[]> set = new HashSet<Integer[]>(combined_list);
		combined_list.clear();
		combined_list.addAll(set);
		// remove duplicates
		
		
		
		// for each combination of OH poisition, get all combination of sugar
		Sugars sugar = new Sugars();
		Integer[] sugar_list = sugar.GetSugarIndexSet();
		HashMap<Integer,String> sugarmap = sugar.SugarGroup();
		HashMap<Integer,Integer> sugerindex = sugar.SugarIndex();
		
		
		if (sugar_list.length < OH_list.size()) {
			// when numbers of OH number is larger than number of sugar, 
			// add permutated all possible add-on for missing sugars
			// sugar_list will be mutated 
			// algorithm for mutating sugar_list
			int extra_space = OH_list.size() - sugar_list.length;
			Utility extrace_space_util = new Utility();
			ArrayList<Integer[]> extract_combination = extrace_space_util.GetAllPossibleSet(sugar_list,extra_space);
			// remove the array that is smaller than required space
			for (Integer[] tmp : extract_combination) {
				if (extra_space == tmp.length) {
					Integer [] new_sugar_list = new Integer[tmp.length + extra_space];
					System.arraycopy(sugar_list, 0, new_sugar_list, 0, sugar_list.length);
					for(int i = 0; i < tmp.length; i++) {
						new_sugar_list[sugar_list.length+i] = tmp[i];
					}
					
					// do the permutations
					for(Integer[] l : combined_list) {
						ArrayList<Integer[]> temp_sugar_list = combine_util.GetAllPossibleSet(sugar_list,l.length); // temp set will be list of all possible of order of sugars
						int[] oh_list = new int[l.length];
						for(int i=0;i<l.length;i++) {oh_list[i]=l[i];}
						for(Integer[] sugar_comb : temp_sugar_list) {
							// sugar_comb contain the index of selected sugar
							String[] smiles = new String[sugar_comb.length];
							int[] smiles_index = new int[sugar_comb.length];
							
							for (int i=0;i<smiles.length;i++) {smiles[i] = sugarmap.get(sugar_comb[i]);}
							for (int i=0;i<sugar_list.length;i++) {smiles_index[i] = sugerindex.get(sugar_comb[i]);}
							IAtomContainer mole_man = AddSugarGroupPermu(mole,oh_list,smiles,smiles_index);
							result_smiles.add(Utility.SplitContainer(smiGen.create(mole_man)));
							
						}
						
					}
					
				}
				
			}
			
		}else {
			// (sugar_list.length > list.size())
			// don't need to worry about previous situation
			// combined_list is the combination of all possible order of OH
			// sugar_list is the combination of sugars
			
			// for each combined_list (list of array of OH position)
			// add each sugar one at time
			for(Integer[] l : combined_list) {
//				System.out.println("Current combination:: " + Arrays.toString(l));
				
				ArrayList<Integer[]> temp_sugar_list = combine_util.GetAllPossibleSet(sugar_list,l.length); // temp set will be list of all possible of order of sugars
				int[] oh_list = new int[l.length];
				for(int i=0;i<l.length;i++) {oh_list[i]=l[i];}
				for(Integer[] sugar_comb : temp_sugar_list) {
					// sugar_comb contain the index of selected sugar
					String[] smiles = new String[sugar_comb.length];
					int[] smiles_index = new int[sugar_comb.length];
					
					for (int i=0;i<smiles.length;i++) {smiles[i] = sugarmap.get(sugar_comb[i]);}
					for (int i=0;i<smiles_index.length;i++) {smiles_index[i] = sugerindex.get(sugar_comb[i]);}
					IAtomContainer copied_mole = mole.clone();
					IAtomContainer mole_man = AddSugarGroupPermu(copied_mole,oh_list,smiles,smiles_index);
//					System.out.println("PreformTransformation => "+Utility.SplitContainer(smiGen.create(mole_man)));
					result_smiles.add(Utility.SplitContainer(smiGen.create(mole_man)));
				}
				
			}
			
		}
		return result_smiles;
		
		
	}
	
	/**
	 * combine two IAtomContainer, they must have charge on the atom of interests
	 * 
	 * order: add sugar to mole
	 * @param mole
	 * @param sugar
	 * @return
	 * @throws CDKException 
	 */
	public IAtomContainer CombineContainers(IAtomContainer mole, double[] charge_table, IAtomContainerSet sugar) throws CDKException {

		for(int c = 0; c < charge_table.length; c++) {
			mole.add(sugar.getAtomContainer(c));
			int index_1 = new Integer(0);
			int index_2 = new Integer(0);
			double charge = charge_table[c];
			for (int k = 0; k < mole.getAtomCount(); k++) {
				
				IAtom atoms = mole.getAtom(k);
				if (atoms.getCharge() != null) {					
					// find sugar molecule mark
					// with garantee of number of sugar can't exceed the number of OH, add sugar one by one
					// need to charge for order
					if (atoms.getCharge() == 1.0) {
						index_1 = k;
						atoms.setCharge(null);
					} 
					
					// find original molecule mark
					else if (atoms.getCharge() == charge) {
						List<IAtom> nearestAtom = mole.getConnectedAtomsList(atoms);
						for (int i = 0; i < nearestAtom.size(); i++) {
							if(nearestAtom.get(i).getSymbol().equals("C")) {
								int real_connect = mole.indexOf(nearestAtom.get(i));
								index_2 = real_connect;
								mole.removeBond(atoms, nearestAtom.get(i));
							}
						}
						
						atoms.setCharge(null);
					}
				}
				
				// add bond / make connection
				
				// remove sugar's charge
				
			}
			mole.addBond(index_1, index_2, IBond.Order.SINGLE);
			for (int k = 0; k < mole.getAtomCount(); k++) {
				IAtom atoms = mole.getAtom(k);
				if (atoms.getCharge() != null) {
					if(atoms.getCharge() == 1.0) {
						atoms.setCharge(null);
					}
					
				}
			}
			
		}
		
		// remove the charges
		for (int k = 0; k < mole.getAtomCount(); k++) {
			IAtom atoms = mole.getAtom(k);
			if (atoms.getCharge() != null) {
				atoms.setCharge(null);
			}
		}
		
		return mole;
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
//			System.out.println("One Sugar");
			String[] new_sugar_smiles = new String[mole_index.length];
			int[] new_sugar_index = new int[mole_index.length];
			for(int i=0; i < mole_index.length; i++) {
				new_sugar_smiles[i] = sugar_smiles[0];
				new_sugar_index[i] = sugar_index[0];
				
			}
//			System.out.println(String.format("new_sugar_smiles length => %d; new_sugar_index length => %d;", 
//					new_sugar_smiles.length,new_sugar_index.length));
			// add one sugar at time
			// increment charges start from two
			// keep all sugar charges one
			// connect the sugar and compound with order
			double [] charge_table = new double[mole_index.length];
			
			IAtomContainerSet sugar_mole_set = builder.newInstance(IAtomContainerSet.class);
			
			for(int i=0; i < mole_index.length;i++) {
				double charge = (double) i+2;
//				System.out.println("mole_index[i] => " + mole_index[i]);
				mole.getAtom(mole_index[i]).setCharge(charge);
				charge_table[i] = charge;
				sugar_mole_set.addAtomContainer(Utility.TransformSmilesToContainer(new_sugar_smiles[i],new_sugar_index[i]));
			}
			
//			for (int k = 0; k < mole.getAtomCount(); k++) {
//				IAtom atoms = mole.getAtom(k);
//				System.out.println("all atom charge => " + atoms.getCharge());
//			}
			try {
				output = CombineContainers(mole,charge_table,sugar_mole_set);
//				System.out.println("AddSugarGroupPermu => "+smiGen.create(mole));
				// after CombineContainers, all charges should be removed;
			} catch (InvalidSmilesException e) {
				e.printStackTrace();
			} catch (CDKException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			
		}else if(sugar_smiles.length > 1) {
			// following the order
			// add whatever how many OH contains.
			// when sugar length is less than OH index
//			System.out.println("More Than One Sugar");
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
				sugar_mole_set.addAtomContainer(Utility.TransformSmilesToContainer(new_sugar_smiles[i],new_sugar_index[i]));
			}
			for(int i=sugar_smiles.length-1; i < mole_index.length; i++) {
				new_sugar_smiles[i] = last_smiles;
				new_sugar_index[i] = last_sugar_index;
				sugar_mole_set.addAtomContainer(Utility.TransformSmilesToContainer(last_smiles,last_sugar_index));
			}
			
			double [] charge_table = new double[mole_index.length];
			
			for(int i=0; i < mole_index.length;i++) {
				double charge = (double) i+2;
				mole.getAtom(mole_index[i]).setCharge(charge);
				charge_table[i] = charge;				
//				System.out.println("new_sugar_smiles[i] =>" + new_sugar_smiles[i]);
			}
			
			try {
				output = CombineContainers(mole,charge_table,sugar_mole_set);
//				System.out.println("AddSugarGroupPermu => "+smiGen.create(mole));
				// after CombineContainers, all charges should be removed;
			} catch (InvalidSmilesException e) {
				e.printStackTrace();
			} catch (CDKException e) {
				// TODO Auto-generated catch block
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
	public void TransformerFoodbSingleCompound(String smiles,String inchikey) throws IOException, InterruptedException {
		
		
		PlantMat fdb = new PlantMat();
		IChemObjectBuilder builder = SilentChemObjectBuilder.getInstance(); // access singleton instance of factory
		SmilesParser smiParser = new SmilesParser(builder);
		ArrayList<String> Exception = new ArrayList<String>(); // should store inchikey
		
		try {
			ArrayList<String> transformedSmiles = fdb.PreformTransformation(smiParser.parseSmiles(smiles));			
			ArrayList<String> non_duplicate_transformedSmiles = RemoveDupliation(transformedSmiles);
			
			File file = new File(String.format("%s/generatedfolder/%s.csv", current_dir,inchikey));
			FileWriter outputfile = new FileWriter(file); 
			CSVWriter writer = new CSVWriter(outputfile);
			

	        FileWriter fw = new FileWriter(String.format("%s/generatedfolder/%s.sdf", current_dir,inchikey), true);
	        SDFWriter sdfwriter = new SDFWriter(fw);			        
	        
	        for(int i = 0; i < non_duplicate_transformedSmiles.size(); i++) {
				writer.writeNext(new String[]{non_duplicate_transformedSmiles.get(i)});
				IAtomContainer mole = Utility.Generate2DCoordinate(smiParser.parseSmiles(non_duplicate_transformedSmiles.get(i)));
				if (mole != null) {
					sdfwriter.write(mole);
				}
				
			}
	        
			writer.close();
			sdfwriter.close();
	        fw.close();
			
			
			
		} catch (InvalidSmilesException e) {
			// TODO Auto-generated catch block
			 e.printStackTrace();
			Exception.add(inchikey);
		} catch (Exception e) {
			// TODO Auto-generated catch block
			 e.printStackTrace();
			Exception.add(inchikey);
		}
		
		if(Exception.size()>0) {
			File file = new File(String.format("%s/generatedfolder/compound_transformation_Exception.csv",current_dir));
			FileWriter outputfile = new FileWriter(file); 
			CSVWriter writer = new CSVWriter(outputfile);
			for(int i = 0; i < Exception.size(); i++) {
				writer.writeNext(new String[] {Exception.get(i)});				
			}
			writer.close();
		}
		
		
	}
	
	/**
	 * parse and transformer all compounds from foodb
	 * get all transformerable compounds and save them to directory
	 * then get these file add to new food
	 * @param data_file
	 * @throws IOException 
	 * @throws InterruptedException 
	 */
	public void TransformerFoodbCompound(String data_file) throws IOException, InterruptedException {
		
//		final CSVParser parser = new CSVParserBuilder().withSeparator('\t').build();
//		final CSVReader reader = new CSVReaderBuilder(new StringReader(data_file)).withCSVParser(parser).build();
		PlantMat fdb = new PlantMat();
		IChemObjectBuilder builder = SilentChemObjectBuilder.getInstance(); // access singleton instance of factory
		SmilesParser smiParser = new SmilesParser(builder);
		Reader reader = Files.newBufferedReader(Paths.get(data_file));
        CSVReader csvReader = new CSVReader(reader);
        ArrayList<String> transformed_holder = new ArrayList<String>();
        ArrayList<String> Exception = new ArrayList<String>(); // should store inchikey
        String[] tmp;
		while ((tmp = csvReader.readNext()) != null) {
			System.out.println(tmp[0]);
			String[] line = tmp[0].split("\t");
			
			String food_id = line[0].replace("\"", "");
			String compound_smiles = line[1].replace("\"", "");
			String compound_inchikey = line[2].replace("\"", "");
			
			if(transformed_holder.contains(compound_inchikey)) {
				// already did transformation
				// just extract the file that contain all the file
				continue;
			}
			else {
				transformed_holder.add(compound_inchikey);
				try {
					ArrayList<String> transformedSmiles = fdb.PreformTransformation(smiParser.parseSmiles(compound_smiles));
//					System.out.println(String.format("number of transformation => %d", transformedSmiles.size()));
					
					ArrayList<String> non_duplicate_transformedSmiles = RemoveDupliation(transformedSmiles);
//					System.out.println(String.format("number of transformation after remove=> %d", non_duplicate_transformedSmiles.size()));
					
					File file = new File(String.format("%s/generatedfolder/%s.csv", current_dir,compound_inchikey));
					FileWriter outputfile = new FileWriter(file); 
					CSVWriter writer = new CSVWriter(outputfile);
					

			        FileWriter fw = new FileWriter(String.format("%s/generatedfolder/%s.sdf", current_dir,compound_inchikey), true);
			        SDFWriter sdfwriter = new SDFWriter(fw);			        
			        
			        for(int i = 0; i < non_duplicate_transformedSmiles.size(); i++) {
						writer.writeNext(new String[]{non_duplicate_transformedSmiles.get(i)});
						IAtomContainer mole = Utility.Generate2DCoordinate(smiParser.parseSmiles(non_duplicate_transformedSmiles.get(i)));
						if (mole != null) {
							sdfwriter.write(mole);
						}
						
					}
			        
					writer.close();
					sdfwriter.close();
			        fw.close();
					
					
					
				} catch (InvalidSmilesException e) {
					// TODO Auto-generated catch block
					 e.printStackTrace();
					Exception.add(compound_inchikey);
				} catch (Exception e) {
					// TODO Auto-generated catch block
					 e.printStackTrace();
					Exception.add(compound_inchikey);
				}
				
			}
			
//			TimeUnit.SECONDS.sleep(2);
			break;
		}
		
		if(Exception.size()>0) {
			File file = new File(String.format("%s/generatedfolder/compound_transformation_Exception.csv",current_dir));
			FileWriter outputfile = new FileWriter(file); 
			CSVWriter writer = new CSVWriter(outputfile);
			for(int i = 0; i < Exception.size(); i++) {
				writer.writeNext(new String[] {Exception.get(i)});				
			}
			writer.close();
		}
		
		
		csvReader.close();
	}
	
	
	/**
	 * remove duplication based on the inchikey
	 * @param transformedSmiles
	 * @return
	 * @throws CDKException 
	 */
	public ArrayList<String> RemoveDupliation(ArrayList<String> transformedSmiles) throws CDKException {
		ArrayList<String> non_dup = new ArrayList<String>();
		ArrayList<String> non_dup_smiles = new ArrayList<String>();
		InChIGeneratorFactory factory = InChIGeneratorFactory.getInstance();
		IChemObjectBuilder builder = SilentChemObjectBuilder.getInstance();
		SmilesParser smiParser = new SmilesParser(builder);
		
		for(int i = 0; i < transformedSmiles.size(); i++) {
			try {
				InChIGenerator gen = factory.getInChIGenerator(smiParser.parseSmiles(transformedSmiles.get(i)));
				String inchikey = gen.getInchiKey();
				if (non_dup.contains(inchikey)) {
					continue;
				}else {
					non_dup.add(inchikey);
					non_dup_smiles.add(transformedSmiles.get(i));
				}
			} catch (InvalidSmilesException e) {
				// TODO Auto-generated catch block
				
				e.printStackTrace();
				continue;
			} catch (CDKException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
				continue;
			}
		}
		
		return non_dup_smiles;
	}
	
	
	
	
	/**
	 * 
	 * @param core Number of CPU core to utilize
	 * @param filename The file path. The file should contain smiles
	 * @param exception_time The time limit for each job
	 * @throws IOException 
	 */
	public void runPlantMatMultiThreading(int core, String file_name, int exception_time) throws IOException {
		HashMap<String,Future<Boolean>> maps = new HashMap<String,Future<Boolean>>();
		
		int numberCore = Integer.valueOf(core);
		
		if(file_name == null) {
			System.out.println("Set the file name (superbio;allHuman)");
			System.exit(0);
		}
		
		System.out.println(String.format("Current exception time is %d", exception_time));
		
		String current_dir = System.getProperty("user.dir");
		
		System.out.println("number of core: " + Integer.toString(numberCore));
		ExecutorService executor = Executors.newFixedThreadPool(numberCore);
		
		CSVReader csvReader = new CSVReader(new FileReader(String.format("%s/generatedfolder/%s", current_dir,file_name)));
		CSVWriter StatusWriter = new CSVWriter(new FileWriter(String.format("%s/generatedfolder/%s", current_dir,"PlantMatException.csv")));
//		csvReader.readNext(); // skip header
		String[] tmp;
		while ((tmp = csvReader.readNext()) != null) {
			String[] line = tmp[0].split("\t");
			String compound_smiles = line[1].replace("\"", "");
			String compound_inchikey = line[2].replace("\"", "");
			Future<Boolean> future = executor.submit(new CallablePlantMat(compound_smiles,compound_inchikey));
            maps.put(compound_inchikey, future);
		}
		
		
		for(String key : maps.keySet()) {
			try {
				// https://docs.oracle.com/javase/6/docs/api/java/util/concurrent/Future.html#get%28long,%20java.util.concurrent.TimeUnit%29
				boolean result = maps.get(key).get(exception_time, TimeUnit.SECONDS);
				if(!result) { StatusWriter.writeNext(new String[] {key,"ProgramException"});}
				
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
		
		PlantMat fdb = new PlantMat();
		fdb.runPlantMatMultiThreading(Integer.valueOf(args[0]), args[1], Integer.valueOf(args[2]));
		
	}
	
	
	
	/**
	 * concurrent class
	 * @author xuan
	 *
	 */
	private static class CallablePlantMat implements Callable<Boolean> {
		private String Smiles = new String();
		private String inchikey = new String();
		
		public CallablePlantMat(String Smiles, String inchikey)
		{
			this.Smiles = Smiles;
			this.inchikey = inchikey;             
		
		}          
		
		public Boolean call() throws Exception{
			try {
				PlantMat planamat = new PlantMat();
				// molecule validation
				planamat.TransformerFoodbSingleCompound(this.Smiles, this.inchikey);
				// molecule selections from transformered metabolites
				return true;
			}
			catch(Exception e) {
				return false;
			}
			
			 
	              
		}
	}
	
	



}

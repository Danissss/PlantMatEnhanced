package WishartLab.PlantGlycosider;

import java.util.ArrayList;
import java.util.List;

import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import ambit2.smarts.query.SmartsPatternCDK;
/**
 * 
 * @author xuan
 * This class is not used; 
 * This class is only for storing used code;
 */
public class LegacyFunction {
	
	
	// String smirks = "[H:3][#8:2]-[#6:1]>>[H][#8]C([H])([H])[C@]1([H])[#8][C@@]([H:3])([#8:2][C:1]([H])([H])[H])[C@@]([H])([#8][H])[C@]([H])([#8][H])[C@@]1([H])[#8][H]";
	
	// monosaccharides
	// hexoses; n = 6
	final static String d_fructose = "C1[C@H]([C@H]([C@@H](C(O1)(CO)O)O)O)O";
	final static String d_glucose = "C([C@@H]1[C@H]([C@@H]([C@H](C(O1)O)O)O)O)O";
	final static String d_galactose = "C([C@@H]1[C@@H]([C@@H]([C@H](C(O1)O)O)O)O)O";
	final static int d_fructose_index = 7;
	final static int d_glucose_index = 11;
	final static int d_galactose_index = 11;
	// pentoses; n = 5;
//	final static String d_ribose = "C1[C@H]([C@H]([C@H](C(O1)O)O)O)O";
//	final static String l_ribose = "C([C@@H]([C@@H]([C@@H](C=O)O)O)O)O"; // not in ring format
//	final static String d_glucose = "C([C@@H]1[C@H]([C@@H]([C@H](C(O1)O)O)O)O)O";
//	final static String d_galactose = "C([C@@H]1[C@@H]([C@@H]([C@H](C(O1)O)O)O)O)O";
//	final static int d_ribose_index = 7;
//	final static int l_ribose_index = 0;
//	final static int d_galactose_index = 11;
	
	// disaccharides 
	// delphidin dihexose 
	final static String sucrose = "C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O[C@]2([C@H]([C@@H]([C@H](O2)CO)O)O)CO)O)O)O)O";
	final static String maltose = "C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O[C@@H]2[C@H](O[C@@H]([C@@H]([C@H]2O)O)O)CO)O)O)O)O";
	final static String lactose = "C([C@@H]1[C@@H]([C@@H]([C@H]([C@@H](O1)O[C@@H]2[C@H](OC([C@@H]([C@H]2O)O)O)CO)O)O)O)O";
	final static int sucrose_index = 14;
	final static int maltose_index = 18;
	final static int lactose_index = 18;
	
	
	
	
	/**
	 * Given the transformer molecule and sugar type, and site of changes
	 * do biotransformation 
	 * legacy function
	 * 
	 * @param mole
	 * @param enzyme_name
	 * @param som_site
	 * @return
	 * @throws Exception
	 */
	public static IAtomContainer AddSugarGroup(IAtomContainer mole, String sugar_name, int som_site)
			throws Exception {
		// System.out.println("-----------------------------------------------------------");
		
		SmilesParser sp = new SmilesParser(SilentChemObjectBuilder.getInstance());
		IAtomContainer sugar_mole = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainer.class);
		if(sugar_name.equals("sucrose")) {
			// convert smiles to IAtomContainer 
			// add charge to the right atom index
			// when combine two IAtomContainers, 
			sugar_mole = sp.parseSmiles(sucrose);
			sugar_mole.getAtom(sucrose_index).setCharge(1.0);
		}
		else if (sugar_name.equals("maltose")) {
			sugar_mole = sp.parseSmiles(maltose);
			sugar_mole.getAtom(maltose_index).setCharge(1.0);
		}
		else if (sugar_name.equals("lactose")) {
			sugar_mole = sp.parseSmiles(lactose);
			sugar_mole.getAtom(lactose_index).setCharge(1.0);
		}
		else if (sugar_name.equals("d_fructose")) {
			sugar_mole = sp.parseSmiles(d_fructose);
			sugar_mole.getAtom(d_fructose_index).setCharge(1.0);
		}
		else if (sugar_name.equals("d_glucose")) {
			sugar_mole = sp.parseSmiles(d_glucose);
			sugar_mole.getAtom(d_glucose_index).setCharge(1.0);
		}
		else if (sugar_name.equals("d_galactose")) {
			sugar_mole = sp.parseSmiles(d_galactose);
			sugar_mole.getAtom(d_galactose_index).setCharge(1.0);
		}

		
		// below is adding step
		mole.getAtom(som_site).setCharge(2.0);
		mole.add(sugar_mole);

		int index_1 = new Integer(0);
		int index_2 = new Integer(0);
		for (int k = 0; k < mole.getAtomCount(); k++) {
			IAtom atoms = mole.getAtom(k);
			if (atoms.getCharge() != null) {
				
				
				// find sugar molecule mark
				if (atoms.getCharge() == 1.0) {
					index_1 = k;
					atoms.setCharge(null);
					
					
				} 
				
				// find original molecule mark
				else if (atoms.getCharge() == 2.0) {
					List<IAtom> nearestAtom = mole.getConnectedAtomsList(atoms);
					for (int i = 0; i < nearestAtom.size(); i++) {
						// System.out.println(nearestAtom.get(i).getSymbol());
						if(nearestAtom.get(i).getSymbol().equals("C")) {
							int real_connect = mole.indexOf(nearestAtom.get(i));
							// System.out.println(real_connect);
							index_2 = real_connect;
							mole.removeBond(atoms, nearestAtom.get(i));
							// mole.addBond(index_1, index_2, IBond.Order.SINGLE);
						}
					}
					
					// index_2 = k;
					atoms.setCharge(null);
				}
			}
		}

		// add bond / make connection
		mole.addBond(index_1, index_2, IBond.Order.SINGLE);

		// remove the charges
		for (int k = 0; k < mole.getAtomCount(); k++) {
			IAtom atoms = mole.getAtom(k);
			if (atoms.getCharge() != null) {
				// System.out.println(atoms.getCharge());
				atoms.setCharge(null);
				// System.out.println(atoms.getCharge());
			}
		}
		 
		return mole;
	}
	
	/**
	 * Detect if the compound has OH group for attaching sugar compound
	 * Sugar's stereochemistry won'y lost
	 * @param mole
	 * @throws Exception
	 */
	public ArrayList<String> PreformTransformationToString(IAtomContainer mole, String sugar) throws Exception {
		SmilesGenerator smiGen = new SmilesGenerator(SmiFlavor.Isomeric);
		SmartsPatternCDK smarts = new SmartsPatternCDK();
		ArrayList<String> transformered_list = new ArrayList<String>();
		smarts.setSmarts("[OH]");
		
		// smarts.hasSMARTSPattern(mole) return number of encountered smarts
		if(smarts.hasSMARTSPattern(mole) != 0) {
			List<List<Integer>> matach_atom_ind = smarts.getUniqueMatchingAtoms(mole); // two matching group will be list.size() = 2
			for(int i = 0; i< matach_atom_ind.size(); i++) {
				List<Integer> tmp = matach_atom_ind.get(i);			// usually the first index is the som;
				int som = tmp.get(0);
				IAtomContainer mole_copy = mole.clone();
				
				// do transformation for every type of sugar
				IAtomContainer mol = AddSugarGroup(mole_copy, sugar,som);
				String unsplited_mol = smiGen.create(mol);
				
				String[] splited_mol = unsplited_mol.split("\\.");
				if (splited_mol.length == 2) {
					if (splited_mol[0].length() > splited_mol[1].length()) {
						System.out.println(splited_mol[0]);
						transformered_list.add(splited_mol[0]);
					}
					else {
						System.out.println(splited_mol[1]);
						transformered_list.add(splited_mol[1]);
					}
				}

			}
		}
		return transformered_list;
		
	}
	

	
	
	/**
	 * do transformation
	 * 
	 * @param mol
	 * @param smirks
	 * @return
	 */
//	public static IAtomContainerSet ApplyTransformation(IAtomContainer mol, String smirks) {
//		IChemObjectBuilder builder = SilentChemObjectBuilder.getInstance(); // access singleton instance of factory
//																			// class
//		SMIRKSManager smrkMan = new SMIRKSManager(builder); // constructor
//		SMIRKSReaction smirksReaction = smrkMan.parse(smirks); // parses a SMIRKS string
//
//		AtomContainerManipulator.convertImplicitToExplicitHydrogens(mol); // adds explicit hydrogens
//		try {
//			IAtomContainerSet products = smrkMan.applyTransformationWithSingleCopyForEachPos(mol, null, smirksReaction);
//			return products;
//		} catch (CDKException e) {
//			// TODO Auto-generated catch block
//			e.printStackTrace();
//			return null;
//		} catch (Exception e) {
//			// TODO Auto-generated catch block
//			e.printStackTrace();
//			return null;
//		}
//
//	}
	
	/**
	 * print the product for debug
	 * 
	 * @param products
	 * @param smiles_input
	 * @throws CDKException
	 */
	public static void print_atomcontainer_list(IAtomContainerSet products, String smiles_input) throws CDKException {
		CDKHydrogenAdder adder = CDKHydrogenAdder.getInstance(SilentChemObjectBuilder.getInstance());
		if (products != null) {
			System.out.println(String.format("Input SMILES: %s", smiles_input));
			System.out.println("Product Count: " + products.getAtomContainerCount()); // prints count of products
			for (IAtomContainer product : products.atomContainers()) { // each iteration of products.atomContainers()
																		// stored in the variable |product|
				adder.addImplicitHydrogens(product);
				IAtomContainer noHydrogenProduct = AtomContainerManipulator.suppressHydrogens(product);
				SmilesGenerator smiGen = new SmilesGenerator(SmiFlavor.Isomeric);
				System.out.println("SMILES: " + smiGen.create(noHydrogenProduct)); // outputs isomeric SMILES of product
			}
		}
		System.out.println("=======================================================");
	}
	
	
	
	public static IAtomContainer TransformSmilesToContainer(String sugar_name, int suger_index) throws InvalidSmilesException {
		SmilesParser sp = new SmilesParser(SilentChemObjectBuilder.getInstance());
		IAtomContainer sugar_mole = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainer.class);
		if(sugar_name.equals("sucrose")) {
			// convert smiles to IAtomContainer 
			// add charge to the right atom index
			// when combine two IAtomContainers, 
			sugar_mole = sp.parseSmiles(sucrose);
			sugar_mole.getAtom(sucrose_index).setCharge(1.0);
		}
		else if (sugar_name.equals("maltose")) {
			sugar_mole = sp.parseSmiles(maltose);
			sugar_mole.getAtom(maltose_index).setCharge(1.0);
		}
		else if (sugar_name.equals("lactose")) {
			sugar_mole = sp.parseSmiles(lactose);
			sugar_mole.getAtom(lactose_index).setCharge(1.0);
		}
		else if (sugar_name.equals("d_fructose")) {
			sugar_mole = sp.parseSmiles(d_fructose);
			sugar_mole.getAtom(d_fructose_index).setCharge(1.0);
		}
		else if (sugar_name.equals("d_glucose")) {
			sugar_mole = sp.parseSmiles(d_glucose);
			sugar_mole.getAtom(d_glucose_index).setCharge(1.0);
		}
		else if (sugar_name.equals("d_galactose")) {
			sugar_mole = sp.parseSmiles(d_galactose);
			sugar_mole.getAtom(d_galactose_index).setCharge(1.0);
		}
		
		return sugar_mole;
	}
	
	
	
//	/**
//	 * parse and transformer all compounds from foodb
//	 * get all transformerable compounds and save them to directory
//	 * then get these file add to new food
//	 * @param data_file
//	 * @throws IOException 
//	 * @throws InterruptedException 
//	 */
//	public void TransformerFoodbCompound(String data_file) throws IOException, InterruptedException {
//		
////		final CSVParser parser = new CSVParserBuilder().withSeparator('\t').build();
////		final CSVReader reader = new CSVReaderBuilder(new StringReader(data_file)).withCSVParser(parser).build();
//		PlantGlycosider fdb = new PlantGlycosider();
//		IChemObjectBuilder builder = SilentChemObjectBuilder.getInstance(); // access singleton instance of factory
//		SmilesParser smiParser = new SmilesParser(builder);
//		Reader reader = Files.newBufferedReader(Paths.get(data_file));
//        CSVReader csvReader = new CSVReader(reader);
//        ArrayList<String> transformed_holder = new ArrayList<String>();
//        ArrayList<String> Exception = new ArrayList<String>(); // should store inchikey
//        String[] tmp;
//		while ((tmp = csvReader.readNext()) != null) {
//			System.out.println(tmp[0]);
//			String[] line = tmp[0].split("\t");
//			
////			String food_id = line[0].replace("\"", "");
//			String compound_smiles = line[1].replace("\"", "");
//			String compound_inchikey = line[2].replace("\"", "");
//			
//			if(transformed_holder.contains(compound_inchikey)) {
//				// already did transformation
//				// just extract the file that contain all the file
//				continue;
//			}
//			else {
//				transformed_holder.add(compound_inchikey);
//				try {
//					ArrayList<String> transformedSmiles = fdb.PreformTransformation(smiParser.parseSmiles(compound_smiles));					
//					IAtomContainerSet uniqueContainer = RemoveDupliation(transformedSmiles);
////					System.out.println(String.format("number of transformation after remove=> %d", non_duplicate_transformedSmiles.size()));
//					
//
//			        FileWriter fw = new FileWriter(String.format("%s/generatedfolder/%s.sdf", current_dir,compound_inchikey), true);
//			        SDFWriter sdfwriter = new SDFWriter(fw);			        
//			        
//			        for(int i = 0; i < uniqueContainer.getAtomContainerCount(); i++) {
//						IAtomContainer mole = Utility.Generate2DCoordinate(uniqueContainer.getAtomContainer(i));
//						if (mole != null) {
//							sdfwriter.write(mole);
//						}
//						
//					}
//			        
//					sdfwriter.close();
//			        fw.close();
//					
//					
//					
//				} catch (InvalidSmilesException e) {
//					 e.printStackTrace();
//					Exception.add(compound_inchikey);
//				} catch (Exception e) {
//					 e.printStackTrace();
//					Exception.add(compound_inchikey);
//				}
//				
//			}
//			
////			TimeUnit.SECONDS.sleep(2);
//			break;
//		}
//		
//		if(Exception.size()>0) {
//			File file = new File(String.format("%s/generatedfolder/compound_transformation_Exception.csv",current_dir));
//			FileWriter outputfile = new FileWriter(file); 
//			CSVWriter writer = new CSVWriter(outputfile);
//			for(int i = 0; i < Exception.size(); i++) {
//				writer.writeNext(new String[] {Exception.get(i)});				
//			}
//			writer.close();
//		}
//		
//		
//		csvReader.close();
//	}
	
	
	
// legacy main
//	ArrayList<String> smiles_list = new ArrayList<String>();
//	smiles_list.add("CO");
//	smiles_list.add("OC1=CC=CC2=C1C=CC=C2");
//	smiles_list.add("CC(=C)C(O)=C");
//	smiles_list.add("OC1=CC=CC=C1O");
//	smiles_list.add("OC1=CC=CC=C1C1=CC=CC=C1O");
//	
//	for (int i = 0; i < smiles_list.size(); i++) {
//		IAtomContainer reactant = smiParser.parseSmiles(smiles_list.get(i));
//		if (Utility.is_terpenoids(reactant)) {
//			System.out.println("terpenoids");
//			IAtomContainerSet products = ApplyTransformation(reactant, smirks);
//			print_atomcontainer_list(products, smiles_list.get(i));
//		}
//		if (Utility.is_polyphenol(reactant)) {
//			System.out.println("polyphenol");
//			IAtomContainerSet products = ApplyTransformation(reactant, smirks);
//			print_atomcontainer_list(products, smiles_list.get(i));
//		}
//	}
//
}

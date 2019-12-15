package WishartLab.PlantGlycosider;


import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

import org.apache.commons.collections.map.LinkedMap;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.io.SDFWriter;
import org.openscience.cdk.io.iterator.IteratingSDFReader;
import org.openscience.cdk.qsar.DescriptorEngine;
import org.openscience.cdk.qsar.IDescriptor;
import org.openscience.cdk.qsar.IMolecularDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.ALOGPDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.APolDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.AromaticAtomsCountDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.AromaticBondsCountDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.AtomCountDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.AutocorrelationDescriptorCharge;
import org.openscience.cdk.qsar.descriptors.molecular.AutocorrelationDescriptorMass;
import org.openscience.cdk.qsar.descriptors.molecular.AutocorrelationDescriptorPolarizability;
import org.openscience.cdk.qsar.descriptors.molecular.BCUTDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.BPolDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.BondCountDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.CarbonTypesDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.ChiChainDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.ChiClusterDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.ChiPathClusterDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.ChiPathDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.EccentricConnectivityIndexDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.FragmentComplexityDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.HBondAcceptorCountDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.HBondDonorCountDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.KappaShapeIndicesDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.KierHallSmartsDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.LargestChainDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.LargestPiSystemDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.LongestAliphaticChainDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.PetitjeanNumberDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.PetitjeanShapeIndexDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.RotatableBondsCountDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.RuleOfFiveDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.TPSADescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.VAdjMaDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.WeightDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.WeightedPathDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.WienerNumbersDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.XLogPDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.ZagrebIndexDescriptor;
import org.openscience.cdk.qsar.result.DoubleArrayResult;
import org.openscience.cdk.qsar.result.DoubleResult;
import org.openscience.cdk.qsar.result.IDescriptorResult;
import org.openscience.cdk.qsar.result.IntegerArrayResult;
import org.openscience.cdk.qsar.result.IntegerResult;

import com.opencsv.CSVParser;
import com.opencsv.CSVParserBuilder;
import com.opencsv.CSVReader;
import com.opencsv.CSVReaderBuilder;
import com.opencsv.CSVWriter;


import weka.classifiers.trees.RandomForest;
import weka.core.Attribute;
import weka.core.DenseInstance;
import weka.core.Instances;
import weka.core.SerializationHelper;
import weka.core.converters.ArffSaver;



public class CreateSolubilityModel {
	
	
	private String current_dir = System.getProperty("user.dir");
	private List<String> classNames = DescriptorEngine.getDescriptorClassNameByPackage("org.openscience.cdk.qsar.descriptors.molecular",
            null);
	private DescriptorEngine descriptoEngine = new DescriptorEngine(classNames, null);
	private String modelName = "LogSPredictor2019_11_27";

	
	/**
	 * 
	 * @param mole
	 * @throws Exception
	 * @author xuan
	 */
	public String generateAllCDKMolecularDescriptors(IAtomContainer mole) throws Exception {
		
		List<IDescriptor> descriptors = 	descriptoEngine.getDescriptorInstances();
		
		String temp = "";
		for (IDescriptor desc : descriptors) {
			try {
				IDescriptorResult res = ((IMolecularDescriptor) desc).calculate(mole).getValue();
				if (res instanceof IntegerResult) {
					temp = temp+","+res.toString();
				} else if (res instanceof DoubleResult) {
					temp = temp+","+res.toString();
				} else if (res instanceof DoubleArrayResult) {
					temp = temp+","+res.toString();
				} else if (res instanceof IntegerArrayResult) {
					temp = temp+","+res.toString();
				} else
					throw new IllegalStateException(
							"Unknown idescriptor result value for '" + desc + "' : " + res.getClass());
			} catch (Throwable e) {
				System.err.println("Could not compute cdk feature " + desc);

			}
			
		}
		return temp;
		
	}
	
	
	/**
	 * get only desired descriptors 
	 * similar to List<IDescriptor> descriptors = descriptoEngine.getDescriptorInstances();
	 * @return
	 * @throws CDKException 
	 */
	public ArrayList<IDescriptor> getSelectedDescriptors() throws CDKException{
		ArrayList<IDescriptor> descriptors = new ArrayList<IDescriptor>();
		descriptors.add(new ALOGPDescriptor());
//		descriptors.add(new BCUTDescriptor());
		descriptors.add(new FragmentComplexityDescriptor());
		descriptors.add(new APolDescriptor());
		descriptors.add(new AromaticAtomsCountDescriptor());
		descriptors.add(new AromaticBondsCountDescriptor());
		descriptors.add(new AtomCountDescriptor());
		descriptors.add(new AutocorrelationDescriptorCharge());
		descriptors.add(new AutocorrelationDescriptorMass());
		descriptors.add(new AutocorrelationDescriptorPolarizability());
		descriptors.add(new BondCountDescriptor());
		descriptors.add(new BPolDescriptor());
		descriptors.add(new CarbonTypesDescriptor());
		descriptors.add(new ChiChainDescriptor());
		descriptors.add(new ChiClusterDescriptor());
		descriptors.add(new ChiPathDescriptor());
		descriptors.add(new ChiPathClusterDescriptor());
		descriptors.add(new EccentricConnectivityIndexDescriptor());
		descriptors.add(new HBondDonorCountDescriptor());
		descriptors.add(new HBondAcceptorCountDescriptor());
		descriptors.add(new KierHallSmartsDescriptor());
		descriptors.add(new KappaShapeIndicesDescriptor());
		descriptors.add(new LargestChainDescriptor());
		descriptors.add(new LargestPiSystemDescriptor());
		descriptors.add(new RuleOfFiveDescriptor());
//		descriptors.add(new LongestAliphaticChainDescriptor());
		descriptors.add(new PetitjeanNumberDescriptor());
		descriptors.add(new PetitjeanShapeIndexDescriptor());
		descriptors.add(new RotatableBondsCountDescriptor());
		descriptors.add(new TPSADescriptor());
		descriptors.add(new VAdjMaDescriptor());
		descriptors.add(new WeightDescriptor());
		descriptors.add(new WeightedPathDescriptor());
		descriptors.add(new WienerNumbersDescriptor());
		descriptors.add(new XLogPDescriptor());
		descriptors.add(new ZagrebIndexDescriptor());
		
		return descriptors;
	}
	/**
	 * getSelectedDescriptors() => either get selected descriptors or get all descriptors 
	 * @param mole
	 * @throws Exception
	 * @author xuan
	 */
	public String generateCDKMolecularDescriptors(IAtomContainer mole, ArrayList<IDescriptor> descriptors) throws Exception {
		// IDescriptorResult res = alogp.calculate(mole).getValue();
		// temp = temp + ","+ alogp.calculate(mole).getValue().toString();
		
		String temp = "";
		for(IDescriptor desc : descriptors) {
			IDescriptorResult res = ((IMolecularDescriptor) desc).calculate(mole).getValue();
			
			temp = temp + "," + res.toString();
		}

		
		return temp;
		
	}
	
	
	/**
	 * get all descriptor name plus the LogS (class attribute)
	 * getSelectedDescriptors() => either get selected descriptors or get all descriptors 
	 * @return
	 * @throws CDKException
	 */
	public String[] getDescriptorNames() throws CDKException {
		
		ArrayList<IDescriptor> descriptors = getSelectedDescriptors();
		String attribute_name = "";
		for (IDescriptor desc : descriptors) {
			String[] desr_name = desc.getDescriptorNames();

			String str = String.join(",", desr_name);
			attribute_name = attribute_name + str + ",";	
		}
		attribute_name = attribute_name + "LogS";
		String[] attribute_name_array = attribute_name.split(",");
		
		return attribute_name_array;
	}
	
	
	/**
	 * get all descriptor name plus the LogS (class attribute)
	 * @return
	 * @throws CDKException
	 */
	public ArrayList<Attribute> getAttributeNames(String[] descriptorName) throws CDKException {
		ArrayList<Attribute> attributeNames = new ArrayList<Attribute>();
		for(int i = 0; i < descriptorName.length; i++) {
			Attribute tmp = new Attribute(descriptorName[i]);
			attributeNames.add(tmp);
		}
		return attributeNames;
		
	}
	
	/**
	 * 
	 * 
	 * @param mole
	 * @param logS
	 * @return
	 * @throws CDKException
	 */
	public double[] getDescriptorValues(IAtomContainer mole, double logS, ArrayList<IDescriptor> descriptors) {
		
		double[] empty = new double[0];
		try {
			String values = generateCDKMolecularDescriptors(mole,descriptors);
			String[] tmp = values.split(",");
			String[] valuesArray = Arrays.copyOfRange(tmp, 1, tmp.length);
			double[] descriptorValues = new double[valuesArray.length+1];
			for(int i = 0; i < valuesArray.length; i++) {
				descriptorValues[i] = Double.valueOf(valuesArray[i]);
			}
			
			descriptorValues[valuesArray.length] = logS;
			
			return descriptorValues;
			
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		return empty;
		
		
	}
	
	
	
	/**
	 * convert csv to instances
	 * row[0] = compound name; row[1] = smiles (from chemspider); row[2] = logS data;
	 * @param csvfile
	 * @return
	 * @throws IOException
	 * @throws CDKException 
	 */
	public Instances getDataset(String file) throws IOException, CDKException {
		
		
		
		final CSVParser parser = new CSVParserBuilder().withSeparator('\t').withIgnoreQuotations(true).build();
		final CSVReaderBuilder builder = new CSVReaderBuilder(new FileReader(file));
		final CSVReader csvReader = builder.withCSVParser(parser).build();
		ArrayList<IDescriptor> descriptors = getSelectedDescriptors();
		String[] descriptorName = getDescriptorNames();
		ArrayList<Attribute> attributeNames = getAttributeNames(descriptorName);
		Instances dataRaw = new Instances("TestInstances", attributeNames , 0);

		csvReader.readNext(); // skip the header
		String[] tmp;
		while ((tmp = csvReader.readNext()) != null) {
			String compound_smiles = tmp[1];
			Double logSexp = Double.valueOf(tmp[2]);
			IAtomContainer mole = CheminformaticUtility.parseSmilesToContainer(compound_smiles);
			double[] descValue = getDescriptorValues(mole, logSexp,descriptors);
			dataRaw.add(new DenseInstance(1.0, descValue));
			
		}
		
		
		return dataRaw;
		
	}
	
	/**
	 * convert csv to instances
	 * row[0] = compound name; row[1] = smiles (from chemspider); row[2] = logS data;
	 * @param csvfile
	 * @return
	 * @throws IOException
	 * @throws CDKException 
	 */
	public Instances getDatasetFromSDF(String file) throws IOException, CDKException {
		
		ArrayList<IDescriptor> descriptors = getSelectedDescriptors();
		File sdfFile = new File(file);
		IteratingSDFReader reader = new IteratingSDFReader(
				new FileInputStream(sdfFile), DefaultChemObjectBuilder.getInstance());
		
		String[] descriptorName = getDescriptorNames();
		ArrayList<Attribute> attributeNames = getAttributeNames(descriptorName);
		Instances dataRaw = new Instances("TestInstances", attributeNames , 0);
		
		while (reader.hasNext()) {
			IAtomContainer molecule = reader.next();
			String logs = molecule.getProperty("logS");
			double[] descValue = getDescriptorValues(molecule, Double.valueOf(logs),descriptors);
			dataRaw.add(new DenseInstance(1.0, descValue));
		}
	
		reader.close();
		return dataRaw;
		
	}
	
	/**
	 * 
	 * @param mole
	 * @return
	 * @throws Exception 
	 * @throws FileNotFoundException 
	 */
	public double predictLogS(IAtomContainer mole) throws FileNotFoundException, Exception {
		
		String[] descriptorName = getDescriptorNames();
		ArrayList<Attribute> attributeNames = getAttributeNames(descriptorName);
		Instances dataRaw = new Instances("TestInstances", attributeNames , 0);
		ArrayList<IDescriptor> descriptors = getSelectedDescriptors();
		
		double[] descValue = getDescriptorValues(mole, 0.0, descriptors);
		dataRaw.add(new DenseInstance(1.0, descValue));
		
		RandomForest classifier = (RandomForest) SerializationHelper.read(
				new FileInputStream(String.format("%s/model/%s.model", current_dir,modelName)));
		dataRaw.setClass(dataRaw.attribute(dataRaw.numAttributes()-1));
		double result = classifier.classifyInstance(dataRaw.get(0));
		
		return result;
	}
	
	/**
	 * 
	 * @param mole
	 * @return
	 * @throws Exception 
	 * @throws FileNotFoundException 
	 */
	public double predictLogS(IAtomContainer mole, ArrayList<Attribute> attributeNames,
			ArrayList<IDescriptor> descriptors, RandomForest classifier ) throws FileNotFoundException, Exception {
		
		
		Instances dataRaw = new Instances("TestInstances", attributeNames , 0);
		double[] descValue = getDescriptorValues(mole, 0.0, descriptors);
		dataRaw.add(new DenseInstance(1.0, descValue));
		dataRaw.setClass(dataRaw.attribute(dataRaw.numAttributes()-1));
		double result = classifier.classifyInstance(dataRaw.get(0));
		
		return result;
	}
	
	/**
	 * 
	 * @param mole
	 * @return
	 * @throws FileNotFoundException
	 * @throws Exception
	 */
	public double[] predictLogS(IAtomContainerSet moleset) throws FileNotFoundException, Exception {
		
		ArrayList<IDescriptor> descriptors = getSelectedDescriptors();
		String[] descriptorName = getDescriptorNames();
		ArrayList<Attribute> attributeNames = getAttributeNames(descriptorName);
		Instances dataRaw = new Instances("TestInstances", attributeNames , 0);
		for(int i = 0; i < moleset.getAtomContainerCount(); i++) {
			double[] descValue = getDescriptorValues(moleset.getAtomContainer(i), 0.0, descriptors);
			dataRaw.add(new DenseInstance(1.0, descValue));
		
		}
		dataRaw.setClass(dataRaw.attribute(dataRaw.numAttributes()-1));
		RandomForest classifier = (RandomForest) SerializationHelper.read(
				new FileInputStream(String.format("%s/model/%s.model", current_dir,modelName)));
		double[] result = new double[moleset.getAtomContainerCount()];
		for(int i = 0; i < dataRaw.size(); i++) {
			double tmp = classifier.classifyInstance(dataRaw.get(i));
			System.out.println(tmp);
			result[i] = tmp;
		}
		
		return result;
	}
	
	/**
	 * sampling dataset for upper logS and lower logS (phytohub, etc.)
	 * @param csvfile
	 * @throws Exception 
	 */
	public void sampleDataset(String sampleFile) throws Exception {
		String csvfile = String.format("%s/generatedfolder/%s", current_dir,sampleFile);
		final CSVParser parser = new CSVParserBuilder().withSeparator(',').withIgnoreQuotations(false).build();
		final CSVReaderBuilder builder = new CSVReaderBuilder(new FileReader(csvfile));
		final CSVReader csvReader = builder.withCSVParser(parser).build();
		IChemObjectBuilder chembuilder = DefaultChemObjectBuilder.getInstance();
		IAtomContainerSet moleset = chembuilder.newInstance(IAtomContainerSet.class);
		
		FileWriter fw = new FileWriter(String.format("%s/generatedfolder/%s", current_dir,"PhytoHubCompoundsInSDF.sdf"), true);
        SDFWriter sdfwriter = new SDFWriter(fw);
        
        
        
		csvReader.readNext(); // skip the header
//		System.out.println(Arrays.toString(csvReader.readNext()));
		String[] tmp;
		while ((tmp = csvReader.readNext()) != null) {
			String phytohubid = tmp[0];
			String compound_inchi = tmp[1];
//			System.out.println(compound_inchi);
			IAtomContainer mole = CheminformaticUtility.getAtomContainerFromInChI(compound_inchi);
		
			mole.setProperty("PhytohubID", phytohubid);
			sdfwriter.write(mole); // save the sdf to view the sugars
			moleset.addAtomContainer(mole);
			
		}
		sdfwriter.close();
		
		
		
		String[] descriptorName = getDescriptorNames();
		ArrayList<Attribute> attributeNames = getAttributeNames(descriptorName);
		ArrayList<IDescriptor> descriptors = getSelectedDescriptors();
		RandomForest classifier = (RandomForest) SerializationHelper.read(
				new FileInputStream(String.format("%s/model/%s.model", current_dir,modelName)));
		double[] result = new double[moleset.getAtomContainerCount()];
		for (int i = 0; i < moleset.getAtomContainerCount(); i++) {
			try {
				String phytohubid = moleset.getAtomContainer(i).getProperty("PhytohubID");
				System.out.println(phytohubid);
				double tmpresult = predictLogS(moleset.getAtomContainer(i),attributeNames,descriptors,classifier);
				result[i] = tmpresult;
				moleset.getAtomContainer(i).removeAllElements();
				
			}catch(Exception e) {
				continue;
			}		
			// System.gc();
		}
		final CSVWriter writer = new CSVWriter(new FileWriter(String.format("%s/generatedfolder/%s", current_dir,"PhytoHubSampledLogS.csv")));
		List<String[]> resultWritable = new ArrayList<String[]>();
		for (int i = 0; i < result.length; i++) {
			resultWritable.add(new String[] {String.valueOf(result[i])});
		}
		
		writer.writeAll(resultWritable);
		writer.close();
		
	}
	
	
	/**
	 * input calculated logS and return is reasonable
	 * upper_logS and lower_logS are obtained from sampling phytohub data.
	 * removed the outler both for upper and lower
	 * @param mole
	 * @return
	 */
	public boolean isReasonableLogS(double logS) {
		double upper_logS = -1.0;
		double lower_logS = -6.0;
		if(logS >= lower_logS && logS <= upper_logS) {
			return true;
		}else{
			return false;
		}
	}
	
	/**
	 * 
	 * @throws Exception
	 */
	public void buildLogSModel() throws Exception {
		CreateSolubilityModel dewf = new CreateSolubilityModel();
		
//		Instances dataset = dewf.getDataset(String.format("%s/generatedfolder/%s", current_dir, 
//				"LogSData.tsv"));
		Instances dataset = dewf.getDatasetFromSDF(String.format("%s/generatedfolder/%s", current_dir, 
				"training_set_everything.sdf"));
		ArffSaver arffsave = new ArffSaver();
		File outputFile = new File(String.format("%s/generatedfolder/%s.arff", current_dir,"LogSTraining_set_everything"));
		arffsave.setInstances(dataset);
		arffsave.setFile(outputFile);
		arffsave.writeBatch();
	}
}


// paper used following descriptors
//org.openscience.cdk.qsar.descriptors.molecular.ALOGPDescriptor@16e7dcfd
//org.openscience.cdk.qsar.descriptors.molecular.BCUTDescriptor@7219ec67
//org.openscience.cdk.qsar.descriptors.molecular.FragmentComplexityDescriptor@61dd025
//org.openscience.cdk.qsar.descriptors.molecular.APolDescriptor@3d121db3
//org.openscience.cdk.qsar.descriptors.molecular.AromaticAtomsCountDescriptor@3b07a0d6
//org.openscience.cdk.qsar.descriptors.molecular.AromaticBondsCountDescriptor@11a9e7c8
//org.openscience.cdk.qsar.descriptors.molecular.AtomCountDescriptor@3901d134
//org.openscience.cdk.qsar.descriptors.molecular.AutocorrelationDescriptorCharge@14d3bc22
//org.openscience.cdk.qsar.descriptors.molecular.AutocorrelationDescriptorMass@12d4bf7e
//org.openscience.cdk.qsar.descriptors.molecular.AutocorrelationDescriptorPolarizability@4c1d9d4b
//org.openscience.cdk.qsar.descriptors.molecular.BondCountDescriptor@45018215
//org.openscience.cdk.qsar.descriptors.molecular.BPolDescriptor@65d6b83b
//org.openscience.cdk.qsar.descriptors.molecular.CarbonTypesDescriptor@d706f19
//org.openscience.cdk.qsar.descriptors.molecular.ChiChainDescriptor@30b7c004
//org.openscience.cdk.qsar.descriptors.molecular.ChiClusterDescriptor@79efed2d
//org.openscience.cdk.qsar.descriptors.molecular.ChiPathDescriptor@27ae2fd0
//org.openscience.cdk.qsar.descriptors.molecular.ChiPathClusterDescriptor@2928854b
//org.openscience.cdk.qsar.descriptors.molecular.EccentricConnectivityIndexDescriptor@2f177a4b
//org.openscience.cdk.qsar.descriptors.molecular.HBondDonorCountDescriptor@77167fb7
//org.openscience.cdk.qsar.descriptors.molecular.HBondAcceptorCountDescriptor@15b204a1
//org.openscience.cdk.qsar.descriptors.molecular.KierHallSmartsDescriptor@6973bf95 // smarts descriptor => fragment based descriptor
//org.openscience.cdk.qsar.descriptors.molecular.KappaShapeIndicesDescriptor@6ce139a4
//org.openscience.cdk.qsar.descriptors.molecular.LargestChainDescriptor@2ddc8ecb
//org.openscience.cdk.qsar.descriptors.molecular.LargestPiSystemDescriptor@229d10bd
//org.openscience.cdk.qsar.descriptors.molecular.RuleOfFiveDescriptor@25d250c6
//org.openscience.cdk.qsar.descriptors.molecular.LongestAliphaticChainDescriptor@33afa13b
//org.openscience.cdk.qsar.descriptors.molecular.PetitjeanNumberDescriptor@18ce0030
//org.openscience.cdk.qsar.descriptors.molecular.PetitjeanShapeIndexDescriptor@4445629 // only extract topoShape
//org.openscience.cdk.qsar.descriptors.molecular.RotatableBondsCountDescriptor@45b9a632
//org.openscience.cdk.qsar.descriptors.molecular.TPSADescriptor@6b26e945
//org.openscience.cdk.qsar.descriptors.molecular.VAdjMaDescriptor@54c562f7
//org.openscience.cdk.qsar.descriptors.molecular.WeightDescriptor@318ba8c8
//org.openscience.cdk.qsar.descriptors.molecular.WeightedPathDescriptor@6dbb137d
//org.openscience.cdk.qsar.descriptors.molecular.WienerNumbersDescriptor@43301423
//org.openscience.cdk.qsar.descriptors.molecular.XLogPDescriptor@2f112965
//org.openscience.cdk.qsar.descriptors.molecular.ZagrebIndexDescriptor@1a04f701

//ALOGPDescriptor alogp = new ALOGPDescriptor();
//BCUTDescriptor bcut   = new BCUTDescriptor();
//FragmentComplexityDescriptor fcd = new FragmentComplexityDescriptor();
//APolDescriptor apol = new APolDescriptor();
//AromaticAtomsCountDescriptor armc = new AromaticAtomsCountDescriptor();
//AromaticBondsCountDescriptor abc = new AromaticBondsCountDescriptor();
//AtomCountDescriptor atomc = new AtomCountDescriptor();
//AutocorrelationDescriptorCharge acdc = new AutocorrelationDescriptorCharge();














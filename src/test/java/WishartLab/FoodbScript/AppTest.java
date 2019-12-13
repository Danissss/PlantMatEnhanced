package WishartLab.FoodbScript;

import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IChemObjectBuilder;

import WishartLab.PlantGlycosider.CheminformaticUtility;
import WishartLab.PlantGlycosider.GlycosiderUtility;
import WishartLab.PlantGlycosider.PostValidation;
import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

/**
 * Unit test for simple App.
 */
public class AppTest 
    extends TestCase
{
    /**
     * Create the test case
     *
     * @param testName name of the test case
     */
    public AppTest( String testName )
    {
        super( testName );
    }

    /**
     * @return the suite of tests being tested
     */
    public static Test suite()
    {
        return new TestSuite( AppTest.class );
    }

    /**
     * Rigourous Test :-)
     */
    public void testApp()
    {
        assertTrue( true );
    }
    
    	
    public void testCombineContainers() throws CDKException {
    		IChemObjectBuilder builder = DefaultChemObjectBuilder.getInstance();
		IAtomContainer mole = CheminformaticUtility.parseSmilesToContainer("CC1=C(O)C2=C3C(=C1)C1=CC=C(C=C1OC3(OC1=C2C=CC(O)=C1)C1=CC=CC=C1)C1=CC2=C(O1)C=CC=C2");
		Integer[] mole_index = new Integer[] {1};
		double [] charge_table = new double[mole_index.length];
		
		for(int i=0; i < mole_index.length;i++) {
			double charge = (double) i+2;
			mole.getAtom(mole_index[i]).setCharge(charge);
			charge_table[i] = charge;				
		}
		
		
		IAtomContainerSet sugar_mole_set = builder.newInstance(IAtomContainerSet.class);
		sugar_mole_set.addAtomContainer(CheminformaticUtility.TransformSugarSmilesToContainer("OC1COC(O)C(O)C1O",5));
		
    }
    
    
    public void testPostValidation() throws CDKException {
    		IChemObjectBuilder builder = DefaultChemObjectBuilder.getInstance();
		IAtomContainerSet mole = builder.newInstance(IAtomContainerSet.class);
		mole.addAtomContainer(CheminformaticUtility.parseSmilesToContainer("CC1=C(O)C2=C3C(=C1)C1=CC=C(C=C1OC3(OC1=C2C=CC(O)=C1)C1=CC=CC=C1)C1=CC2=C(O1)C=CC=C2"));
		PostValidation pv = new PostValidation();
//		pv.validateCompound(mole);
    }
    
}

package org.openscience.cdk.qsar.descriptors.molecular;

import org.junit.Assert;
import org.junit.Test;
import org.openscience.cdk.Molecule;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IMolecule;
import org.openscience.cdk.io.IChemObjectReader.Mode;
import org.openscience.cdk.io.MDLReader;
import org.openscience.cdk.qsar.IMolecularDescriptor;
import org.openscience.cdk.qsar.result.DoubleResult;
import org.openscience.cdk.tools.GenerateFragments;

import java.io.InputStream;
import java.util.List;

/**
 * @author      chhoppe from EUROSCREEN
 * @cdk.module  test-qsarmolecular
 */
public class FragmentComplexityDescriptorTest extends MolecularDescriptorTest {
	
	public  FragmentComplexityDescriptorTest() {}

	   
	@Test public void test1FragmentComplexityDescriptor() throws Exception {
		IMolecularDescriptor descriptor = new FragmentComplexityDescriptor();
		String filename = "data/mdl/murckoTest1.mol";
    	//System.out.println("\nFragmentComplexityTest: " + filename);
    	InputStream ins = this.getClass().getClassLoader().getResourceAsStream(filename);
    	GenerateFragments gf=new GenerateFragments();
    	double Complexity=0;
        	MDLReader reader = new MDLReader(ins, Mode.STRICT);
        	Molecule mol = (Molecule)reader.read(new Molecule());
        	gf.generateMurckoFragments(mol,false,false,4);
        	List setOfFragments=gf.getMurckoFrameworks();
        	for (int i=0;i<setOfFragments.size();i++){
        		addExplicitHydrogens((IMolecule) setOfFragments.get(i));
        		Complexity=((DoubleResult)descriptor.calculate((IAtomContainer) setOfFragments.get(i)).getValue()).doubleValue();
        		//System.out.println("Complexity:"+Complexity);
        	}
    	 Assert.assertEquals(659.00,Complexity , 0.01);
	}
	
	@Test
    public void test2FragmentComplexityDescriptor() throws Exception {
		IMolecularDescriptor descriptor = new FragmentComplexityDescriptor();
		String filename = "data/mdl/murckoTest10.mol";
    	//System.out.println("\nFragmentComplexityTest: " + filename);
    	InputStream ins = this.getClass().getClassLoader().getResourceAsStream(filename);
    	GenerateFragments gf=new GenerateFragments();
    	double Complexity=0;
    	MDLReader reader = new MDLReader(ins, Mode.STRICT);
    	Molecule mol = (Molecule)reader.read(new Molecule());
    	gf.generateMurckoFragments(mol,false,false,4);
    	List setOfFragments=gf.getMurckoFrameworks();
    	for (int i=0;i<setOfFragments.size();i++){
    		addExplicitHydrogens((IMolecule) setOfFragments.get(i));
    		Complexity=((DoubleResult)descriptor.calculate((IAtomContainer) setOfFragments.get(i)).getValue()).doubleValue();
    		//System.out.println("Complexity:"+Complexity);
    	}
    	Assert.assertEquals(544.01,Complexity , 0.01);
	}
	
	
}

 

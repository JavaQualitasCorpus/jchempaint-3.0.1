/* 
 * Copyright (C) 2004-2007  The Chemistry Development Kit (CDK) project
 * 
 * Contact: cdk-devel@lists.sourceforge.net
 * 
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA. 
 */
package org.openscience.cdk.tools;

import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.List;

import org.junit.Assert;
import org.junit.Test;
import org.openscience.cdk.CDKTestCase;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IMolecule;
import org.openscience.cdk.interfaces.IRingSet;
import org.openscience.cdk.io.ISimpleChemObjectReader;
import org.openscience.cdk.io.MDLV2000Reader;
import org.openscience.cdk.nonotify.NNMolecule;
import org.openscience.cdk.ringsearch.SSSRFinder;
import org.openscience.cdk.smiles.SmilesParser;

/**
 * TestSuite that runs all QSAR tests.
 *
 * @cdk.module test-extra
 */
public class GenerateFragmentsTest extends CDKTestCase{
	
	private GenerateFragments gf = new GenerateFragments();

	@Test public void testGenerateMurckoFragments1() throws Exception {
	    	String filename = "data/mdl/murckoTest1.mol";
	    	//logger.debug("\nMurckoTesting: " + filename);
	    	InputStream ins = this.getClass().getClassLoader().getResourceAsStream(filename);
	    	ISimpleChemObjectReader reader = new MDLV2000Reader(ins);
	    	IMolecule mol = (IMolecule)reader.read(new NNMolecule());
	    	gf.generateMurckoFragments(mol,true,true,4);
	    	//logger.debug("Murcko Fragments generated");
	    	String[] smiles=gf.getMurckoFrameworksAsSmileArray();
	    	for (int i =0;i<smiles.length;i++){
	    		//logger.debug("MF"+i+" :"+smiles[i]);
	    	}
	    	Assert.assertEquals(1,smiles.length);
	    	Assert.assertEquals("C1=CC=C(C=C1)CCC2=CC=CC=C2",smiles[0]);
	}

	@Test public void testGenerateMurckoFragments2() throws Exception {
    	String filename = "data/mdl/murckoTest2.mol";
    	//logger.debug("\nMurckoTesting: " + filename);
    	InputStream ins = this.getClass().getClassLoader().getResourceAsStream(filename);
    	ISimpleChemObjectReader reader = new MDLV2000Reader(ins);
    	IMolecule mol = (IMolecule)reader.read(new NNMolecule());
    	gf.generateMurckoFragments(mol,true,true,4);
    	//logger.debug("Murcko Fragments generated");
    	String[] smiles=gf.getMurckoFrameworksAsSmileArray();
    	for (int i =0;i<smiles.length;i++){
    		//logger.debug("MF"+i+" :"+smiles[i]);
    	}
    	Assert.assertEquals(1,smiles.length);
    	Assert.assertEquals("C1CCC(C1)CCC2CC3=CC=CC=C3(C2)",smiles[0]);
	}

	@Test public void testGenerateMurckoFragments3() throws Exception {
    	String filename = "data/mdl/murckoTest3.mol";
    	//logger.debug("\nMurckoTesting: " + filename);
    	InputStream ins = this.getClass().getClassLoader().getResourceAsStream(filename);
    	ISimpleChemObjectReader reader = new MDLV2000Reader(ins);
    	IMolecule mol = (IMolecule)reader.read(new NNMolecule());
    	gf.generateMurckoFragments(mol,true,true,4);
    	//logger.debug("Murcko Fragments generated");
    	String[] smiles=gf.getMurckoFrameworksAsSmileArray();
    	Assert.assertEquals(1,smiles.length);
    	Assert.assertEquals("C=1C=CC=2CC=CC=2(C=1)", smiles[0]);
	}

	@Test public void testGenerateMurckoFragments4() throws Exception {
    	String filename = "data/mdl/murckoTest4.mol";
    	//logger.debug("\nMurckoTesting: " + filename);
    	InputStream ins = this.getClass().getClassLoader().getResourceAsStream(filename);
    	ISimpleChemObjectReader reader = new MDLV2000Reader(ins);
    	IMolecule mol = (IMolecule)reader.read(new NNMolecule());
    	gf.generateMurckoFragments(mol,true,true,4);
    	//logger.debug("Murcko Fragments generated");
    	String[] smiles=gf.getMurckoFrameworksAsSmileArray();
    	Assert.assertEquals(3,smiles.length);
    	boolean found=false;
    	for (int i =0;i<smiles.length;i++){
//    		System.out.println("MF"+i+" :"+smiles[i]);
    		if (smiles[i].equals("C1=CC=C(C=C1)CCC3CCC(CCC2=CC=CC=C2)C3")){
    			found=true;
    		}
    	}
    	Assert.assertTrue(found);
    	//Assert.assertEquals("c1ccc(cc1)CCC2CCC(C2)C4C4(c3ccccc3)",smiles[2]);

    	/*String[] rings=gf.getRingFragmentsAsSmileArray();
        	for (int i =0;i<rings.length;i++){
        		System.out.println("RF"+i+" :"+smiles[i]);
        	}*/
    	String[] linker=gf.getLinkerFragmentsAsSmileArray();
    	for (int i =0;i<linker.length;i++){
    		//logger.debug("LF"+i+" :"+linker[i]);
    	}
	}
	
	@Test public void testGenerateMurckoFragments5() throws Exception {
    	String filename = "data/mdl/murckoTest5.mol";
    	//logger.debug("\nMurckoTesting: " + filename);
    	InputStream ins = this.getClass().getClassLoader().getResourceAsStream(filename);
    	ISimpleChemObjectReader reader = new MDLV2000Reader(ins);
    	IMolecule mol = (IMolecule)reader.read(new NNMolecule());
    	gf.generateMurckoFragments(mol,true,true,4);
    	//logger.debug("Murcko Fragments generated");
    	String[] smiles=gf.getMurckoFrameworksAsSmileArray();
    	for (int i =0;i<smiles.length;i++){
    		//logger.debug("MF"+i+" :"+smiles[i]);
    	}
    	Assert.assertEquals(0,smiles.length);
	}
	
	@Test public void testGenerateMurckoFragments6() throws Exception {
    	String filename = "data/mdl/murckoTest6.mol";
    	//logger.debug("\nMurckoTesting: " + filename);
    	InputStream ins = this.getClass().getClassLoader().getResourceAsStream(filename);
    	ISimpleChemObjectReader reader = new MDLV2000Reader(ins);
    	IMolecule mol = (IMolecule)reader.read(new NNMolecule());
    	gf.generateMurckoFragments(mol,true,true,4);
    	//logger.debug("Murcko Fragments generated");
    	String[] smiles=gf.getMurckoFrameworksAsSmileArray();
    	Assert.assertEquals(1,smiles.length);
    	Assert.assertEquals("NC(CC1=CC=CC=C1)C(=O)C2=CC=CC=C2", smiles[0]);
	}
	
	@Test public void testGenerateMurckoFragments7() throws Exception {
    	//logger.debug("\nMurckoTesting 7");
    	SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        String smile="OC(CC[n+]1c(Nc2ccccc2)scc1c3ccccc3)(P(=O)([O-])[O-])P(=O)([O-])[O-]";//ZINK5
        IMolecule mol = sp.parseSmiles(smile); 
      
    	gf.generateMurckoFragments(mol,true,true,4);
    	//logger.debug("Murcko Fragments generated");
    	String[] smiles=gf.getMurckoFrameworksAsSmileArray();
    	Assert.assertEquals(3,smiles.length);
    	for (int i =0;i<smiles.length;i++){
    		//logger.debug("MF"+i+" :"+smiles[i]);
    	}
    	//Assert.assertEquals("NC3(Cc1ccccc1)(C3(=O)(c2ccccc2))",smiles[0]);
	}
	
	@Test public void testGenerateMurckoFragments8() throws Exception {
    	//logger.debug("\nMurckoTesting 8");
    	SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        String smile="c2ccc(Cc1ccccc1)cc2";
        IMolecule mol = sp.parseSmiles(smile); 
      
    	gf.generateMurckoFragments(mol,true,true,4);
    	//logger.debug("Murcko Fragments generated");
    	String[] smiles=gf.getMurckoFrameworksAsSmileArray();
    	Assert.assertEquals(1,smiles.length);
    	for (int i =0;i<smiles.length;i++){
    		//logger.debug("MF"+i+" :"+smiles[i]);
    	}
    	Assert.assertEquals("c1=cc=c(c=c1)Cc2=cc=cc=c2",smiles[0]);
	}
	
	@Test public void testGenerateMurckoFragments9() throws Exception {
    	//logger.debug("\nMurckoTesting 9");
    	SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        String smile="c2ccc(c1ccccc1)cc2";
        IMolecule mol = sp.parseSmiles(smile); 
      
    	gf.generateMurckoFragments(mol,true,true,4);
    	//logger.debug("Murcko Fragments generated");
    	String[] smiles=gf.getMurckoFrameworksAsSmileArray();
    	Assert.assertEquals(1,smiles.length);
    	for (int i =0;i<smiles.length;i++){
    		//logger.debug("MF"+i+" :"+smiles[i]);
    	}
    	Assert.assertEquals("c1=cc=c(c=c1)c2=cc=cc=c2",smiles[0]);
	}
	
	// The next test creates an invalid SMILES for the fragment
	/*@Test public void testGenerateMurckoFragments10() throws Exception {
    	//logger.debug("\nMurckoTesting 10");
    	SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        String smile="Cc1nn(C)cc1[C@H]2[C@H](C(=O)N)C(=O)C[C@@](C)(O)[C@@H]2C(=O)N";//ZINK19
                      
        IMolecule mol = sp.parseSmiles(smile); 
        
       	GenerateFragments gf=new GenerateFragments();
    	try {
        	
        	gf.generateMurckoFragments(mol,true,true);
        	System.out.println("Murcko Fragments generated");
        	String[] smiles=gf.getMurckoFrameworksAsSmileArray();
        	for (int i =0;i<smiles.length;i++){
        		//logger.debug("MF"+i+" :"+smiles[i]);
        	}
        	//Assert.assertEquals("NC3(Cc1ccccc1)(C3(=O)(c2ccccc2))",smiles[0]);
        	//Assert.assertEquals(1,smiles.length);
        }catch (Exception e){
        	System.out.println("Error in testGenerateMurckoFragments6:");
        	e.printStackTrace();
        }
	}*/
	
	@Test public void testGenerateMurckoFragments11() throws Exception {
    	String filename = "data/mdl/murckoTest7.mol";
    	//logger.debug("\nMurckoTesting: " + filename);
    	InputStream ins = this.getClass().getClassLoader().getResourceAsStream(filename);
    	ISimpleChemObjectReader reader = new MDLV2000Reader(ins);
    	IMolecule mol = (IMolecule)reader.read(new NNMolecule());
    	gf.generateMurckoFragments(mol,true,true,4);
    	//logger.debug("Murcko Fragments generated");
    	String[] smiles=gf.getMurckoFrameworksAsSmileArray();
    	Assert.assertEquals(3,smiles.length);
    	boolean found=false;
    	for (int i =0;i<smiles.length;i++){
//    		System.out.println("MF"+i+" :"+smiles[i]);
    		if (smiles[i].equals("C1=CC=C(C=C1)C2CCCC(C2)C3=CC=CC=C3")){
    			found=true;
    		}
    	}
    	Assert.assertTrue(found);
	}
	
	@Test public void testGenerateMurckoFragments12() throws Exception {
    	String filename = "data/mdl/murckoTest8.mol";
    	//logger.debug("\nMurckoTesting: " + filename);
    	InputStream ins = this.getClass().getClassLoader().getResourceAsStream(filename);
    	ISimpleChemObjectReader reader = new MDLV2000Reader(ins);
    	IMolecule mol = (IMolecule)reader.read(new NNMolecule());
    	gf.generateMurckoFragments(mol,true,true,4);
    	//logger.debug("Murcko Fragments generated");
    	String[] smiles=gf.getMurckoFrameworksAsSmileArray();
    	Assert.assertEquals(3,smiles.length);
    	boolean found=false;
    	for (int i =0;i<smiles.length;i++){
    		if (smiles[i].equals("C1=CC=C(C=C1)C2CC(C2)C3=CC=CC=C3")){
    			found=true;
    		}
//    		System.out.println("MF"+i+" :"+smiles[i]);
    	}
    	Assert.assertTrue(found);
	}

	//without add explicit hydrogen thetest fails due to problems with smile generator 
	@Test public void testGenerateMurckoFragments13() throws Exception {
		//logger.debug("\nMurckoTesting 13");
		SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
		String smile="Oc1cc2ccccn2c1C(=O)OCCN3CCCCC3";//MDDR 31 
                  
		IMolecule mol = sp.parseSmiles(smile); 
		addExplicitHydrogens(mol);
    
		gf.generateMurckoFragments(mol,false,false,4);
		String[] smiles=gf.getMurckoFrameworksAsSmileArray();
		Assert.assertEquals(1,smiles.length);
		Assert.assertEquals("C1CCN(CC1)CCOCc2ccc3ccccn23", smiles[0]);
	}
	
	//same as test 13
	//add test for linkers
	@Test public void testGenerateMurckoFragments14() throws Exception {
		//logger.debug("\nMurckoTesting 14");
		SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
		String smile="C(c1ccc(cc1)c2ccccc2)n3cnc4cccnc34";//MDDR 52 
                  
		IMolecule mol = sp.parseSmiles(smile); 
   
		addExplicitHydrogens(mol);
    
		gf.generateMurckoFragments(mol,false,false,4);
		String[] smiles=gf.getMurckoFrameworksAsSmileArray();
		String[] linkers=gf.getLinkerFragmentsAsSmileArray();
		int found=0;
		Assert.assertEquals(3,linkers.length);
		for (int i =0;i<linkers.length;i++){
			//System.out.println("linkers"+i+" :"+linkers[i]);
			if (linkers[i].equals("nCc1=cc=c(c)c=c1") || linkers[i].equals("cCn") || linkers[i].equals("cc")){
				found++;
			}
		}
		Assert.assertEquals(3,found);
		
		Assert.assertEquals(3,smiles.length);
		boolean fnd=false;
		for (int i =0;i<smiles.length;i++){
			//System.out.println("MF"+i+" :"+smiles[i]);
			//if (smiles[i].equals("c1=cc=c(c=c1)c2=cc=c(c=c2)Cn4c=nc=3c=cc=nc=34")){
			if (smiles[i].equals("c1=cc=c(c=c1)c2=cc=c(c=c2)Cn4c=nc3=cc=cn=c34")){
				fnd=true;
			}
		}
		Assert.assertTrue(fnd);
	}
	
	//check for spiro ring systems
	@Test public void testGenerateMurckoFragments15() throws Exception {
    	String filename = "data/mdl/murckoTest9.mol";
    	//logger.debug("\nMurckoTesting: " + filename);
    	InputStream ins = this.getClass().getClassLoader().getResourceAsStream(filename);
    	ISimpleChemObjectReader reader = new MDLV2000Reader(ins);
    	IMolecule mol = (IMolecule)reader.read(new NNMolecule());
    	gf.generateMurckoFragments(mol,true,true,4);
    	//logger.debug("Murcko Fragments generated");
    	String[] smiles=gf.getMurckoFrameworksAsSmileArray();
    	for (int i =0;i<smiles.length;i++){
    		//logger.debug("MF"+i+" :"+smiles[i]);
    	}
    	Assert.assertEquals("C1CCC2(CC1)(CCCC2)",smiles[0]);
    	Assert.assertEquals(1,smiles.length);
	}
	
	@Test public void testGenerateMurckoFragments16() throws Exception {
    	String filename = "data/mdl/murckoTest10.mol";
    	//logger.debug("\nMurckoTesting: " + filename);
    	InputStream ins = this.getClass().getClassLoader().getResourceAsStream(filename);
    	ISimpleChemObjectReader reader = new MDLV2000Reader(new InputStreamReader(ins));
    	IMolecule mol = (IMolecule)reader.read(new NNMolecule());
    	gf.generateMurckoFragments(mol,true,true,4);
    	//logger.debug("Murcko Fragments generated");
    	String[] smiles=gf.getMurckoFrameworksAsSmileArray();
    	Assert.assertEquals(1,smiles.length);
    	for (int i =0;i<smiles.length;i++){
//    		System.out.println("MF"+i+" :"+smiles[i]);
    	}
    	Assert.assertEquals("O=C2CC=CN2(CC1=CC=CC=C1)",smiles[0]);
	}

	@Test public void testGenerateMurckoFragments17() throws Exception {
    	String filename = "data/mdl/murckoTest11.mol";
    	//logger.debug("\nMurckoTesting: " + filename);
    	InputStream ins = this.getClass().getClassLoader().getResourceAsStream(filename);
    	ISimpleChemObjectReader reader = new MDLV2000Reader(new InputStreamReader(ins));
    	IMolecule mol = (IMolecule)reader.read(new NNMolecule());
    	//System.out.println(mol.toString());
    	//This strange, when the next line is comment out the test will fail
    	mol.toString();
    	gf.generateMurckoFragments(mol,true,true,4);
    	//logger.debug("Murcko Fragments generated");
    	String[] smiles=gf.getMurckoFrameworksAsSmileArray();
    	Assert.assertEquals(1,smiles.length);
    	for (int i =0;i<smiles.length;i++){
    		//System.out.println("MF"+i+" :"+smiles[i]);
    	}
    	Assert.assertEquals("C1CCC2=CC=CC=C2(C1)",smiles[0]);
	}

    /**
     * @cdk.bug 2729120
     */
    @Test
    public void testFourRings() throws Exception {
        String smiles = "C1(c2ccccc2)(CC(CC1)CCc1ccccc1)CC1C=CC=C1";
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IMolecule mol = sp.parseSmiles(smiles);
    	gf.generateMurckoFragments(mol,true,true,5);        
        String[] mfrags = gf.getMurckoFrameworksAsSmileArray();
        Assert.assertEquals(6, mfrags.length);

        boolean found = false;
        for (String frag : mfrags) {
            if (frag.equals("C=1C=CC(C=1)CC2CCCC2")) {
                found = true;
                break;
            }
        }
        Assert.assertTrue(found);
        found = false;
        for (String frag : mfrags) {
            if (frag.equals("C1CCC(CC1)CCC2CC(CC2)CC3C=CC=C3")) {
                found = true;
                break;
            }
        }
        Assert.assertTrue(found);

    }

    /**
     * @cdk.bug 1848591
     */
    @Test public void testFragmentWithThreeRings() throws Exception {
        String filename = "data/mdl/bugtest.mol";
        InputStream ins = this.getClass().getClassLoader().getResourceAsStream(filename);
        ISimpleChemObjectReader reader = new MDLV2000Reader(ins);
        IMolecule mol = (IMolecule)reader.read(new NNMolecule());
        gf.generateMurckoFragments(mol,true,true,4);
        List frameworks = gf.getMurckoFrameworks();
        Assert.assertEquals(1,frameworks.size());
        SSSRFinder sssr = new SSSRFinder((IAtomContainer)frameworks.get(0));
        IRingSet ringSet = sssr.findSSSR();
        Assert.assertEquals(3, ringSet.getAtomContainerCount());
    }
}
	            

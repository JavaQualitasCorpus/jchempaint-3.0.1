/* $Revision$ $Author$ $Date$
 * 
 *  Copyright (C) 2007  Miguel Rojasch <miguelrojasch@users.sf.net>
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
package org.openscience.cdk.tools.manipulator;

import java.io.IOException;
import java.util.List;

import org.junit.Assert;
import org.junit.Test;
import org.openscience.cdk.CDKTestCase;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.config.IsotopeFactory;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.formula.MolecularFormula;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.interfaces.IElement;
import org.openscience.cdk.interfaces.IIsotope;
import org.openscience.cdk.interfaces.IMolecularFormula;
import org.openscience.cdk.nonotify.NoNotificationChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;

/**
 * Checks the functionality of the MolecularFormulaManipulator.
 *
 * @cdk.module test-formula
 */
public class MolecularFormulaManipulatorTest extends CDKTestCase {

	private final static  IChemObjectBuilder builder = NoNotificationChemObjectBuilder.getInstance();
	private IsotopeFactory ifac;

	/**
	 *  Constructor for the MolecularFormulaManipulatorTest object.
	 *
	 */
	public MolecularFormulaManipulatorTest(){
		
		super();
		try {
			ifac = IsotopeFactory.getInstance(builder);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}


    /**
	 * A unit test suite for JUnit.
	 *
	 * @return    The test suite
	 */
    @Test 
    public void testGetAtomCount_IMolecularFormula() {
    	
        IMolecularFormula formula = new MolecularFormula();
        formula.addIsotope(builder.newIsotope("C") );
        formula.addIsotope( builder.newIsotope("H") ,3);

        Assert.assertEquals(2, formula.getIsotopeCount());
        
        Assert.assertEquals(4, MolecularFormulaManipulator.getAtomCount(formula));
    }
    

    /**
	 * A unit test suite for JUnit.
	 *
	 * @return    The test suite
	 */
    @Test 
    public void testGetElementCount_IMolecularFormula_IElement() {
        IMolecularFormula formula = new MolecularFormula();
    	IIsotope carb = builder.newIsotope("C");
    	IIsotope flu = builder.newIsotope("F");
    	IIsotope h1 = builder.newIsotope("H");
    	IIsotope h2 = builder.newIsotope("H");
    	h2.setExactMass(2.014101778);
    	formula.addIsotope( carb ,2);
    	formula.addIsotope( flu );
    	formula.addIsotope( h1 ,3);
    	formula.addIsotope( h2 ,4);
    	
        Assert.assertEquals(10, MolecularFormulaManipulator.getAtomCount(formula));
        Assert.assertEquals(4, formula.getIsotopeCount());
        Assert.assertEquals(3, formula.getIsotopeCount(h1));
        Assert.assertEquals(4, formula.getIsotopeCount(h2));
        
        Assert.assertEquals(2, MolecularFormulaManipulator.getElementCount(formula,builder.newElement(carb)));
        Assert.assertEquals(1, MolecularFormulaManipulator.getElementCount(formula,builder.newElement(flu)));
        Assert.assertEquals(7, MolecularFormulaManipulator.getElementCount(formula,builder.newElement(h1)));
    }
    /**
	 * A unit test suite for JUnit. Not null.
	 *
	 * @return    The test suite
	 */
	@Test 
    public void testGetIsotopes_IMolecularFormula_IElement() {
		IMolecularFormula formula = new MolecularFormula();
    	IIsotope carb = builder.newIsotope("C");
    	IIsotope flu = builder.newIsotope("F");
    	IIsotope h1 = builder.newIsotope("H");
    	IIsotope h2 = builder.newIsotope("H");
    	h2.setExactMass(2.014101778);
    	formula.addIsotope( carb ,1);
    	formula.addIsotope( flu );
    	formula.addIsotope( h1 ,1);
    	formula.addIsotope( h2 ,2);
    	
    	List<IIsotope> isotopes = MolecularFormulaManipulator.getIsotopes(formula, builder.newElement("H"));
		Assert.assertEquals(2,isotopes.size());
	}
	/**
	 * A unit test suite for JUnit. Not null.
	 *
	 * @return    The test suite
	 */
	@Test 
    public void testContainsElement_IMolecularFormula_IElement() {
		IMolecularFormula formula = new MolecularFormula();
    	IIsotope carb = builder.newIsotope("C");
    	IIsotope flu = builder.newIsotope("F");
    	IIsotope h1 = builder.newIsotope("H");
    	IIsotope h2 = builder.newIsotope("H");
    	h2.setExactMass(2.014101778);
    	formula.addIsotope( carb ,1);
    	formula.addIsotope( flu );
    	formula.addIsotope( h1 ,1);
    	formula.addIsotope( h2 ,2);
    	
		Assert.assertTrue(MolecularFormulaManipulator.containsElement(formula, builder.newElement("C")));
		Assert.assertTrue(MolecularFormulaManipulator.containsElement(formula, builder.newElement("H")));
		Assert.assertTrue(MolecularFormulaManipulator.containsElement(formula, builder.newElement("F")));
	}
	
    /**
	 * A unit test suite for JUnit. Not null.
	 *
	 * @return    The test suite
	 */
	@Test 
    public void testGetString_IMolecularFormula_Empty() {
		String stringMF = MolecularFormulaManipulator.getString(new MolecularFormula());
		Assert.assertNotNull(stringMF);
		Assert.assertEquals("",stringMF);
	}
	/**
	 * A unit test suite for JUnit. Not null.
	 *
	 * @return    The test suite
	 */
	@Test 
    public void testGetString_IMolecularFormula_arrayString_boolean() {
		IMolecularFormula formula = new MolecularFormula();
		formula.addIsotope(builder.newIsotope("C"), 2);
		formula.addIsotope(builder.newIsotope("H"), 2);
		Assert.assertEquals("C2H2",MolecularFormulaManipulator.getString(formula));
		
		String[] newOrder = new String[2];
		newOrder[0] = "H";
		newOrder[1] = "C";
		
		Assert.assertEquals("H2C2",MolecularFormulaManipulator.getString(formula,newOrder,true));
		
	}
	/**
	 * A unit test suite for JUnit. Not null.
	 *
	 * @return    The test suite
	 */
	@Test 
    public void testPutInOrder_arrayString_IMolecularFormula() {
		IMolecularFormula formula = new MolecularFormula();
		formula.addIsotope(builder.newIsotope("C"), 2);
		formula.addIsotope(builder.newIsotope("H"), 2);
		
		String[] newOrder = new String[2];
		newOrder[0] = "H";
		newOrder[1] = "C";
		
		List<IIsotope> list = MolecularFormulaManipulator.putInOrder(newOrder, formula);
		Assert.assertEquals("H",list.get(0).getSymbol());
		Assert.assertEquals("C",list.get(1).getSymbol());
		
		newOrder = new String[2];
		newOrder[0] = "C";
		newOrder[1] = "H";
		
		list = MolecularFormulaManipulator.putInOrder(newOrder, formula);
		Assert.assertEquals("C",list.get(0).getSymbol());
		Assert.assertEquals("H",list.get(1).getSymbol());
		
	}
    /**
	 * A unit test suite for JUnit. Not null.
	 *
	 * @return    The test suite
	 */
	@Test 
    public void testGetString__String_IMolecularFormula() 	{
		Assert.assertNotNull(MolecularFormulaManipulator.getMolecularFormula("C10H16", new MolecularFormula()));
		Assert.assertNotNull(MolecularFormulaManipulator.getMolecularFormula("C10H16", builder));
	}
	
	/**
	 * A unit test suite for JUnit.
	 *
	 * @return    The test suite
	 */
	@Test 
    public void testGetString_IMolecularFormula()	{
		IMolecularFormula mf1 = new MolecularFormula();
		mf1.addIsotope(builder.newIsotope("C"),10);
		mf1.addIsotope(builder.newIsotope("H"),16);
		
		Assert.assertEquals("C10H16", MolecularFormulaManipulator.getString(mf1));
		
		IMolecularFormula mf2 = new MolecularFormula();
		mf2.addIsotope(builder.newAtom("H"),16);
		mf2.addIsotope(builder.newAtom("C"),10);
		
		Assert.assertEquals("C10H16", MolecularFormulaManipulator.getString(mf2));
		
		Assert.assertEquals(MolecularFormulaManipulator.getString(mf2), MolecularFormulaManipulator.getString(mf1));
		
	}

    /**
     * @cdk.bug 2276507
     */
    @Test
    public void testBug2276507() {
        IMolecularFormula mf1 = new MolecularFormula();
        mf1.addIsotope(builder.newIsotope("C"), 1);
        mf1.addIsotope(builder.newIsotope("H"), 4);

        Assert.assertEquals("CH4", MolecularFormulaManipulator.getString(mf1));        
    }

    /**
     * 
     */
    @Test
    public void testGetString_IMolecularFormula_boolean() {
        IMolecularFormula mf1 = new MolecularFormula();
        mf1.addIsotope(builder.newIsotope("C"), 1);
        mf1.addIsotope(builder.newIsotope("H"), 4);

        Assert.assertEquals("C1H4", MolecularFormulaManipulator.getString(mf1,true));        
    }

    @Test
    public void testGetString_Isotopes()	{
		IMolecularFormula mf1 = new MolecularFormula();
		mf1.addIsotope(builder.newIsotope("C",12),9);
		mf1.addIsotope(builder.newIsotope("C",13),1);
		mf1.addIsotope(builder.newIsotope("H"),16);
		
		Assert.assertEquals("C10H16", MolecularFormulaManipulator.getString(mf1));
	}

	/**
	 * A unit test suite for JUnit.
	 *
	 * @return    The test suite
	 */
	@Test 
    public void testGetMolecularFormula_String_IChemObjectBuilder()	{
		IMolecularFormula molecularFormula = MolecularFormulaManipulator.getMolecularFormula("C10H16", builder);
		
		Assert.assertEquals(26, MolecularFormulaManipulator.getAtomCount(molecularFormula));
		Assert.assertEquals(2,  molecularFormula.getIsotopeCount());
		
	}
	
	
	/**
	 * A unit test suite for JUnit. 
	 *
	 * @return    The test suite
	 */
	@Test 
    public void testGetMolecularFormula_String_IMolecularFormula()	{

		IMolecularFormula mf1 = new MolecularFormula();
		mf1.addIsotope(builder.newIsotope("C"),10);
		mf1.addIsotope(builder.newIsotope("H"),16);
		
		Assert.assertEquals(26, MolecularFormulaManipulator.getAtomCount(mf1));
		Assert.assertEquals(2, mf1.getIsotopeCount());
		
		IMolecularFormula mf2 = MolecularFormulaManipulator.getMolecularFormula("C11H17",mf1);
		
		
		Assert.assertEquals(54, MolecularFormulaManipulator.getAtomCount(mf2));
		Assert.assertEquals(2, mf2.getIsotopeCount());
	}
	
	@Test 
    public void testGetMajorIsotopeMolecularFormula_String_IChemObjectBuilder() throws Exception {
		IMolecularFormula mf2 = MolecularFormulaManipulator.getMajorIsotopeMolecularFormula("C11H17", builder);
		
		Assert.assertEquals(28, MolecularFormulaManipulator.getAtomCount(mf2));
		Assert.assertEquals(2, mf2.getIsotopeCount());
		IIsotope carbon = IsotopeFactory.getInstance(builder).getMajorIsotope("C");
		IIsotope hydrogen = IsotopeFactory.getInstance(builder).getMajorIsotope("H");
		double totalMass = carbon.getExactMass()*11;
		totalMass += hydrogen.getExactMass()*17;
		Assert.assertEquals(totalMass, MolecularFormulaManipulator.getTotalExactMass(mf2), 0.0000001);
	}
	
    /**
	 * A unit test suite for JUnit.
	 *
	 * @return    The test suite
	 */
    @Test 
    public void testRemoveElement_IMolecularFormula_IElement() {
		IMolecularFormula formula = new MolecularFormula();
		formula.addIsotope(builder.newIsotope("C"),1);
		IIsotope fl = builder.newIsotope("F");
		IIsotope hy2 = builder.newIsotope("H");
		IIsotope hy1 = builder.newIsotope("H");
		hy2.setExactMass(2.014101778);
		formula.addIsotope(fl,1);
		formula.addIsotope(hy1,2);
		formula.addIsotope(hy2,1);
        
		Assert.assertEquals(4, formula.getIsotopeCount());
        
		formula = MolecularFormulaManipulator.removeElement(formula,builder.newElement("F"));

		Assert.assertEquals(3, formula.getIsotopeCount());
        Assert.assertEquals(4, MolecularFormulaManipulator.getAtomCount(formula));
        

        formula = MolecularFormulaManipulator.removeElement(formula,builder.newElement("H"));
        
        Assert.assertEquals(1, MolecularFormulaManipulator.getAtomCount(formula));
        Assert.assertEquals(1, formula.getIsotopeCount());
        
    }
	/**
     * A unit test suite for JUnit. Test total Exact Mass.
     *
     */
    @Test 
    public void testGetTotalExactMass_IMolecularFormula() throws Exception{

		IMolecularFormula formula = new MolecularFormula();
		IIsotope carb = builder.newIsotope("C");
    	carb.setExactMass(12.00);
    	IIsotope cl = builder.newIsotope("Cl");
        cl.setExactMass(34.96885268);
        
        formula.addIsotope(carb);
        formula.addIsotope(cl);
    	
    	double totalExactMass = MolecularFormulaManipulator.getTotalExactMass(formula);

        Assert.assertEquals(46.96885268,totalExactMass,0.000001);
    }
	/**
     * A unit test suite for JUnit. Test total Exact Mass.
     *
     * @throws IOException
     * @throws ClassNotFoundException
     * @throws CDKException
     */
    @Test 
    public void testGetTotalExactMassWithCharge_IMolecularFormula() throws Exception{

		IMolecularFormula formula = MolecularFormulaManipulator.getMajorIsotopeMolecularFormula("CH5O", builder);
        
    	double totalExactMass = MolecularFormulaManipulator.getTotalExactMass(formula);
        Assert.assertEquals(33.034040,totalExactMass,0.0001);
        
        formula.setCharge(1);
    	double totalExactMass2 = MolecularFormulaManipulator.getTotalExactMass(formula);
        Assert.assertEquals(33.03349,totalExactMass2,0.0001);
    }
	/**
     * A unit test suite for JUnit. Test total Exact Mass.
     *
     * @throws IOException
     * @throws ClassNotFoundException
     * @throws CDKException
     */
    @Test 
    public void testGetTotalExactMassWithChargeNeg_IMolecularFormula() throws IOException, ClassNotFoundException {

		IMolecularFormula formula = MolecularFormulaManipulator.getMajorIsotopeMolecularFormula("H2PO4", builder);
        formula.setCharge(-1);
    	double totalExactMass2 = MolecularFormulaManipulator.getTotalExactMass(formula);
        Assert.assertEquals(96.96961875390926,totalExactMass2,0.0001);
    }
    @Test 
    public void testGetNaturalExactMass_IMolecularFormula() throws Exception {
		IMolecularFormula formula = new MolecularFormula();
        formula.addIsotope(builder.newIsotope("C"));
        formula.addIsotope(builder.newIsotope("Cl"));
    	
        double expectedMass = 0.0;
        expectedMass += IsotopeFactory.getInstance(builder).getNaturalMass(builder.newElement("C"));
        expectedMass += IsotopeFactory.getInstance(builder).getNaturalMass(builder.newElement("Cl"));
        
    	double totalExactMass = MolecularFormulaManipulator.getNaturalExactMass(formula);
        Assert.assertEquals(expectedMass, totalExactMass, 0.000001);
    }
    @Test 
    public void testGetTotalMassNumber_IMolecularFormula() throws Exception {
		IMolecularFormula formula = new MolecularFormula();
        formula.addIsotope(builder.newIsotope("C"));
        formula.addIsotope(builder.newIsotope("O"));
    	
        double totalExactMass = MolecularFormulaManipulator.getTotalMassNumber(formula);
        Assert.assertEquals(28, totalExactMass, 0.000001);
    }
    @Test
    public void testGetMajorIsotopeMass_IMolecularFormula() throws Exception {
        IMolecularFormula formula = new MolecularFormula();
        formula.addIsotope(builder.newIsotope("C"));
        formula.addIsotope(builder.newIsotope("H"), 4);

        double expectedMass = 0.0;
        expectedMass += IsotopeFactory.getInstance(builder).getMajorIsotope("C").getExactMass();
        expectedMass += 4.0*IsotopeFactory.getInstance(builder).getMajorIsotope("H").getExactMass();

        double totalExactMass = MolecularFormulaManipulator.getMajorIsotopeMass(formula);
        Assert.assertEquals(expectedMass, totalExactMass, 0.000001);
    }

    /**
     * A unit test suite for JUnit. Test total Exact Mass. It is 
     * necessary to have added the corresponding isotope before to calculate
     * the exact mass.
     *
     */
    @Test 
    public void testBug_1944604() throws Exception{

		IMolecularFormula formula = new MolecularFormula();
		IIsotope carb = builder.newIsotope("C");
    	
        formula.addIsotope(carb);
    	
        Assert.assertEquals("C1",MolecularFormulaManipulator.getString(formula,true));
        
    	double totalExactMass = MolecularFormulaManipulator.getTotalExactMass(formula);
    	
        Assert.assertEquals(12.0,totalExactMass,0.000001);
    }
    /**
     * A unit test suite for JUnit. Test total natural abundance.
     *
     */
    @Test 
    public void testGetTotalNaturalAbundance_IMolecularFormula() throws Exception{

		IMolecularFormula formula = new MolecularFormula();
		IIsotope carb = builder.newIsotope("C");
    	carb.setNaturalAbundance(98.93);
    	IIsotope cl = builder.newIsotope("Cl");
        cl.setNaturalAbundance(75.78);
        formula.addIsotope(carb);
        formula.addIsotope(cl);
    	
        double totalAbudance = MolecularFormulaManipulator.getTotalNaturalAbundance(formula);

        Assert.assertEquals(0.74969154,totalAbudance,0.000001);
    }
    
    /**
     * A unit test suite for JUnit. Test total natural abundance.
     *
     */
    @Test 
    public void testGetTotalNaturalAbundance_IMolecularFormula2() throws Exception{

		IMolecularFormula formula1 = new MolecularFormula();
		IIsotope br1 = builder.newIsotope("Br");
    	br1.setNaturalAbundance(49.31);
    	IIsotope br2 = builder.newIsotope("Br");
        br2.setNaturalAbundance(50.69);
        formula1.addIsotope(br1);
        formula1.addIsotope(br2);
    	
        Assert.assertEquals(2,formula1.getIsotopeCount(),0.000001);
        double totalAbudance = MolecularFormulaManipulator.getTotalNaturalAbundance(formula1);
        Assert.assertEquals(0.24995235,totalAbudance,0.000001);
    }
        
    @Test 
    public void testGetTotalNaturalAbundance_IMolecularFormula3() throws Exception{
        IMolecularFormula formula2 = new MolecularFormula();
		IIsotope br1 = builder.newIsotope("Br");
    	br1.setNaturalAbundance(50.69);
    	IIsotope br2 = builder.newIsotope("Br");
        br2.setNaturalAbundance(50.69);
        formula2.addIsotope(br1);
        formula2.addIsotope(br2);

        Assert.assertEquals(1,formula2.getIsotopeCount(),0.000001);
        double totalAbudance = MolecularFormulaManipulator.getTotalNaturalAbundance(formula2);

        Assert.assertEquals(0.25694761,totalAbudance,0.000001);
    }
    @Test 
    public void testGetTotalNaturalAbundance_IMolecularFormula4() throws Exception{
        IMolecularFormula formula2 = new MolecularFormula();
		IIsotope br1 = builder.newIsotope("Br");
    	br1.setNaturalAbundance(50.69);
        formula2.addIsotope(br1);
        formula2.addIsotope(br1);

        Assert.assertEquals(1,formula2.getIsotopeCount());
        double totalAbudance = MolecularFormulaManipulator.getTotalNaturalAbundance(formula2);

        Assert.assertEquals(0.25694761,totalAbudance,0.000001);
    }
    /**
	 * A unit test suite for JUnit.
	 *
	 * @return    The test suite
	 */
    @Test 
    public void testGetDBE_IMolecularFormula()throws Exception{
    	IMolecularFormula formula = new MolecularFormula();
    	formula.addIsotope(builder.newIsotope("C"),10);
    	formula.addIsotope(builder.newIsotope("H"),22);
		
    	Assert.assertEquals(0.0,MolecularFormulaManipulator.getDBE(formula), 0.01);
    	
    	formula = new MolecularFormula();
    	formula.addIsotope(builder.newIsotope("C"),10);
    	formula.addIsotope(builder.newIsotope("H"),16);
		
    	Assert.assertEquals(3.0,MolecularFormulaManipulator.getDBE(formula), 0.01);
    	

    	formula = new MolecularFormula();
    	formula.addIsotope(builder.newIsotope("C"),10);
    	formula.addIsotope(builder.newIsotope("H"),16);
    	formula.addIsotope(builder.newIsotope("O"));
		
    	Assert.assertEquals(3.0,MolecularFormulaManipulator.getDBE(formula), 0.01);
    	

    	formula = new MolecularFormula();
    	formula.addIsotope(builder.newIsotope("C"),10);
    	formula.addIsotope(builder.newIsotope("H"),19);
    	formula.addIsotope(builder.newIsotope("N"));

    	Assert.assertEquals(2.0,MolecularFormulaManipulator.getDBE(formula), 0.01);

    }
    /**
	 * A unit test suite for JUnit.
	 *
	 * @return    The test suite
	 */
    @Test 
    public void testGetHTML_IMolecularFormula() {
    	MolecularFormula formula = new MolecularFormula();
    	formula.addIsotope(builder.newIsotope("C"),8);
    	formula.addIsotope(builder.newIsotope("H"),10);
    	formula.addIsotope(builder.newIsotope("Cl"),2);
    	formula.addIsotope(builder.newIsotope("O"),2);
    	
        Assert.assertEquals("C<sub>8</sub>H<sub>10</sub>O<sub>2</sub>Cl<sub>2</sub>", MolecularFormulaManipulator.getHTML(formula));
    }
    /**
	 * A unit test suite for JUnit.
	 *
	 * @return    The test suite
	 */
    @Test 
    public void testGetHTML_IMolecularFormula_boolean_boolean() {
    	MolecularFormula formula = new MolecularFormula();
    	formula.addIsotope(builder.newIsotope("C"),10);
    	
        Assert.assertEquals("C<sub>10</sub>", MolecularFormulaManipulator.getHTML(formula,true,false));
        formula.setCharge(1);
        Assert.assertEquals("C<sub>10</sub><sup>1+</sup>", MolecularFormulaManipulator.getHTML(formula,true,false));
        formula.setCharge(formula.getCharge() - 2);
        Assert.assertEquals("C<sub>10</sub><sup>1-</sup>", MolecularFormulaManipulator.getHTML(formula,true,false));
    }
    /**
	 * A unit test suite for JUnit.
	 *
	 * @return    The test suite
	 */
    @Test 
    public void testGetHTML_IMolecularFormulaWhitIsotope() {
    	MolecularFormula formula = new MolecularFormula();
    	formula.addIsotope(ifac.getMajorIsotope("C"),2);
    	formula.addIsotope(ifac.getMajorIsotope("H"),6);
        Assert.assertEquals("<sup>12</sup>C<sub>2</sub><sup>1</sup>H<sub>6</sub>", MolecularFormulaManipulator.getHTML(formula,false,true));
    }

    /**
	 * A unit test suite for JUnit.
	 *
	 * @return    The test suite
	 */
    @Test 
    public void testGetHTML_IMolecularFormulaWhitIsotopeAndCharge() {
    	MolecularFormula formula = new MolecularFormula();
    	formula.addIsotope(ifac.getMajorIsotope("C"),2);
    	formula.addIsotope(ifac.getMajorIsotope("H"),6);
    	formula.setCharge(1);
        Assert.assertEquals("<sup>12</sup>C<sub>2</sub><sup>1</sup>H<sub>6</sub><sup>1+</sup>", MolecularFormulaManipulator.getHTML(formula,true,true));
    }
    /**
	 * A unit test suite for JUnit.
	 *
	 * @return    The test suite
	 */
	@Test 
    public void testGetMolecularFormula_IAtomContainer(){
		IAtomContainer ac = builder.newAtomContainer();
		ac.addAtom(builder.newAtom("C"));
		ac.addAtom(builder.newAtom("C"));
		ac.addAtom(builder.newAtom("H"));
		ac.addAtom(builder.newAtom("H"));
		ac.addAtom(builder.newAtom("H"));
		ac.addAtom(builder.newAtom("H"));

		IMolecularFormula mf1 = MolecularFormulaManipulator.getMolecularFormula(ac);		

		IMolecularFormula mf2 = new MolecularFormula();
		mf2.addIsotope(builder.newIsotope("C"),2);
		mf2.addIsotope(builder.newIsotope("H"),4);
		

		Assert.assertEquals(MolecularFormulaManipulator.getAtomCount(mf2), MolecularFormulaManipulator.getAtomCount(mf1));
		Assert.assertEquals(mf2.getIsotopeCount(), mf1.getIsotopeCount());
		IElement elemC = builder.newElement("C");
		IElement elemH = builder.newElement("H");
		Assert.assertEquals(mf2.getIsotopeCount(builder.newIsotope(elemC)), mf1.getIsotopeCount(builder.newIsotope(elemC)));
		Assert.assertEquals(mf2.getIsotopeCount(builder.newIsotope(elemH)), mf1.getIsotopeCount(builder.newIsotope(elemH)));
		Assert.assertEquals(MolecularFormulaManipulator.getElementCount(mf2, elemC), MolecularFormulaManipulator.getElementCount(mf1,elemC));
		Assert.assertEquals(MolecularFormulaManipulator.getElementCount(mf2, elemH), MolecularFormulaManipulator.getElementCount(mf1,elemH));
				
	}

    /**
	 * A unit test suite for JUnit.
	 *
	 * @return    The test suite
	 */
	@Test 
    public void testGetMolecularFormula_IAtomContainer_withCharge(){
		IAtomContainer ac = builder.newAtomContainer();
		ac.addAtom(builder.newAtom("C"));
		ac.getAtom(0).setFormalCharge(1);
		ac.addAtom(builder.newAtom("C"));
		ac.addAtom(builder.newAtom("H"));
		ac.addAtom(builder.newAtom("H"));
		ac.addAtom(builder.newAtom("H"));

		IMolecularFormula mf1 = MolecularFormulaManipulator.getMolecularFormula(ac);		

		Assert.assertEquals(1, mf1.getCharge(),0.000);
	}
	/**
	 * A unit test suite for JUnit.
	 *
	 * @return    The test suite
	 */
	@Test 
    public void testGetMolecularFormula_IAtomContainer_IMolecularFormula(){
		IAtomContainer ac = builder.newAtomContainer();
		ac.addAtom(builder.newAtom("C"));
		ac.addAtom(builder.newAtom("C"));
		ac.addAtom(builder.newAtom("H"));
		ac.addAtom(builder.newAtom("H"));
		ac.addAtom(builder.newAtom("H"));
		ac.addAtom(builder.newAtom("H"));

		IMolecularFormula mf1 = MolecularFormulaManipulator.getMolecularFormula(ac,new MolecularFormula());		


		IMolecularFormula mf2 = new MolecularFormula();
		mf2.addIsotope(builder.newIsotope("C"),2);
		mf2.addIsotope(builder.newIsotope("H"),4);

		Assert.assertEquals(MolecularFormulaManipulator.getAtomCount(mf2), MolecularFormulaManipulator.getAtomCount(mf1));
		Assert.assertEquals(mf2.getIsotopeCount(), mf1.getIsotopeCount());
		IElement elemC = builder.newElement("C");
		IElement elemH = builder.newElement("H");
		Assert.assertEquals(mf2.getIsotopeCount(builder.newIsotope(elemC)), mf1.getIsotopeCount(builder.newIsotope(elemC)));
		Assert.assertEquals(mf2.getIsotopeCount(builder.newIsotope(elemH)), mf1.getIsotopeCount(builder.newIsotope(elemH)));
		Assert.assertEquals(MolecularFormulaManipulator.getElementCount(mf2, elemC), MolecularFormulaManipulator.getElementCount(mf1,elemC));
		Assert.assertEquals(MolecularFormulaManipulator.getElementCount(mf2, elemH), MolecularFormulaManipulator.getElementCount(mf1,elemH));
						
	}
	
	/**
	 * A unit test suite for JUnit.
	 *
	 * @return    The test suite
	 */
	@Test 
    public void testGetMolecularFormula_IAtomContainerIMolecularFormula_2(){
		IAtomContainer ac = builder.newAtomContainer();
		ac.addAtom(builder.newAtom("C"));
		ac.addAtom(builder.newAtom("C"));
		ac.addAtom(builder.newAtom("H"));
		ac.addAtom(builder.newAtom("H"));
		ac.addAtom(builder.newAtom("H"));
		ac.addAtom(builder.newAtom("H"));

		IMolecularFormula mf0 = new MolecularFormula();
		mf0.addIsotope(builder.newIsotope("C"),2);
		mf0.addIsotope(builder.newIsotope("H"),5);
		
		IMolecularFormula mf1 = MolecularFormulaManipulator.getMolecularFormula(ac,mf0);		

		IMolecularFormula mf2 = new MolecularFormula();
		mf2.addIsotope(builder.newIsotope("C"),4);
		mf2.addIsotope(builder.newIsotope("H"),9);
		
		Assert.assertEquals(MolecularFormulaManipulator.getAtomCount(mf2), MolecularFormulaManipulator.getAtomCount(mf1));
		Assert.assertEquals(mf2.getIsotopeCount(), mf1.getIsotopeCount());
		IElement elemC = builder.newElement("C");
		IElement elemH = builder.newElement("H");
		Assert.assertEquals(mf2.getIsotopeCount(builder.newIsotope(elemC)), mf1.getIsotopeCount(builder.newIsotope(elemC)));
		Assert.assertEquals(mf2.getIsotopeCount(builder.newIsotope(elemH)), mf1.getIsotopeCount(builder.newIsotope(elemH)));
		Assert.assertEquals(MolecularFormulaManipulator.getElementCount(mf2, elemC), MolecularFormulaManipulator.getElementCount(mf1,elemC));
		Assert.assertEquals(MolecularFormulaManipulator.getElementCount(mf2, elemH), MolecularFormulaManipulator.getElementCount(mf1,elemH));
				
	}
	

	/**
	 * A unit test suite for JUnit.
	 *
	 * @return    The test suite
	 */
	@Test 
    public void testGetAtomContainer_IMolecularFormula(){
		

		IMolecularFormula mf2 = new MolecularFormula();
		mf2.addIsotope(builder.newIsotope("C"),2);
		mf2.addIsotope(builder.newIsotope("H"),4);
		
		IAtomContainer ac = MolecularFormulaManipulator.getAtomContainer(mf2);		


		Assert.assertEquals(6, ac.getAtomCount());
				
	}
	
	/**
	 * A unit test suite for JUnit.
	 *
	 * @return    The test suite
	 */
	@Test 
    public void testGetAtomContainer_IMolecularFormula_IAtomContainer(){

		IMolecularFormula mf2 = new MolecularFormula();
		mf2.addIsotope(builder.newIsotope("C"),2);
		mf2.addIsotope(builder.newIsotope("H"),4);
		
		IAtomContainer ac = MolecularFormulaManipulator.getAtomContainer(mf2,builder.newAtomContainer());		


		Assert.assertEquals(6, ac.getAtomCount());
				
	}
	
	
	/**
	 * A unit test suite for JUnit.
	 *
	 * @return    The test suite
	 */
	@Test 
    public void testMolecularFormulaIAtomContainer_to_IAtomContainer2(){
		IAtomContainer ac = builder.newAtomContainer();
		ac.addAtom(builder.newAtom("C"));
		ac.addAtom(builder.newAtom("C"));
		ac.addAtom(builder.newAtom("H"));
		ac.addAtom(builder.newAtom("H"));
		ac.addAtom(builder.newAtom("H"));
		ac.addAtom(builder.newAtom("H"));

		IMolecularFormula mf2 = new MolecularFormula();
		mf2.addIsotope(builder.newIsotope("C"),2);
		mf2.addIsotope(builder.newIsotope("H"),4);
		
		IAtomContainer ac2 = MolecularFormulaManipulator.getAtomContainer(mf2,builder.newAtomContainer());		


		Assert.assertEquals(ac2.getAtomCount(), ac2.getAtomCount());
		Assert.assertEquals(ac2.getAtom(0).getSymbol(),ac2.getAtom(0).getSymbol());
		Assert.assertEquals(ac2.getAtom(5).getSymbol(),ac2.getAtom(5).getSymbol());
				
	}
	
	/**
	 * A unit test suite for JUnit.
	 *
	 * @return    The test suite
	 */
	@Test 
    public void testElements_IMolecularFormula(){
		
		IMolecularFormula formula = new MolecularFormula();
    	formula.addIsotope(builder.newIsotope("C"),1);
    	formula.addIsotope(builder.newIsotope("H"),2);
    	
    	IIsotope br1 = builder.newIsotope("Br");
    	br1.setNaturalAbundance(50.69);
        formula.addIsotope(br1);
    	IIsotope br2 = builder.newIsotope("Br");
        br2.setNaturalAbundance(50.69);
        formula.addIsotope(br2);
        
        List<IElement> elements = MolecularFormulaManipulator.elements(formula);
        
        Assert.assertEquals(5, MolecularFormulaManipulator.getAtomCount(formula));
        Assert.assertEquals(3, elements.size());
	}
	
	/**
	 * A unit test suite for JUnit.
	 *
	 * @return    The test suite
	 */
	@Test 
    public void testCompare_Charge(){
		
		IMolecularFormula formula1 = new MolecularFormula();
    	formula1.addIsotope(builder.newIsotope("C"),1);
    	formula1.addIsotope(builder.newIsotope("H"),2);
    	
    	IMolecularFormula formula2 = new MolecularFormula();
    	formula2.addIsotope(builder.newIsotope("C"),1);
    	formula2.addIsotope(builder.newIsotope("H"),2);
    	
    	IMolecularFormula formula3 = new MolecularFormula();
    	formula3.addIsotope(builder.newIsotope("C"),1);
    	formula3.addIsotope(builder.newIsotope("H"),2);
    	formula3.setCharge(0);
    	
    	Assert.assertTrue(MolecularFormulaManipulator.compare(formula1, formula2));
    	Assert.assertFalse(MolecularFormulaManipulator.compare(formula1, formula3));
        
	}
	/**
	 * A unit test suite for JUnit.
	 *
	 * @return    The test suite
	 */
	@Test 
    public void testCompare_NumberIsotope(){
		
		IMolecularFormula formula1 = new MolecularFormula();
    	formula1.addIsotope(builder.newIsotope("C"),1);
    	formula1.addIsotope(builder.newIsotope("H"),2);
    	
    	IMolecularFormula formula2 = new MolecularFormula();
    	formula2.addIsotope(builder.newIsotope("C"),1);
    	formula2.addIsotope(builder.newIsotope("H"),2);
    	
    	IMolecularFormula formula3 = new MolecularFormula();
    	formula3.addIsotope(builder.newIsotope("C"),1);
    	formula3.addIsotope(builder.newIsotope("H"),3);
    	
    	Assert.assertTrue(MolecularFormulaManipulator.compare(formula1, formula2));
    	Assert.assertFalse(MolecularFormulaManipulator.compare(formula1, formula3));
        
	}
	/**
	 * A unit test suite for JUnit.
	 *
	 * @return    The test suite
	 */
	@Test 
    public void testCompare_IMolecularFormula_IMolecularFormula(){
		
		IMolecularFormula formula1 = new MolecularFormula();
    	formula1.addIsotope(builder.newIsotope("C"),1);
    	formula1.addIsotope(builder.newIsotope("H"),2);
    	
    	IMolecularFormula formula2 = new MolecularFormula();
    	formula2.addIsotope(builder.newIsotope("C"),1);
    	formula2.addIsotope(builder.newIsotope("H"),2);
    	
    	IMolecularFormula formula3 = new MolecularFormula();
    	formula3.addIsotope(builder.newIsotope("C"),1);
    	IIsotope hyd = builder.newIsotope("H");
    	hyd.setExactMass(2.002334234);
    	formula3.addIsotope(hyd,2);
    	
    	Assert.assertTrue(MolecularFormulaManipulator.compare(formula1, formula2));
    	Assert.assertFalse(MolecularFormulaManipulator.compare(formula1, formula3));
        
	}
	/**
	 * A unit test suite for JUnit.
	 *
	 * @return    The test suite
	 */
	@Test public void testGetHeavyElements_IMolecularFormula() {
		IMolecularFormula formula = new MolecularFormula();
		formula.addIsotope(builder.newIsotope("C"),10);
		formula.addIsotope(builder.newIsotope("H"),16);
        Assert.assertEquals(1, MolecularFormulaManipulator.getHeavyElements(formula).size());
    }

	/**
	 * A unit test suite for JUnit.
	 *
	 * @return    The test suite
	 */
    @Test public void testGetHeavyElements_IMolecularFormula_2() {
    	IMolecularFormula formula = MolecularFormulaManipulator.getMolecularFormula("CH3OH", builder);
        Assert.assertEquals(2, MolecularFormulaManipulator.getHeavyElements(formula).size());
    }

	/**
	 * A unit test suite for JUnit. Not null.
	 *
	 * @return    The test suite
	 */
	@Test 
    public void testGenerateOrderEle() {
		String[] listElements = new String[]{
			    "C", "H", "O", "N", "Si", "P", "S", "F", "Cl",
			    "Br", "I", "Sn", "B", "Pb", "Tl", "Ba", "In", "Pd",
			    "Pt", "Os", "Ag", "Zr", "Se", "Zn", "Cu", "Ni", "Co", 
			    "Fe", "Cr", "Ti", "Ca", "K", "Al", "Mg", "Na", "Ce",
			    "Hg", "Au", "Ir", "Re", "W", "Ta", "Hf", "Lu", "Yb", 
			    "Tm", "Er", "Ho", "Dy", "Tb", "Gd", "Eu", "Sm", "Pm",
			    "Nd", "Pr", "La", "Cs", "Xe", "Te", "Sb", "Cd", "Rh", 
			    "Ru", "Tc", "Mo", "Nb", "Y", "Sr", "Rb", "Kr", "As", 
			    "Ge", "Ga", "Mn", "V", "Sc", "Ar", "Ne", "Be", "Li", 
			    "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", 
			    "Th", "Pa", "U", "Np", "Pu","R"};
		
		String[] listGenerated = MolecularFormulaManipulator.generateOrderEle();
		Assert.assertEquals(listElements.length,listGenerated.length);
		
		for(int i = 0 ; i < listElements.length; i++)
			Assert.assertEquals(listElements[i],listGenerated[i]);
		
	}
	/**
	 * A unit test suite for JUnit. Not null.
	 * TODO: REACT: Introduce method
	 * 
	 * @return    The test suite
     * @cdk.bug 2672696
	 */
    @Test
    public void testGetHillString_IMolecularFormula() {
        IMolecularFormula formula = MolecularFormulaManipulator.getMolecularFormula("CH3OH", builder);
        String listGenerated = MolecularFormulaManipulator.getHillString(formula);
        Assert.assertEquals("CH4O", listGenerated);

        formula = MolecularFormulaManipulator.getMolecularFormula("CH3CH2Br", builder);
        listGenerated = MolecularFormulaManipulator.getHillString(formula);
        Assert.assertEquals("C2H5Br", listGenerated);

        formula = MolecularFormulaManipulator.getMolecularFormula("HCl", builder);
        listGenerated = MolecularFormulaManipulator.getHillString(formula);
        Assert.assertEquals("ClH", listGenerated);

        formula = MolecularFormulaManipulator.getMolecularFormula("HBr", builder);
        listGenerated = MolecularFormulaManipulator.getHillString(formula);
        Assert.assertEquals("BrH", listGenerated);
    }

    /**
     * @cdk.bug 1944604
     * @throws InvalidSmilesException
     */
    @Test
    public void testSingleAtomFromSmiles() throws InvalidSmilesException {
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer mol = sp.parseSmiles("C");
        IMolecularFormula mf = MolecularFormulaManipulator.getMolecularFormula(mol);
        double exactMass = MolecularFormulaManipulator.getTotalExactMass(mf);
        Assert.assertEquals(12.0000, exactMass, 0.0001);
    }

    /**
     */
    @Test
    public void testSingleAtom() {
        String formula = "CH4";
        IMolecularFormula mf = MolecularFormulaManipulator.getMolecularFormula(formula, builder);
    	Assert.assertEquals(1, MolecularFormulaManipulator.getIsotopes(mf, mf.getBuilder().newElement("C")).size());
    }

    /**
     */
    @Test
    public void testSimplifyMolecularFormula_String() {
        String formula = "C1H41.H2O";
    	String simplifyMF = MolecularFormulaManipulator.simplifyMolecularFormula(formula);
        Assert.assertEquals("C1H43O", simplifyMF);
    }
    /**
     */
    @Test
    public void testSimplifyMolecularFormula_String2() {
        String formula = "CH41.H2O";
    	String simplifyMF = MolecularFormulaManipulator.simplifyMolecularFormula(formula);
        Assert.assertEquals("CH43O", simplifyMF);
    }
    /**
     */
    @Test
    public void testSimplifygetMF() {
        String formula = "CH4.H2O";
    	IMolecularFormula formula1 = new MolecularFormula();
    	formula1.addIsotope(builder.newIsotope("C"),1);
    	formula1.addIsotope(builder.newIsotope("H"),6);
    	formula1.addIsotope(builder.newIsotope("O"),1);
    	IMolecularFormula ff = MolecularFormulaManipulator.getMolecularFormula(formula, builder);
    	Assert.assertTrue(MolecularFormulaManipulator.compare(formula1, MolecularFormulaManipulator.getMolecularFormula(formula, builder)));
    	Assert.assertEquals("CH6O",MolecularFormulaManipulator.getString(ff));
    }

      /**
       */
      @Test
      public void testSpace() {
          String formula = "C17H21NO. C7H6O3";
          String simplifyMF = MolecularFormulaManipulator.simplifyMolecularFormula(formula);
          Assert.assertEquals("C24H27NO4", simplifyMF);
      }
    /**
     */
    @Test
    public void test0() {
        String formula = "Fe.(C6H11O7)3";
        String simplifyMF = MolecularFormulaManipulator.simplifyMolecularFormula(formula);
        Assert.assertEquals("FeC18H33O21", simplifyMF);
    }
      /**
       */
      @Test
      public void test1() {
          String formula = "(C6H11O7)3.Fe";
          String simplifyMF = MolecularFormulaManipulator.simplifyMolecularFormula(formula);
          Assert.assertEquals("C18H33O21Fe", simplifyMF);
      }

	  /**
	   */
	  @Test
	  public void test2() {
	      String formula = "C14H14N2.2HCl";
	      String simplifyMF = MolecularFormulaManipulator.simplifyMolecularFormula(formula);
	      Assert.assertEquals("C14H16N2Cl2", simplifyMF);
	    }
	  /**
	   */
	  @Test
      public void test3() {
          String formula = "(C27H33N3O8)2.2HNO3.3H2O";
          String simplifyMF = MolecularFormulaManipulator.simplifyMolecularFormula(formula);
          Assert.assertEquals("C54H74N8O25", simplifyMF);
      }

	  /**
	   */
	  @Test
      public void test4() {
          String formula = "(C27H33N3O8)2.2HNO3.3H2O";
          IMolecularFormula formula1 = new MolecularFormula();
          formula1.addIsotope(builder.newIsotope("C"),54);
      	  formula1.addIsotope(builder.newIsotope("H"),74);
      	  formula1.addIsotope(builder.newIsotope("O"),25);
      	  formula1.addIsotope(builder.newIsotope("N"),8);
      	  IMolecularFormula ff = MolecularFormulaManipulator.getMolecularFormula(formula, builder);
      	  Assert.assertTrue(MolecularFormulaManipulator.compare(formula1, MolecularFormulaManipulator.getMolecularFormula(formula, builder)));
      	  Assert.assertEquals("C54H74N8O25", MolecularFormulaManipulator.getString(ff));
      }
	  /**
	   */
	  @Test
      public void test5() {
          String formula = "[SO3]2-";
          IMolecularFormula formula1 = new MolecularFormula();
          formula1.addIsotope(builder.newIsotope("S"),1);
      	  formula1.addIsotope(builder.newIsotope("O"),3);
      	  formula1.setCharge(-2);
      	  IMolecularFormula ff = MolecularFormulaManipulator.getMolecularFormula(formula, builder);
      	  Assert.assertTrue(MolecularFormulaManipulator.compare(formula1, MolecularFormulaManipulator.getMolecularFormula(formula, builder)));
      	  Assert.assertEquals("O3S", MolecularFormulaManipulator.getString(ff));
      	  Assert.assertEquals(-2, ff.getCharge(),0.00001);
      }
	  /**
	   */
	  @Test
      public void test6() {
          String formula = "(CH3)2";
          IMolecularFormula formula1 = new MolecularFormula();
          formula1.addIsotope(builder.newIsotope("C"),2);
      	  formula1.addIsotope(builder.newIsotope("H"),6);
      	  IMolecularFormula ff = MolecularFormulaManipulator.getMolecularFormula(formula, builder);
      	  Assert.assertTrue(MolecularFormulaManipulator.compare(formula1, MolecularFormulaManipulator.getMolecularFormula(formula, builder)));
      	  Assert.assertEquals("C2H6", MolecularFormulaManipulator.getString(ff));
      }

	  /**
	   */
	  @Test
      public void testWithH_Initial() {
          String formula = "HC5H11NO2H";
          IMolecularFormula formula1 = new MolecularFormula();
          formula1.addIsotope(builder.newIsotope("C"),5);
      	  formula1.addIsotope(builder.newIsotope("H"),13);
      	  formula1.addIsotope(builder.newIsotope("N"),1);
      	  formula1.addIsotope(builder.newIsotope("O"),2);
      	  IMolecularFormula ff = MolecularFormulaManipulator.getMolecularFormula(formula, builder);
      	  Assert.assertTrue(MolecularFormulaManipulator.compare(formula1, MolecularFormulaManipulator.getMolecularFormula(formula, builder)));
      	  Assert.assertEquals("C5H13NO2", MolecularFormulaManipulator.getString(ff));
      }
}

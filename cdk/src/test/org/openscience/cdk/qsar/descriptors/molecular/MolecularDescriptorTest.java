/* $Revision$ $Author$ $Date$
 * 
 * Copyright (C) 2007  Egon Willighagen <egonw@users.sf.net>
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
package org.openscience.cdk.qsar.descriptors.molecular;

import org.junit.Assert;
import org.junit.Ignore;
import org.junit.Test;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.*;
import org.openscience.cdk.interfaces.IBond.Order;
import org.openscience.cdk.qsar.DescriptorValue;
import org.openscience.cdk.qsar.IMolecularDescriptor;
import org.openscience.cdk.qsar.descriptors.DescriptorTest;
import org.openscience.cdk.qsar.result.*;
import org.openscience.cdk.tools.diff.AtomContainerDiff;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import javax.vecmath.Point3d;

/**
 * Tests for molecular descriptors.
 *
 * @cdk.module test-qsarmolecular
 */
public abstract class MolecularDescriptorTest extends DescriptorTest {
	
	protected IMolecularDescriptor descriptor;

	public MolecularDescriptorTest() {}

	public void setDescriptor(Class descriptorClass) throws Exception {
		if (descriptor == null) {
			Object descriptor = descriptorClass.newInstance();
			if (!(descriptor instanceof IMolecularDescriptor)) {
				throw new CDKException("The passed descriptor class must be a IMolecularDescriptor");
			}
			this.descriptor = (IMolecularDescriptor)descriptor;
		}
		super.setDescriptor(descriptorClass);
	}

    @Test
    public void testCalculate_IAtomContainer() {
        IAtomContainer mol = null;
        try {
            mol = someoneBringMeSomeWater();
        } catch (Exception e) {
            e.printStackTrace();
            Assert.fail("Error in generating the test molecule");
        }

        DescriptorValue v = null;
        try {
            v = descriptor.calculate(mol);
        } catch (Exception e) {
            Assert.fail("A descriptor must not throw an exception. Exception was:\n"+e.getMessage());
        }
        Assert.assertNotNull(v);
        Assert.assertTrue(
        	"The descriptor did not calculate any value.",
        	0 != v.getValue().length()
        );
    }

    @Test
    public void testCalculate_NoModifications() throws Exception {
        IAtomContainer mol = someoneBringMeSomeWater();
        IAtomContainer clone = (IAtomContainer)mol.clone();
        descriptor.calculate(mol);
        String diff = AtomContainerDiff.diff(clone, mol); 
        Assert.assertEquals(
          "The descriptor must not change the passed molecule in any respect, but found this diff: " + diff,
          0, diff.length()
        );
    }

    /**
	 * Checks if the given labels are consistent.
	 * 
	 * @throws Exception Passed on from calculate.
	 */
    @Test
    public void testLabels() throws Exception {
        IAtomContainer mol = someoneBringMeSomeWater();
        
        DescriptorValue v = descriptor.calculate(mol);
        Assert.assertNotNull(v);
        String[] names = v.getNames();
        Assert.assertNotNull(
        	"The descriptor must return labels using the getNames() method.",
        	names
        );
        Assert.assertNotSame(
        	"At least one label must be given.",
        	0, names.length
        );
        for (int i=0; i<names.length; i++) {
        	Assert.assertNotNull(
        		"A descriptor label may not be null.",
        		names[i]
        	);
        	Assert.assertNotSame(
        		"The label string must not be empty.",
        		0, names[i].length()
        	);
//        	System.out.println("Label: " + names[i]);
        }
        Assert.assertNotNull(v.getValue());
        int valueCount = v.getValue().length();
        Assert.assertEquals(
        	"The number of labels must equals the number of values.",
        	names.length, valueCount
        );
    }

     /**
     * Check if the names obtained directly from the decsriptor without
     * calculation match those obtained from the descriptor value object.
     * Also ensure that the number of actual values matches the length
     * of the names
     */
     @Test
    public void testNamesConsistency() throws Exception {
        IAtomContainer mol = someoneBringMeSomeWater();

        String[] names1 = descriptor.getDescriptorNames();
        DescriptorValue v = descriptor.calculate( mol);
        String[] names2 = v.getNames();

        Assert.assertEquals(names1.length, names2.length);
        Assert.assertArrayEquals(names1, names2);

        int valueCount = v.getValue().length();
        Assert.assertEquals(valueCount, names1.length);
    }

    @Test
    public void testGetDescriptorResultType() throws Exception {
    	IDescriptorResult result = descriptor.getDescriptorResultType();
    	Assert.assertNotNull(
    		"The getDescriptorResultType() must not be null.",
    		result
    	);
    	
    	IAtomContainer mol = someoneBringMeSomeWater();
        DescriptorValue v = descriptor.calculate(mol);
        
    	Assert.assertTrue(
    		"The getDescriptorResultType() is inconsistent with the calculated descriptor results",
    		result.getClass().getName().contains(v.getValue().getClass().getName()) 
    	);
    	Assert.assertEquals(
    		"The specified getDescriptorResultType() length does not match the actually calculated result vector length",
    		v.getValue().length(), result.length()
    	);
    }

    @Test
    public void testTakeIntoAccountImplicitHydrogens() {
        IChemObjectBuilder builder = DefaultChemObjectBuilder.getInstance();
        IMolecule methane1 = builder.newMolecule();
        IAtom c1 = builder.newAtom("C");
        c1.setHydrogenCount(4);
        methane1.addAtom(c1);

        IMolecule methane2 = builder.newMolecule();
        IAtom c2 = builder.newAtom("C");
        methane2.addAtom(c2);
        IAtom h1 = builder.newAtom("H"); methane2.addAtom(h1);
        IAtom h2 = builder.newAtom("H"); methane2.addAtom(h2);
        IAtom h3 = builder.newAtom("H"); methane2.addAtom(h3);
        IAtom h4 = builder.newAtom("H"); methane2.addAtom(h4);
        methane2.addBond(0, 1, Order.SINGLE);
        methane2.addBond(0, 2, Order.SINGLE);
        methane2.addBond(0, 3, Order.SINGLE);
        methane2.addBond(0, 4, Order.SINGLE);

        IDescriptorResult v1 = descriptor.calculate(methane1).getValue();
        IDescriptorResult v2 = descriptor.calculate(methane2).getValue();

        String errorMessage = "("+descriptor.getClass().toString()+") The descriptor does not give the same results depending on whether hydrogens are implicit or explicit.";
        if (v1 instanceof IntegerResult) {
            Assert.assertEquals(errorMessage, ((IntegerResult)v1).intValue(), ((IntegerResult)v2).intValue());
        } else if (v1 instanceof DoubleResult) {
            Assert.assertEquals(errorMessage, ((DoubleResult)v1).doubleValue(), ((DoubleResult)v2).doubleValue(), 0.00001);
        } else if (v1 instanceof BooleanResult) {
            Assert.assertEquals(errorMessage, ((BooleanResult)v1).booleanValue(), ((BooleanResult)v2).booleanValue());
        } else if (v1 instanceof DoubleArrayResult) {
            DoubleArrayResult da1 = (DoubleArrayResult)v1;
            DoubleArrayResult da2 = (DoubleArrayResult)v2;
            for (int i=0; i<da1.length(); i++) {
                Assert.assertEquals(errorMessage, da1.get(i), da2.get(i), 0.00001);
            }
        } else if (v1 instanceof IntegerArrayResult) {
            IntegerArrayResult da1 = (IntegerArrayResult)v1;
            IntegerArrayResult da2 = (IntegerArrayResult)v2;
            for (int i=0; i<da1.length(); i++) {
                Assert.assertEquals(errorMessage, da1.get(i), da2.get(i));
            }
        }
    }

    @Test
    public void testTakeIntoAccountImplicitHydrogensInEthane() {
        IChemObjectBuilder builder = DefaultChemObjectBuilder.getInstance();
        IMolecule ethane1 = builder.newMolecule();
        IAtom c1 = builder.newAtom("C");
        IAtom c2 = builder.newAtom("C");
        c1.setHydrogenCount(3);
        c2.setHydrogenCount(3);
        ethane1.addAtom(c1);
        ethane1.addAtom(c2);
        ethane1.addBond(0,1,IBond.Order.SINGLE);

        IMolecule ethane2 = builder.newMolecule();
        IAtom c3 = builder.newAtom("C");
        IAtom c4 = builder.newAtom("C");
        ethane2.addAtom(c3);
        ethane2.addAtom(c4);

        IAtom h1 = builder.newAtom("H");
        ethane2.addAtom(h1);
        IAtom h2 = builder.newAtom("H");
        ethane2.addAtom(h2);
        IAtom h3 = builder.newAtom("H");
        ethane2.addAtom(h3);

        IAtom h4 = builder.newAtom("H");
        IAtom h5 = builder.newAtom("H");
        IAtom h6 = builder.newAtom("H");
        ethane2.addAtom(h4);
        ethane2.addAtom(h5);
        ethane2.addAtom(h6);

        ethane2.addBond(0, 1, Order.SINGLE);
        ethane2.addBond(0, 2, Order.SINGLE);
        ethane2.addBond(0, 3, Order.SINGLE);
        ethane2.addBond(0, 4, Order.SINGLE);

        ethane2.addBond(1, 5, Order.SINGLE);
        ethane2.addBond(1, 6, Order.SINGLE);
        ethane2.addBond(1, 7, Order.SINGLE);

        IDescriptorResult v1 = descriptor.calculate(ethane1).getValue();
        IDescriptorResult v2 = descriptor.calculate(ethane2).getValue();

        String errorMessage = "(" + descriptor.getClass().toString() + ") The descriptor does not give the same results depending on whether hydrogens are implicit or explicit.";
        if (v1 instanceof IntegerResult) {
            Assert.assertEquals(errorMessage, ((IntegerResult) v1).intValue(), ((IntegerResult) v2).intValue());
        } else if (v1 instanceof DoubleResult) {
            Assert.assertEquals(errorMessage, ((DoubleResult) v1).doubleValue(), ((DoubleResult) v2).doubleValue(), 0.00001);
        } else if (v1 instanceof BooleanResult) {
            Assert.assertEquals(errorMessage, ((BooleanResult) v1).booleanValue(), ((BooleanResult) v2).booleanValue());
        } else if (v1 instanceof DoubleArrayResult) {
            DoubleArrayResult da1 = (DoubleArrayResult) v1;
            DoubleArrayResult da2 = (DoubleArrayResult) v2;
            for (int i = 0; i < da1.length(); i++) {
                Assert.assertEquals(errorMessage, da1.get(i), da2.get(i), 0.00001);
            }
        } else if (v1 instanceof IntegerArrayResult) {
            IntegerArrayResult da1 = (IntegerArrayResult) v1;
            IntegerArrayResult da2 = (IntegerArrayResult) v2;
            for (int i = 0; i < da1.length(); i++) {
                Assert.assertEquals(errorMessage, da1.get(i), da2.get(i));
            }
        }
    }

    @Ignore
    @Test public void testTakeIntoAccountBondHybridization() {
        IChemObjectBuilder builder = DefaultChemObjectBuilder.getInstance();
        IMolecule ethane1 = builder.newMolecule();
        IAtom c1 = builder.newAtom("C");
        IAtom c2 = builder.newAtom("C");
        ethane1.addAtom(c1);
        ethane1.addAtom(c2);
        ethane1.addBond(0, 1, IBond.Order.DOUBLE);

        IMolecule ethane2 = builder.newMolecule();
        IAtom c3 = builder.newAtom("C");
        c3.setHybridization(IAtomType.Hybridization.SP2);
        IAtom c4 = builder.newAtom("C");
        c4.setHybridization(IAtomType.Hybridization.SP2);
        ethane2.addAtom(c3);
        ethane2.addAtom(c4);
        ethane2.addBond(0, 1, IBond.Order.SINGLE);

        IDescriptorResult v1 = descriptor.calculate(ethane1).getValue();
        IDescriptorResult v2 = descriptor.calculate(ethane2).getValue();

        String errorMessage = "("+descriptor.getClass().toString()+") The descriptor does not give the same results depending on whether bond order or atom type are considered.";
        if (v1 instanceof IntegerResult) {
            Assert.assertEquals(errorMessage, ((IntegerResult)v1).intValue(), ((IntegerResult)v2).intValue());
        } else if (v1 instanceof DoubleResult) {
            Assert.assertEquals(errorMessage, ((DoubleResult)v1).doubleValue(), ((DoubleResult)v2).doubleValue(), 0.00001);
        } else if (v1 instanceof BooleanResult) {
            Assert.assertEquals(errorMessage, ((BooleanResult)v1).booleanValue(), ((BooleanResult)v2).booleanValue());
        } else if (v1 instanceof DoubleArrayResult) {
            DoubleArrayResult da1 = (DoubleArrayResult)v1;
            DoubleArrayResult da2 = (DoubleArrayResult)v2;
            for (int i=0; i<da1.length(); i++) {
                Assert.assertEquals(errorMessage, da1.get(i), da2.get(i), 0.00001);
            }
        } else if (v1 instanceof IntegerArrayResult) {
            IntegerArrayResult da1 = (IntegerArrayResult)v1;
            IntegerArrayResult da2 = (IntegerArrayResult)v2;
            for (int i=0; i<da1.length(); i++) {
                Assert.assertEquals(errorMessage, da1.get(i), da2.get(i));
            }
        }
    }

    private IMolecule someoneBringMeSomeWater() throws Exception {
        IMolecule mol = DefaultChemObjectBuilder.getInstance().newMolecule();
        IAtom c1 = DefaultChemObjectBuilder.getInstance().newAtom("O");
        c1.setPoint3d(new Point3d(0.0, 0.0, 0.0));
        IAtom h1 = DefaultChemObjectBuilder.getInstance().newAtom("H");
        h1.setPoint3d(new Point3d(1.0, 0.0, 0.0));
        IAtom h2 = DefaultChemObjectBuilder.getInstance().newAtom("H");
        h2.setPoint3d(new Point3d(-1.0, 0.0, 0.0));
        mol.addAtom(c1);
        mol.addAtom(h1);
        mol.addAtom(h2);
        mol.addBond(0,1,IBond.Order.SINGLE);
        mol.addBond(0,2,IBond.Order.SINGLE);
        AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(mol);
        return mol;
    }
    
}

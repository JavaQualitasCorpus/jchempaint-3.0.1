/* $RCSfile$
 * $Author$
 * $Date$
 * $Revision$
 * 
 * Copyright (C) 2002-2007  The Chemistry Development Kit (CDK) project
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
 * 
 */
package org.openscience.cdk.interfaces;

import java.util.Iterator;

import javax.vecmath.Point2d;
import javax.vecmath.Point3d;

import org.junit.Assert;
import org.junit.Test;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.interfaces.IBond.Order;

/**
 * Checks the functionality of {@link IBond} implementations.
 *
 * @cdk.module test-interfaces
 */
public abstract class AbstractBondTest extends AbstractElectronContainerTest {

    @Test
    public void testCompare_Object() {
    	IBond b = (IBond)newChemObject();
        IAtom c = b.getBuilder().newAtom("C");
        IAtom o = b.getBuilder().newAtom("O");
        b.setAtom(c, 0); b.setAtom(o, 1); b.setOrder(Order.SINGLE);
        
        IBond b2 = (IBond)newChemObject();
        b2.setAtom(c, 0); b2.setAtom(o, 1); b2.setOrder(Order.SINGLE);

        Assert.assertTrue(b.compare(b2));
    }

    @Test
    public void testContains_IAtom() {
    	IBond b = (IBond)newChemObject();
        IAtom c = b.getBuilder().newAtom("C");
        IAtom o = b.getBuilder().newAtom("O");
        b.setAtom(c, 0); b.setAtom(o, 1); b.setOrder(Order.SINGLE);

        Assert.assertTrue(b.contains(c));
        Assert.assertTrue(b.contains(o));
    }

    @Test
    public void testGetAtomCount() {
    	IBond b = (IBond)newChemObject();
        IAtom c = b.getBuilder().newAtom("C");
        IAtom o = b.getBuilder().newAtom("O");
        b.setAtom(c, 0); b.setAtom(o, 1); b.setOrder(Order.SINGLE);

        Assert.assertEquals(2.0, b.getAtomCount(), 0.001);
    }

    @Test
    public void testSetAtoms_arrayIAtom() {
    	IBond b = (IBond)newChemObject();
        IAtom[] atomsToAdd = new IAtom[2];
        atomsToAdd[0] = b.getBuilder().newAtom("C");
        atomsToAdd[1] = b.getBuilder().newAtom("O");

        b.setAtoms(atomsToAdd);

        Assert.assertEquals(2, b.getAtomCount());
        Assert.assertEquals(atomsToAdd[0], b.getAtom(0));
        Assert.assertEquals(atomsToAdd[1], b.getAtom(1));
    }

    @Test
    public void testAtoms() {
    	IBond b = (IBond)newChemObject();
        IAtom c = b.getBuilder().newAtom("C");
        IAtom o = b.getBuilder().newAtom("O");
        b.setAtom(c, 0); b.setAtom(o, 1); b.setOrder(Order.SINGLE);

        Iterator<IAtom> atoms = b.atoms().iterator();
        Assert.assertEquals(2, b.getAtomCount());
        Assert.assertTrue(atoms.hasNext());
        Assert.assertEquals(c, atoms.next());
        Assert.assertTrue(atoms.hasNext());
        Assert.assertEquals(o, atoms.next());
        Assert.assertFalse(atoms.hasNext());
    }

    @Test
    public void testGetAtom_int() {
    	IBond b = (IBond)newChemObject();
        IAtom c = b.getBuilder().newAtom("C");
        IAtom o = b.getBuilder().newAtom("O");
        b.setAtom(c, 0); b.setAtom(o, 1); b.setOrder(Order.SINGLE);

        Assert.assertEquals(c, b.getAtom(0));
        Assert.assertEquals(o, b.getAtom(1));
    }

    @Test
    public void testSetAtom_IAtom_int() {
        IBond b = (IBond)newChemObject();
        IAtom c = b.getBuilder().newAtom("C");
        IAtom o = b.getBuilder().newAtom("O");

        b.setAtom(c, 0);
        b.setAtom(o, 1);

        Assert.assertEquals(c, b.getAtom(0));
        Assert.assertEquals(o, b.getAtom(1));
    }

    @Test
    public void testGetConnectedAtom_IAtom() {
    	IBond b = (IBond)newChemObject();
        IAtom c = b.getBuilder().newAtom("C");
        IAtom o = b.getBuilder().newAtom("O");
        b.setAtom(c, 0); b.setAtom(o, 1); b.setOrder(Order.SINGLE);

        Assert.assertEquals(c, b.getConnectedAtom(o));
        Assert.assertEquals(o, b.getConnectedAtom(c));

        // test default return value
        Assert.assertNull(b.getConnectedAtom(b.getBuilder().newAtom()));
    }

    @Test
    public void testGetConnectedAtoms_IAtom() {
    	IBond b = (IBond)newChemObject();
    	IAtom[] atoms = new IAtom[3];
        atoms[0] = b.getBuilder().newAtom("B");
        atoms[1] = b.getBuilder().newAtom("H");
        atoms[2] = b.getBuilder().newAtom("B");

        b.setAtoms(atoms);
        b.setOrder(IBond.Order.SINGLE); // C=O bond

        IAtom[] connectedAtoms = b.getConnectedAtoms(atoms[1]);
        Assert.assertNotNull(connectedAtoms);
        Assert.assertEquals(2, connectedAtoms.length);
        Assert.assertNotNull(connectedAtoms[0]);
        Assert.assertNotNull(connectedAtoms[1]);

        // test default return value
        connectedAtoms = b.getConnectedAtoms(b.getBuilder().newAtom());
        Assert.assertNull(connectedAtoms);
    }

    @Test
    public void testIsConnectedTo_IBond() {
    	IBond b = (IBond)newChemObject();
        IAtom c1 = b.getBuilder().newAtom("C");
        IAtom o = b.getBuilder().newAtom("O");
        IAtom c2 = b.getBuilder().newAtom("C");
        IAtom c3 = b.getBuilder().newAtom("C");

        IBond b1 = b.getBuilder().newBond(c1, o);
        IBond b2 = b.getBuilder().newBond(o, c2);
        IBond b3 = b.getBuilder().newBond(c2, c3);

        Assert.assertTrue(b1.isConnectedTo(b2));
        Assert.assertTrue(b2.isConnectedTo(b1));
        Assert.assertTrue(b2.isConnectedTo(b3));
        Assert.assertTrue(b3.isConnectedTo(b2));
        Assert.assertFalse(b1.isConnectedTo(b3));
        Assert.assertFalse(b3.isConnectedTo(b1));
    }

    @Test
    public void testGetOrder() {
    	IBond b = (IBond)newChemObject();
        IAtom c = b.getBuilder().newAtom("C");
        IAtom o = b.getBuilder().newAtom("O");
        b.setAtom(c, 0); b.setAtom(o, 1); b.setOrder(Order.DOUBLE);

        Assert.assertEquals(IBond.Order.DOUBLE, b.getOrder());
    }

    @Test
    public void testSetOrder_IBond_Order() {
    	IBond b = (IBond)newChemObject();
        IAtom c = b.getBuilder().newAtom("C");
        IAtom o = b.getBuilder().newAtom("O");
        b.setAtom(c, 0); b.setAtom(o, 1); b.setOrder(Order.DOUBLE);

        Assert.assertEquals(IBond.Order.DOUBLE, b.getOrder());

        b.setOrder(IBond.Order.SINGLE);
        Assert.assertEquals(IBond.Order.SINGLE, b.getOrder());
    }

    @Test
    public void testSetStereo_IBond_Stereo() {
    	IBond b = (IBond)newChemObject();
        IAtom c = b.getBuilder().newAtom("C");
        IAtom o = b.getBuilder().newAtom("O");
        b.setAtom(c, 0); b.setAtom(o, 1); b.setOrder(Order.DOUBLE);
        b.setStereo(IBond.Stereo.DOWN);
        Assert.assertEquals(IBond.Stereo.DOWN, b.getStereo());
        b.setStereo(IBond.Stereo.UP);
        Assert.assertEquals(IBond.Stereo.UP, b.getStereo());
    }

    @Test
    public void testGetStereo() {
    	IChemObject object = newChemObject(); 
        IAtom c = object.getBuilder().newAtom("C");
        IAtom o = object.getBuilder().newAtom("O");

        IBond b = object.getBuilder().newBond(c, o, IBond.Order.DOUBLE, IBond.Stereo.UP);
        Assert.assertEquals(IBond.Stereo.UP, b.getStereo());
    }

    @Test
    public void testGet2DCenter() {
    	IChemObject object = newChemObject(); 
        IAtom o = object.getBuilder().newAtom("O", new Point2d(0.0, 0.0));
        IAtom c = object.getBuilder().newAtom("C", new Point2d(1.0, 1.0));
        IBond b = object.getBuilder().newBond(c, o);

        Assert.assertEquals(0.5, b.get2DCenter().x, 0.001);
        Assert.assertEquals(0.5, b.get2DCenter().y, 0.001);
    }

    @Test
    public void testGet3DCenter() {
    	IChemObject object = newChemObject(); 
        IAtom o = object.getBuilder().newAtom("O", new Point3d(0.0, 0.0, 0.0));
        IAtom c = object.getBuilder().newAtom("C", new Point3d(1.0, 1.0, 1.0));
        IBond b = object.getBuilder().newBond(c, o);

        Assert.assertEquals(0.5, b.get3DCenter().x, 0.001);
        Assert.assertEquals(0.5, b.get3DCenter().y, 0.001);
        Assert.assertEquals(0.5, b.get3DCenter().z, 0.001);
    }

    @Test
    public void testClone() throws Exception {
        IBond bond = (IBond)newChemObject();
        Object clone = bond.clone();
        Assert.assertNotNull(clone);
        Assert.assertTrue(clone instanceof org.openscience.cdk.interfaces.IBond);
    }

    @Test
    public void testClone_IAtom() throws Exception {
    	IChemObject object = newChemObject(); 
        IAtom atom1 = object.getBuilder().newAtom("C");
        IAtom atom2 = object.getBuilder().newAtom("O");
        IBond bond = object.getBuilder().newBond(atom1, atom2);
        IBond clone = (IBond) bond.clone();

        // test cloning of atoms
        Assert.assertNotSame(atom1, clone.getAtom(0));
        Assert.assertNotSame(atom2, clone.getAtom(1));
    }

    @Test
    public void testClone_Order() throws Exception {
    	IChemObject object = newChemObject(); 
        IAtom atom1 = object.getBuilder().newAtom("C");
        IAtom atom2 = object.getBuilder().newAtom("O");
        IBond bond = object.getBuilder().newBond(atom1, atom2, IBond.Order.SINGLE);
        IBond clone = (IBond) bond.clone();

        // test cloning of bond order
        bond.setOrder(IBond.Order.DOUBLE);
        Assert.assertEquals(IBond.Order.SINGLE, clone.getOrder());
    }

    @Test
    public void testClone_Stereo() throws Exception {
    	IChemObject object = newChemObject(); 
        IAtom atom1 = object.getBuilder().newAtom("C");
        IAtom atom2 = object.getBuilder().newAtom("O");
        IBond bond = object.getBuilder().newBond(
        	atom1, atom2, IBond.Order.SINGLE, IBond.Stereo.UP
        );
        IBond clone = (IBond) bond.clone();

        // test cloning of bond order
        bond.setStereo(IBond.Stereo.UP_INVERTED);
        Assert.assertEquals(IBond.Stereo.UP, clone.getStereo());
    }

    /**
     * Test for RFC #9
     */
    @Test
    public void testToString() {
        IBond bond = (IBond)newChemObject();
        String description = bond.toString();
        for (int i = 0; i < description.length(); i++) {
            Assert.assertTrue(description.charAt(i) != '\n');
            Assert.assertTrue(description.charAt(i) != '\r');
        }
    }

    @Test
    public void testMultiCenter1() {
    	IChemObject object = newChemObject(); 
        IAtom atom1 = object.getBuilder().newAtom("C");
        IAtom atom2 = object.getBuilder().newAtom("O");
        IAtom atom3 = object.getBuilder().newAtom("C");

        IBond bond = object.getBuilder().newBond(new IAtom[]{atom1, atom2, atom3});
        Assert.assertEquals(3, bond.getAtomCount());
        Assert.assertEquals(atom1, bond.getAtom(0));
        Assert.assertEquals(atom2, bond.getAtom(1));
        Assert.assertEquals(atom3, bond.getAtom(2));

        Assert.assertEquals(bond.getOrder(), CDKConstants.UNSET);
    }

    @Test
    public void testMultiCenterCompare() {
    	IChemObject object = newChemObject(); 
        IAtom atom1 = object.getBuilder().newAtom("C");
        IAtom atom2 = object.getBuilder().newAtom("O");
        IAtom atom3 = object.getBuilder().newAtom("C");

        IBond bond1 = object.getBuilder().newBond(new IAtom[]{atom1, atom2, atom3});
        IBond bond2 = object.getBuilder().newBond(new IAtom[]{atom1, atom2, atom3});

        Assert.assertTrue(bond1.compare(bond2));

        IAtom atom4 = object.getBuilder().newAtom("C");
        IBond bond3 = object.getBuilder().newBond(new IAtom[]{atom1, atom2, atom4});
        Assert.assertFalse(bond1.compare(bond3));
    }

    @Test
    public void testMltiCenterContains() {
    	IChemObject object = newChemObject(); 
        IAtom atom1 = object.getBuilder().newAtom("C");
        IAtom atom2 = object.getBuilder().newAtom("O");
        IAtom atom3 = object.getBuilder().newAtom("C");
        IAtom atom4 = object.getBuilder().newAtom("C");

        IBond bond1 = object.getBuilder().newBond(new IAtom[]{atom1, atom2, atom3});
        Assert.assertTrue(bond1.contains(atom1));
        Assert.assertTrue(bond1.contains(atom2));
        Assert.assertTrue(bond1.contains(atom3));
        Assert.assertFalse(bond1.contains(atom4));
    }

    @Test
    public void testMultiCenterIterator() {
    	IChemObject object = newChemObject(); 
        IAtom atom1 = object.getBuilder().newAtom("C");
        IAtom atom2 = object.getBuilder().newAtom("O");
        IAtom atom3 = object.getBuilder().newAtom("C");
        IAtom atom4 = object.getBuilder().newAtom("C");

        IBond bond1 = object.getBuilder().newBond(new IAtom[]{atom1, atom2, atom3, atom4});
        Iterator<IAtom> atoms = bond1.atoms().iterator();
        int natom = 0;
        while (atoms.hasNext()) {
            IAtom atom = atoms.next();
            Assert.assertNotNull(atom);
            natom++;
        }
        Assert.assertEquals(4, natom);
    }

    @Test
    public void testMultiCenterConnectedAtoms() {
    	IChemObject object = newChemObject(); 
        IAtom atom1 = object.getBuilder().newAtom("C");
        IAtom atom2 = object.getBuilder().newAtom("O");
        IAtom atom3 = object.getBuilder().newAtom("C");
        IAtom atom4 = object.getBuilder().newAtom("C");

        IBond bond1 = object.getBuilder().newBond(new IAtom[]{atom1, atom2, atom3, atom4});
        Assert.assertEquals(atom2, bond1.getConnectedAtom(atom1));
        Assert.assertNull(bond1.getConnectedAtom(object.getBuilder().newAtom()));

        IAtom[] conAtoms = bond1.getConnectedAtoms(atom1);
        boolean correct = true;
        for (IAtom atom : conAtoms) {
            if (atom == atom1) {
                correct = false;
                break;
            }
        }
        Assert.assertTrue(correct);


        conAtoms = bond1.getConnectedAtoms(atom3);
        correct = true;
        for (IAtom atom : conAtoms) {
            if (atom == atom3) {
                correct = false;
                break;
            }
        }
        Assert.assertTrue(correct);
    }

    @Test
    public void testMultiCenterIsConnectedTo() {
    	IChemObject object = newChemObject(); 
        IAtom atom1 = object.getBuilder().newAtom("C");
        IAtom atom2 = object.getBuilder().newAtom("O");
        IAtom atom3 = object.getBuilder().newAtom("C");
        IAtom atom4 = object.getBuilder().newAtom("C");
        IAtom atom5 = object.getBuilder().newAtom("C");

        IBond bond1 = object.getBuilder().newBond(new IAtom[]{atom1, atom2, atom3});
        IBond bond2 = object.getBuilder().newBond(new IAtom[]{atom2, atom3, atom4});
        IBond bond3 = object.getBuilder().newBond(new IAtom[]{atom2, atom4});
        IBond bond4 = object.getBuilder().newBond(new IAtom[]{atom5, atom4});

        Assert.assertTrue(bond1.isConnectedTo(bond2));
        Assert.assertTrue(bond2.isConnectedTo(bond1));
        Assert.assertTrue(bond1.isConnectedTo(bond3));
        Assert.assertFalse(bond4.isConnectedTo(bond1));
    }
}

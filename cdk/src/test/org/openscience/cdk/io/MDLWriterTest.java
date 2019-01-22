/* $Revision$ $Author$ $Date$
 *
 * Copyright (C) 2004-2007  The Chemistry Development Kit (CDK) project
 * 
 * Contact: cdk-devel@slists.sourceforge.net
 * 
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 * All we ask is that proper credit is given for our work, which includes
 * - but is not limited to - adding the above copyright notice to the beginning
 * of your source code files, and to any copyright notice that you may distribute
 * with programs based on this work.
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
package org.openscience.cdk.io;

import java.io.StringWriter;
import java.util.Properties;

import javax.vecmath.Point2d;
import javax.vecmath.Point3d;

import org.junit.Assert;
import org.junit.BeforeClass;
import org.junit.Test;
import org.openscience.cdk.Atom;
import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.ChemFile;
import org.openscience.cdk.ChemModel;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.Molecule;
import org.openscience.cdk.PseudoAtom;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.interfaces.IMolecule;
import org.openscience.cdk.io.listener.PropertiesListener;
import org.openscience.cdk.templates.MoleculeFactory;

/**
 * TestCase for the writer MDL mol files using one test file.
 *
 * @cdk.module test-io
 *
 * @see org.openscience.cdk.io.MDLWriter
 */
public class MDLWriterTest extends ChemObjectIOTest {

    private static IChemObjectBuilder builder;

    @BeforeClass public static void setup() {
        builder = DefaultChemObjectBuilder.getInstance();
        setChemObjectIO(new MDLWriter());
    }

    @Test public void testAccepts() throws Exception {
    	MDLWriter reader = new MDLWriter();
    	Assert.assertTrue(reader.accepts(ChemFile.class));
    	Assert.assertTrue(reader.accepts(ChemModel.class));
    	Assert.assertTrue(reader.accepts(Molecule.class));
        Assert.assertTrue(reader.accepts(AtomContainer.class));
    }

    /**
     * @cdk.bug 890456
     * @cdk.bug 1524466
     */
    @Test public void testBug890456() throws Exception {
        StringWriter writer = new StringWriter();
        Molecule molecule = new Molecule();
        molecule.addAtom(new PseudoAtom("*"));
        molecule.addAtom(new Atom("C"));
        molecule.addAtom(new Atom("C"));
        
        MDLWriter mdlWriter = new MDLWriter(writer);
        mdlWriter.write(molecule);
        Assert.assertTrue(writer.toString().indexOf("M  END") != -1);
    }

    /**
     * @cdk.bug 1212219
     */
    @Test public void testBug1212219() throws Exception {
        StringWriter writer = new StringWriter();
        Molecule molecule = new Molecule();
        Atom atom = new Atom("C");
        atom.setMassNumber(14);
        molecule.addAtom(atom);
        
        MDLWriter mdlWriter = new MDLWriter(writer);
        mdlWriter.write(molecule);
        String output = writer.toString();
        //logger.debug("MDL output for testBug1212219: " + output);
        Assert.assertTrue(output.indexOf("M  ISO  1   1  14") != -1);
    }
    
    @Test public void testWriteValence() throws Exception {
    	StringWriter writer = new StringWriter();
        Molecule molecule = MoleculeFactory.makeAlphaPinene();
        molecule.getAtom(0).setValency(1);
        molecule.getAtom(1).setValency(0);
        MDLWriter mdlWriter = new MDLWriter(writer);
        mdlWriter.write(molecule);
        String output = writer.toString();
        Assert.assertTrue(output.indexOf("0  0  0  0  0  1  0  0  0  0  0  0") != -1);
        Assert.assertTrue(output.indexOf("0  0  0  0  0 15  0  0  0  0  0  0") != -1);
    }
    
    @Test public void testWriteAtomAtomMapping() throws Exception {
        StringWriter writer = new StringWriter();
        Molecule molecule = MoleculeFactory.makeAlphaPinene();
        molecule.getAtom(0).setProperty(CDKConstants.ATOM_ATOM_MAPPING,1);
        molecule.getAtom(1).setProperty(CDKConstants.ATOM_ATOM_MAPPING,15);
        MDLWriter mdlWriter = new MDLWriter(writer);
        mdlWriter.write(molecule);
        String output = writer.toString();
        Assert.assertTrue(output.indexOf("0  0  0  0  0  0  0  0  0  1  0  0") != -1);
        Assert.assertTrue(output.indexOf("0  0  0  0  0  0  0  0  0 15  0  0") != -1);
    }
    
    /**
     * Test for bug #1778479 "MDLWriter writes empty PseudoAtom label string".
     * When a molecule contains an IPseudoAtom without specifying the atom label
     * the MDLWriter generates invalid output as it prints the zero-length atom
     * label.
     * This was fixed with letting PseudoAtom have a default label of '*'.
     *
     * Author: Andreas Schueller <a.schueller@chemie.uni-frankfurt.de>
     * 
     * @cdk.bug 1778479
     */
    @Test public void testBug1778479() throws Exception {
        StringWriter writer = new StringWriter();
        IMolecule molecule = builder.newMolecule();
        IAtom atom1 = builder.newPseudoAtom();
        IAtom atom2 = builder.newAtom("C");
        IBond bond = builder.newBond(atom1, atom2);
        molecule.addAtom(atom1);
        molecule.addAtom(atom2);
        molecule.addBond(bond);
            
        MDLWriter mdlWriter = new MDLWriter(writer);
        mdlWriter.write(molecule);
        String output = writer.toString();
        Assert.assertEquals("Test for zero length pseudo atom label in MDL file", -1, output.indexOf("0.0000    0.0000    0.0000     0  0  0  0  0  0  0  0  0  0  0  0"));
    }

    @Test public void testNullFormalCharge() throws Exception {
        StringWriter writer = new StringWriter();
        IMolecule molecule = builder.newMolecule();
        IAtom atom = builder.newAtom("C");
        atom.setFormalCharge(null);
        molecule.addAtom(atom);

        MDLWriter mdlWriter = new MDLWriter(writer);
        mdlWriter.write(molecule);
        String output = writer.toString();
        // test ensures that the writer does not throw an exception on
        // null formal charges, so a mere assert on output being non-zero
        // length is enough
        Assert.assertNotNull(output);
        Assert.assertNotSame(0, output.length());
    }

    @Test public void testPrefer3DCoordinateOutput() throws Exception {
        StringWriter writer = new StringWriter();
        IMolecule molecule = builder.newMolecule();
        IAtom atom = builder.newAtom("C");
        atom.setPoint2d(new Point2d(1.0, 2.0));
        atom.setPoint3d(new Point3d(3.0, 4.0, 5.0));
        molecule.addAtom(atom);

        MDLWriter mdlWriter = new MDLWriter(writer);
        mdlWriter.write(molecule);
        mdlWriter.close();
        String output = writer.toString();
        // the current behavior is that if both 2D and 3D coordinates
        // are available, the 3D is outputed, and the 2D not
        Assert.assertTrue(output.contains("3.0"));
        Assert.assertTrue(output.contains("4.0"));
        Assert.assertTrue(output.contains("5.0"));
    }

    @Test public void testForce2DCoordinates() throws Exception {
        StringWriter writer = new StringWriter();
        IMolecule molecule = builder.newMolecule();
        IAtom atom = builder.newAtom("C");
        atom.setPoint2d(new Point2d(1.0, 2.0));
        atom.setPoint3d(new Point3d(3.0, 4.0, 5.0));
        molecule.addAtom(atom);

        MDLWriter mdlWriter = new MDLWriter(writer);
        Properties prop = new Properties();
        prop.setProperty("ForceWriteAs2DCoordinates","true");
        PropertiesListener listener = new PropertiesListener(prop);
        mdlWriter.addChemObjectIOListener(listener);
        mdlWriter.customizeJob();
        mdlWriter.write(molecule);
        mdlWriter.close();
        String output = writer.toString();
        // the current behavior is that if both 2D and 3D coordinates
        // are available, the 3D is outputed, and the 2D not
        Assert.assertTrue(output.contains("1.0"));
        Assert.assertTrue(output.contains("2.0"));
    }

    @Test public void testUndefinedStereo() throws Exception {
      IMolecule mol = MoleculeFactory.makeAlphaPinene();
      mol.getBond(0).setStereo(IBond.Stereo.UP_OR_DOWN);
      mol.getBond(1).setStereo(IBond.Stereo.E_OR_Z);
      StringWriter writer = new StringWriter();
        MDLWriter mdlWriter = new MDLWriter(writer);
        mdlWriter.write(mol);
        String output = writer.toString();
        Assert.assertTrue(output.indexOf("1  2  2  4  0  0  0")>-1);
        Assert.assertTrue(output.indexOf("2  3  1  3  0  0  0")>-1);
    }


}

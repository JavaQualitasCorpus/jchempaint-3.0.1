/* $RCSfile$    
 * $Author$    
 * $Date$    
 * $Revision$
 * 
 * Copyright (C) 1997-2007  The Chemistry Development Kit (CKD) project
 *                    2009  Egon Willighagen <egonw@users.sf.net>
 * 
 * Contact: cdk-devel@lists.sourceforge.net
 * 
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 * All I ask is that proper credit is given for my work, which includes
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
package org.openscience.cdk.fingerprint;

import java.util.BitSet;

import org.junit.Assert;
import org.junit.Test;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.smiles.SmilesParser;

/**
 * @cdk.module test-fingerprint
 */
public class SubstructureFingerprinterTest extends AbstractFingerprinterTest {

    public IFingerprinter getFingerprinter() {
        return new SubstructureFingerprinter();
    }

    @Test
    public void testSize() throws Exception {
        SubstructureFingerprinter fp = new SubstructureFingerprinter();
        Assert.assertEquals(307, fp.getSize());
    }
    
    @Test public void testUserFunctionalGroups() throws Exception {
        String[] smarts = {"c1ccccc1", "[CX4H3][#6]", "[CX2]#[CX2]"};
        IFingerprinter printer = new SubstructureFingerprinter(smarts);
        Assert.assertEquals(3, printer.getSize());

        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer mol1 = sp.parseSmiles("c1ccccc1CCC");
        BitSet fp = printer.getFingerprint(mol1);
		Assert.assertNotNull(fp);

        Assert.assertTrue(fp.get(0));
        Assert.assertTrue(fp.get(1));
        Assert.assertFalse(fp.get(2));

        mol1 = sp.parseSmiles("C=C=C");
        fp = printer.getFingerprint(mol1);
		Assert.assertNotNull(fp);
        Assert.assertFalse(fp.get(0));
        Assert.assertFalse(fp.get(1));
        Assert.assertFalse(fp.get(2));
    }

    @Test
    public void testFingerprint() throws Exception {
        IFingerprinter printer = new SubstructureFingerprinter();
        Assert.assertEquals(307, printer.getSize());

        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer mol1 = sp.parseSmiles("c1ccccc1CCC");
        BitSet fp = printer.getFingerprint(mol1);
        Assert.assertNotNull(fp);
        Assert.assertTrue(fp.get(273));
        Assert.assertTrue(fp.get(0));
        Assert.assertTrue(fp.get(1));
        Assert.assertFalse(fp.get(100));
    }
}


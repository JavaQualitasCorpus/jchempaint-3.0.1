/* $Revision: 5889 $ $Author: egonw $ $Date: 2006-04-06 15:24:58 +0200 (Thu, 06 Apr 2006) $
 * 
 * Copyright (C) 2006-2007  Egon Willighagen <egonw@users.sf.net>
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
package org.openscience.cdk.atomtype;

import java.util.List;

import org.junit.Assert;
import org.junit.Test;
import org.openscience.cdk.Atom;
import org.openscience.cdk.Bond;
import org.openscience.cdk.Molecule;
import org.openscience.cdk.CDKTestCase;
import org.openscience.cdk.interfaces.IAtomType;
import org.openscience.cdk.interfaces.IBond;

/**
 * @cdk.module test-structgen
 */
public class StructGenAtomTypeGuesserTest extends CDKTestCase {

    @Test
    public void testPossibleAtomTypes_IAtomContainer_IAtom() throws java.lang.Exception {
        Molecule mol = new Molecule();
        Atom atom = new Atom("C");
        atom.setHydrogenCount(3);
        Atom atom2 = new Atom("N");
        atom2.setHydrogenCount(2);
        mol.addAtom(atom);
        mol.addAtom(atom2);
        mol.addBond(new Bond(atom, atom2, IBond.Order.SINGLE));

        StructGenAtomTypeGuesser atm = new StructGenAtomTypeGuesser();
        List<IAtomType> matched = atm.possibleAtomTypes(mol, atom);
        Assert.assertNotNull(matched);
        Assert.assertTrue(matched.size() > 0);
        Assert.assertTrue(matched.get(0) instanceof IAtomType);

        Assert.assertEquals("C", ((IAtomType)matched.get(0)).getSymbol());
    }

    @Test
    public void testStructGenAtomTypeGuesser() throws Exception {
        StructGenAtomTypeGuesser matcher = new StructGenAtomTypeGuesser();
        Assert.assertNotNull(matcher);
    }
}
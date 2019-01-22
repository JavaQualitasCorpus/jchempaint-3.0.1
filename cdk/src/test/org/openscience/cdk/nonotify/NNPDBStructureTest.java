/* $RCSfile$
 * $Author$    
 * $Date$    
 * $Revision$
 * 
 *  Copyright (C) 2004-2007  Miguel Rojas <miguel.rojas@uni-koeln.de>
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
package org.openscience.cdk.nonotify;

import org.junit.Assert;
import org.junit.BeforeClass;
import org.junit.Test;
import org.openscience.cdk.interfaces.IPDBStructure;
import org.openscience.cdk.interfaces.AbstractPDBStructureTest;

/**
 * Checks the functionality of the {@link NNPDBStructure}.
 *
 * @cdk.module test-nonotify
 */
public class NNPDBStructureTest extends AbstractPDBStructureTest {

    @BeforeClass public static void setUp() {
        setChemObject(new NNPDBStructure());
    }

	@Test public void testNNPDBStructure() {
		IPDBStructure structure = new NNPDBStructure();
		Assert.assertNotNull(structure);
	}

    @Test public void testGetBuilder() {
        NNPDBStructure structure = new NNPDBStructure();
        Assert.assertTrue(structure.getBuilder() instanceof NoNotificationChemObjectBuilder);
    }

    @Test public void testAddListener_IChemObjectListener() {
        NNChemObjectTestHelper.testAddListener_IChemObjectListener(newChemObject());
    }

}

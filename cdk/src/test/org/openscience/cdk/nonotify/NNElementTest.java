/* $RCSfile$
 * $Author$    
 * $Date$    
 * $Revision$
 * 
 * Copyright (C) 1997-2007  The Chemistry Development Kit (CDK) project
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
import org.openscience.cdk.interfaces.IChemObject;
import org.openscience.cdk.interfaces.IElement;
import org.openscience.cdk.interfaces.AbstractElementTest;
import org.openscience.cdk.interfaces.ITestObjectBuilder;

/**
 * Checks the functionality of {@link NNElement}.
 *
 * @cdk.module test-nonotify
 */
public class NNElementTest extends AbstractElementTest {

    @BeforeClass public static void setUp() {
        setTestObjectBuilder(new ITestObjectBuilder() {
            public IChemObject newTestObject() {
                return new NNElement();
            }
        });
    }

    @Test public void testNNElement() {
        IElement e = new NNElement();
        Assert.assertTrue(e instanceof IChemObject);
    }
    
    @Test public void testNNElement_IElement() {
    	IElement element = new NNElement();
        IElement e = element.getBuilder().newElement(element);
        Assert.assertTrue(e instanceof IChemObject);
    }
    
    @Test public void testNNElement_String() {
        IElement e = new NNElement("C");
        Assert.assertEquals("C", e.getSymbol());
    }
    
    @Test public void testNNElement_String_int() {
        IElement e = new NNElement("H", 1);
        Assert.assertEquals("H", e.getSymbol());
        Assert.assertEquals(1, e.getAtomicNumber().intValue());
    }

    // Overwrite default methods: no notifications are expected!
    
    @Test public void testNotifyChanged() {
        NNChemObjectTestHelper.testNotifyChanged(newChemObject());
    }
    @Test public void testNotifyChanged_IChemObjectChangeEvent() {
        NNChemObjectTestHelper.testNotifyChanged_IChemObjectChangeEvent(newChemObject());
    }
    @Test public void testStateChanged_IChemObjectChangeEvent() {
        NNChemObjectTestHelper.testStateChanged_IChemObjectChangeEvent(newChemObject());
    }
    @Test public void testClone_ChemObjectListeners() throws Exception {
        NNChemObjectTestHelper.testClone_ChemObjectListeners(newChemObject());
    }
    @Test public void testAddListener_IChemObjectListener() {
        NNChemObjectTestHelper.testAddListener_IChemObjectListener(newChemObject());
    }
    @Test public void testGetListenerCount() {
        NNChemObjectTestHelper.testGetListenerCount(newChemObject());
    }
    @Test public void testRemoveListener_IChemObjectListener() {
        NNChemObjectTestHelper.testRemoveListener_IChemObjectListener(newChemObject());
    }
    @Test public void testSetNotification_true() {
        NNChemObjectTestHelper.testSetNotification_true(newChemObject());
    }
}

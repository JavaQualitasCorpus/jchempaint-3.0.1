/* $Revision$ $Author$ $Date$    
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
import org.openscience.cdk.interfaces.IChemObject;

/**
 * Helper class to test the functionality of the NNChemObject.
 *
 * @cdk.module test-nonotify
 */
public class NNChemObjectTestHelper {

    public static void testNotifyChanged(IChemObject chemObject) {
        NNChemObjectListener listener = new NNChemObjectListener();
        chemObject.addListener(listener);
        
        chemObject.setID("Changed");
        Assert.assertFalse(listener.getChanged());
    }

    public static void testNotifyChanged_IChemObjectChangeEvent(IChemObject chemObject) {
        NNChemObjectListener listener = new NNChemObjectListener();
        chemObject.addListener(listener);
        
        chemObject.setID("Changed");
        Assert.assertNull(listener.getEvent());
    }

    public static void testStateChanged_IChemObjectChangeEvent(IChemObject chemObject) {
        NNChemObjectListener listener = new NNChemObjectListener();
        chemObject.addListener(listener);
        
        chemObject.setID("Changed");
        Assert.assertFalse(listener.getChanged());
        
        listener.reset();
        Assert.assertFalse(listener.getChanged());
        chemObject.setProperty("Changed", "Again");
        Assert.assertFalse(listener.getChanged());

        listener.reset();
        Assert.assertFalse(listener.getChanged());
        chemObject.setFlag(3, true);
        Assert.assertFalse(listener.getChanged());
    }

    public static void testClone_ChemObjectListeners(IChemObject chemObject) throws Exception {
        NNChemObjectListener listener = new NNChemObjectListener();
        chemObject.addListener(listener);
        IChemObject chemObject2 = (IChemObject)chemObject.clone();

        // test lack of cloning of listeners
        Assert.assertEquals(0, chemObject.getListenerCount());
        Assert.assertEquals(0, chemObject2.getListenerCount());
    }
    
    public static void testAddListener_IChemObjectListener(IChemObject chemObject) {
        Assert.assertEquals(0, chemObject.getListenerCount());
        NNChemObjectListener listener = new NNChemObjectListener();
        chemObject.addListener(listener);
        Assert.assertEquals(0, chemObject.getListenerCount());
    }
    
    public static void testGetListenerCount(IChemObject chemObject) {
        NNChemObjectListener listener = new NNChemObjectListener();
        chemObject.addListener(listener);
        Assert.assertEquals(0, chemObject.getListenerCount());
    }

    public static void testRemoveListener_IChemObjectListener(IChemObject chemObject) {
        NNChemObjectListener listener = new NNChemObjectListener();
        chemObject.addListener(listener);
        Assert.assertEquals(0, chemObject.getListenerCount());
        chemObject.removeListener(listener);
        Assert.assertEquals(0, chemObject.getListenerCount());
    }
    
    public static void testSetNotification_true(IChemObject chemObject) {
        NNChemObjectListener listener = new NNChemObjectListener();
        chemObject.addListener(listener);
        chemObject.setNotification(true);
        
        chemObject.setID("Changed");
        Assert.assertFalse(listener.getChanged());
    }
}

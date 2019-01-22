/* $Revision$ $Author$ $Date$
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
 */
package org.openscience.cdk.interfaces;

import javax.vecmath.Point2d;
import javax.vecmath.Point3d;

import org.junit.Assert;
import org.junit.Test;
import org.openscience.cdk.CDKTestCase;

/**
 * Checks the functionality of {@link IChemObjectBuilder} implementations.
 *
 * @cdk.module test-interfaces
 */
public abstract class AbstractChemObjectBuilderTest extends CDKTestCase {

	private static IChemObject rootObject;

	public static IChemObject getRootObject() {
		return rootObject;
	}

	public static void setRootObject(IChemObject rootObject) {
		AbstractChemObjectBuilderTest.rootObject = rootObject;
	}

	@Test public void testNewAminoAcid() {
		Object object = rootObject.getBuilder().newAminoAcid();
		Assert.assertNotNull(object);
		Assert.assertTrue(object instanceof IChemObject);
		Assert.assertTrue(object instanceof IAminoAcid);
	}
	 
	@Test public void testNewAtom() {
		Object object = rootObject.getBuilder().newAtom();
		Assert.assertNotNull(object);
		Assert.assertTrue(object instanceof IChemObject);
		Assert.assertTrue(object instanceof IAtom);
	}
	 
	@Test public void testNewAtom_IElement() {
		IElement element = rootObject.getBuilder().newElement();
		Object object = rootObject.getBuilder().newAtom(element);
		Assert.assertNotNull(object);
		Assert.assertTrue(object instanceof IChemObject);
		Assert.assertTrue(object instanceof IAtom);
	}
	 
	@Test public void testNewAtom_String() {
		Object object = rootObject.getBuilder().newAtom("C");
		Assert.assertNotNull(object);
		Assert.assertTrue(object instanceof IChemObject);
		Assert.assertTrue(object instanceof IAtom);
	}
	
	@Test public void testNewAtom_String_Point2d() {
		Object object = rootObject.getBuilder().newAtom("C", new Point2d(1.0, 2.0));
		Assert.assertNotNull(object);
		Assert.assertTrue(object instanceof IChemObject);
		Assert.assertTrue(object instanceof IAtom);
	}
	
	@Test public void testNewAtom_String_Point3d() {
		Object object = rootObject.getBuilder().newAtom("C", new Point3d(1.0, 2.0, 3.0));
		Assert.assertNotNull(object);
		Assert.assertTrue(object instanceof IChemObject);
		Assert.assertTrue(object instanceof IAtom);
	}
	
	@Test public void testNewAtomContainer() {
		Object object = rootObject.getBuilder().newAtomContainer();
		Assert.assertNotNull(object);
		Assert.assertTrue(object instanceof IChemObject);
		Assert.assertTrue(object instanceof IAtomContainer);
	}
	
	@Test public void testNewAtomContainer_int_int_int_int() {
		Object object = rootObject.getBuilder().newAtomContainer(10,10,0,0);
		Assert.assertNotNull(object);
		Assert.assertTrue(object instanceof IChemObject);
		Assert.assertTrue(object instanceof IAtomContainer);
	}
	
	@Test public void testNewAtomContainer_IAtomContainer() {
		Object object = rootObject.getBuilder().newAtomContainer(rootObject.getBuilder().newAtomContainer());
		Assert.assertNotNull(object);
		Assert.assertTrue(object instanceof IChemObject);
		Assert.assertTrue(object instanceof IAtomContainer);
	}
	
	@Test public void testNewAtomParity_IAtom_IAtom_IAtom_IAtom_IAtom_int() {
		Object object = rootObject.getBuilder().newAtomParity(
				rootObject.getBuilder().newAtom(),
				rootObject.getBuilder().newAtom(),
				rootObject.getBuilder().newAtom(),
				rootObject.getBuilder().newAtom(),
				rootObject.getBuilder().newAtom(),
				1
		);
		Assert.assertNotNull(object);
		Assert.assertTrue(object instanceof IAtomParity);
	}
	
	@Test public void testNewAtomType_String() {
		Object object = rootObject.getBuilder().newAtomType("Carom");
		Assert.assertNotNull(object);
		Assert.assertTrue(object instanceof IChemObject);
		Assert.assertTrue(object instanceof IAtomType);
	}
	
	@Test public void testNewAtomType_IElement() {
		IElement element = rootObject.getBuilder().newElement();
		Object object = rootObject.getBuilder().newAtomType(element);
		Assert.assertNotNull(object);
		Assert.assertTrue(object instanceof IChemObject);
		Assert.assertTrue(object instanceof IAtomType);
	}
	
	@Test public void testNewAtomType_String_String() {
		Object object = rootObject.getBuilder().newAtomType("Carom");
		Assert.assertNotNull(object);
		Assert.assertTrue(object instanceof IChemObject);
		Assert.assertTrue(object instanceof IAtomType);
	}
	
	@Test public void testNewBioPolymer() {
		Object object = rootObject.getBuilder().newBioPolymer();
		Assert.assertNotNull(object);
		Assert.assertTrue(object instanceof IChemObject);
		Assert.assertTrue(object instanceof IBioPolymer);
	}
	
	@Test public void testNewBond() {
		Object object = rootObject.getBuilder().newBond();
		Assert.assertNotNull(object);
		Assert.assertTrue(object instanceof IChemObject);
		Assert.assertTrue(object instanceof IBond);
	}
	
	@Test public void testNewBond_IAtom_IAtom() {
		Object object = rootObject.getBuilder().newBond(
			rootObject.getBuilder().newAtom(),
			rootObject.getBuilder().newAtom()
		);
		Assert.assertNotNull(object);
		Assert.assertTrue(object instanceof IChemObject);
		Assert.assertTrue(object instanceof IBond);
	}
	
	@Test public void testNewBond_IAtom_IAtom_IBond_Order() {
		Object object = rootObject.getBuilder().newBond(
			rootObject.getBuilder().newAtom(),
			rootObject.getBuilder().newAtom(),
			IBond.Order.SINGLE
		);
		Assert.assertNotNull(object);
		Assert.assertTrue(object instanceof IChemObject);
		Assert.assertTrue(object instanceof IBond);
	}
	
	@Test public void testNewBond_IAtom_IAtom_IBond_Order_IBond_Stereo() {
		Object object = rootObject.getBuilder().newBond(
			rootObject.getBuilder().newAtom(),
			rootObject.getBuilder().newAtom(),
			IBond.Order.SINGLE,
			IBond.Stereo.UP
		);
		Assert.assertNotNull(object);
		Assert.assertTrue(object instanceof IChemObject);
		Assert.assertTrue(object instanceof IBond);
	}

    @Test
    public void testNewBond_arrayIAtom() {
        Object object = rootObject.getBuilder().newBond(
                new IAtom[]{
                        rootObject.getBuilder().newAtom(),
                        rootObject.getBuilder().newAtom(),
                        rootObject.getBuilder().newAtom(),
                        rootObject.getBuilder().newAtom()}
        );
        Assert.assertNotNull(object);
        Assert.assertTrue(object instanceof IChemObject);
        Assert.assertTrue(object instanceof IBond);
    }

     @Test
    public void testNewBond_arrayIAtom_IBond_Order() {
        Object object = rootObject.getBuilder().newBond(
                new IAtom[]{
                        rootObject.getBuilder().newAtom(),
                        rootObject.getBuilder().newAtom(),
                        rootObject.getBuilder().newAtom(),
                        rootObject.getBuilder().newAtom()},
                IBond.Order.SINGLE
        );
        Assert.assertNotNull(object);
        Assert.assertTrue(object instanceof IChemObject);
        Assert.assertTrue(object instanceof IBond);
    }


    @Test public void testNewChemFile() {
		Object object = rootObject.getBuilder().newChemFile();
		Assert.assertNotNull(object);
		Assert.assertTrue(object instanceof IChemObject);
		Assert.assertTrue(object instanceof IChemFile);
	}
	
	@Test public void testNewChemModel() {
		Object object = rootObject.getBuilder().newChemModel();
		Assert.assertNotNull(object);
		Assert.assertTrue(object instanceof IChemObject);
		Assert.assertTrue(object instanceof IChemModel);
	}

	@Test public void testNewChemObject() {
		Object object = rootObject.getBuilder().newChemObject();
		Assert.assertNotNull(object);
		Assert.assertTrue(object instanceof IChemObject);
	}

	@Test public void testNewChemObject_IChemObject() {
		IChemObject chemObject = rootObject.getBuilder().newChemObject();
		Object object = rootObject.getBuilder().newChemObject(chemObject);
		Assert.assertNotNull(object);
		Assert.assertTrue(object instanceof IChemObject);
	}

	@Test public void testNewChemSequence() {
		Object object = rootObject.getBuilder().newChemSequence();
		Assert.assertNotNull(object);
		Assert.assertTrue(object instanceof IChemObject);
		Assert.assertTrue(object instanceof IChemSequence);
	}
	
	@Test public void testNewCrystal() {
		Object object = rootObject.getBuilder().newCrystal();
		Assert.assertNotNull(object);
		Assert.assertTrue(object instanceof IChemObject);
		Assert.assertTrue(object instanceof ICrystal);
	}
	
	@Test public void testNewCrystal_IAtomContainer() {
		Object object = rootObject.getBuilder().newCrystal(
			rootObject.getBuilder().newAtomContainer()
		);
		Assert.assertNotNull(object);
		Assert.assertTrue(object instanceof IChemObject);
		Assert.assertTrue(object instanceof ICrystal);
	}
	
	@Test public void testNewElectronContainer() {
		Object object = rootObject.getBuilder().newElectronContainer();
		Assert.assertNotNull(object);
		Assert.assertTrue(object instanceof IChemObject);
		Assert.assertTrue(object instanceof IElectronContainer);
	}

	@Test public void testNewElement() {
		Object object = rootObject.getBuilder().newElement();
		Assert.assertNotNull(object);
		Assert.assertTrue(object instanceof IChemObject);
		Assert.assertTrue(object instanceof IElement);
	}
	
	@Test public void testNewElement_IElement() {
		IElement element = rootObject.getBuilder().newElement();
		Object object = rootObject.getBuilder().newElement(element);
		Assert.assertNotNull(object);
		Assert.assertTrue(object instanceof IChemObject);
		Assert.assertTrue(object instanceof IElement);
	}
	
	@Test public void testNewElement_String() {
		Object object = rootObject.getBuilder().newElement("C");
		Assert.assertNotNull(object);
		Assert.assertTrue(object instanceof IChemObject);
		Assert.assertTrue(object instanceof IElement);
	}
	
	@Test public void testNewElement_String_int() {
		Object object = rootObject.getBuilder().newElement("C", 6);
		Assert.assertNotNull(object);
		Assert.assertTrue(object instanceof IChemObject);
		Assert.assertTrue(object instanceof IElement);
	}

	@Test public void testNewIsotope_int_String_double_double() {
		Object object = rootObject.getBuilder().newIsotope(
			12, "C", 12.001, 100.0
		);
		Assert.assertNotNull(object);
		Assert.assertTrue(object instanceof IChemObject);
		Assert.assertTrue(object instanceof IIsotope);
	}
	
	@Test public void testNewIsotope_int_String_int_double_double() {
		Object object = rootObject.getBuilder().newIsotope(
			12, "C", 6, 12.001, 100.0
		);
		Assert.assertNotNull(object);
		Assert.assertTrue(object instanceof IChemObject);
		Assert.assertTrue(object instanceof IIsotope);
	}
	
	@Test public void testNewIsotope_IElement() {
		IElement element = rootObject.getBuilder().newElement();
		Object object = rootObject.getBuilder().newIsotope(element);
		Assert.assertNotNull(object);
		Assert.assertTrue(object instanceof IChemObject);
		Assert.assertTrue(object instanceof IIsotope);
	}
	
	@Test public void testNewIsotope_String() {
		Object object = rootObject.getBuilder().newIsotope("N");
		Assert.assertNotNull(object);
		Assert.assertTrue(object instanceof IChemObject);
		Assert.assertTrue(object instanceof IIsotope);
	}
	
	@Test public void testNewIsotope_String_int() {
		Object object = rootObject.getBuilder().newIsotope("N", 5);
		Assert.assertNotNull(object);
		Assert.assertTrue(object instanceof IChemObject);
		Assert.assertTrue(object instanceof IIsotope);
	}

	@Test public void testNewLonePair() {
		Object object = rootObject.getBuilder().newLonePair();
		Assert.assertNotNull(object);
		Assert.assertTrue(object instanceof IChemObject);
		Assert.assertTrue(object instanceof ILonePair);
	}	

	@Test public void testNewLonePair_IAtom() {
		Object object = rootObject.getBuilder().newLonePair(
			rootObject.getBuilder().newAtom()
		);
		Assert.assertNotNull(object);
		Assert.assertTrue(object instanceof IChemObject);
		Assert.assertTrue(object instanceof ILonePair);
	}	

    @Test public void testNewMapping_IChemObject_IChemObject() {
		Object object = rootObject.getBuilder().newMapping(
			rootObject.getBuilder().newAtom(), 
            rootObject.getBuilder().newAtom()
        );
		Assert.assertNotNull(object);
		Assert.assertTrue(object instanceof IMapping);
	}
    
	@Test public void testNewMolecule() {
		Object object = rootObject.getBuilder().newMolecule();
		Assert.assertNotNull(object);
		Assert.assertTrue(object instanceof IChemObject);
		Assert.assertTrue(object instanceof IMolecule);
	}	

	@Test public void testNewMolecule_int_int_int_int() {
		Object object = rootObject.getBuilder().newMolecule(5,5,1,1);
		Assert.assertNotNull(object);
		Assert.assertTrue(object instanceof IChemObject);
		Assert.assertTrue(object instanceof IMolecule);
	}	

	@Test public void testNewMolecule_IAtomContainer() {
		Object object = rootObject.getBuilder().newMolecule(
			rootObject.getBuilder().newAtomContainer()
		);
		Assert.assertNotNull(object);
		Assert.assertTrue(object instanceof IChemObject);
		Assert.assertTrue(object instanceof IMolecule);
	}	

	@Test public void testNewMonomer() {
		Object object = rootObject.getBuilder().newMonomer();
		Assert.assertNotNull(object);
		Assert.assertTrue(object instanceof IChemObject);
		Assert.assertTrue(object instanceof IMonomer);
	}	

	@Test public void testNewPolymer() {
		Object object = rootObject.getBuilder().newPolymer();
		Assert.assertNotNull(object);
		Assert.assertTrue(object instanceof IChemObject);
		Assert.assertTrue(object instanceof IPolymer);
	}	
	
	@Test public void testNewPDBAtom_IElement() {
		IElement element = rootObject.getBuilder().newElement();
		Object object = rootObject.getBuilder().newPDBAtom(element);
		Assert.assertNotNull(object);
		Assert.assertTrue(object instanceof IChemObject);
		Assert.assertTrue(object instanceof IPDBAtom);
	}	
	
	@Test public void testNewPDBAtom_String() {
		Object object = rootObject.getBuilder().newPDBAtom("C");
		Assert.assertNotNull(object);
		Assert.assertTrue(object instanceof IChemObject);
		Assert.assertTrue(object instanceof IPDBAtom);
	}	
	
	@Test public void testNewPDBAtom_String_Point3D() {
		Object object = rootObject.getBuilder().newPDBAtom("C", new Point3d(1.0, 2.0, 3.0));
		Assert.assertNotNull(object);
		Assert.assertTrue(object instanceof IChemObject);
		Assert.assertTrue(object instanceof IPDBAtom);
	}
	
	@Test public void testNewPDBPolymer() {
		Object object = rootObject.getBuilder().newPDBPolymer();
		Assert.assertNotNull(object);
		Assert.assertTrue(object instanceof IChemObject);
		Assert.assertTrue(object instanceof IPDBPolymer);
	}
	
	@Test public void testNewPDBStructure() {
		Object object = rootObject.getBuilder().newPDBStructure();
		Assert.assertNotNull(object);
		Assert.assertTrue(object instanceof IChemObject);
		Assert.assertTrue(object instanceof IPDBStructure);
	}
	
	@Test public void testNewPDBMonomer() {
		Object object = rootObject.getBuilder().newPDBMonomer();
		Assert.assertNotNull(object);
		Assert.assertTrue(object instanceof IChemObject);
		Assert.assertTrue(object instanceof IPDBMonomer);
	}

	
	@Test public void testNewPseudoAtom() {
		Object object = rootObject.getBuilder().newPseudoAtom();
		Assert.assertNotNull(object);
		Assert.assertTrue(object instanceof IChemObject);
		Assert.assertTrue(object instanceof IPseudoAtom);
	}

	@Test public void testNewPseudoAtom_IElement() {
		IElement element = rootObject.getBuilder().newElement();
		Object object = rootObject.getBuilder().newPseudoAtom(element);
		Assert.assertNotNull(object);
		Assert.assertTrue(object instanceof IChemObject);
		Assert.assertTrue(object instanceof IPseudoAtom);
	}	

	@Test public void testNewPseudoAtom_IAtom() {
		Object object = rootObject.getBuilder().newPseudoAtom(
			rootObject.getBuilder().newAtom()
		);
		Assert.assertNotNull(object);
		Assert.assertTrue(object instanceof IChemObject);
		Assert.assertTrue(object instanceof IPseudoAtom);
	}	

	@Test public void testNewPseudoAtom_String() {
		Object object = rootObject.getBuilder().newPseudoAtom("Glu178");
		Assert.assertNotNull(object);
		Assert.assertTrue(object instanceof IChemObject);
		Assert.assertTrue(object instanceof IPseudoAtom);
	}	

	@Test public void testNewPseudoAtom_String_Point2d() {
		Object object = rootObject.getBuilder().newPseudoAtom("Glue178", new Point2d(1.0,2.0));
		Assert.assertNotNull(object);
		Assert.assertTrue(object instanceof IChemObject);
		Assert.assertTrue(object instanceof IPseudoAtom);
	}	

	@Test public void testNewPseudoAtom_String_Point3d() {
		Object object = rootObject.getBuilder().newPseudoAtom("Glue178", new Point3d(1.0,2.0,3.0));
		Assert.assertNotNull(object);
		Assert.assertTrue(object instanceof IChemObject);
		Assert.assertTrue(object instanceof IPseudoAtom);
	}	

	@Test public void testNewPDBAtom_String_Point3d() {
		Object object = rootObject.getBuilder().newPDBAtom("CA", new Point3d(1.0,2.0,3.0));
		Assert.assertNotNull(object);
		Assert.assertTrue(object instanceof IChemObject);
		Assert.assertTrue(object instanceof IPDBAtom);
	}	

	@Test public void testNewReaction() {
		Object object = rootObject.getBuilder().newReaction();
		Assert.assertNotNull(object);
		Assert.assertTrue(object instanceof IChemObject);
		Assert.assertTrue(object instanceof IReaction);
	}

	@Test public void testNewRing() {
		Object object = rootObject.getBuilder().newRing();
		Assert.assertNotNull(object);
		Assert.assertTrue(object instanceof IChemObject);
		Assert.assertTrue(object instanceof IRing);
	}	

	@Test public void testNewRing_int() {
		Object object = rootObject.getBuilder().newRing(5);
		Assert.assertNotNull(object);
		Assert.assertTrue(object instanceof IChemObject);
		Assert.assertTrue(object instanceof IRing);
	}	

	@Test public void testNewRing_int_String() {
		Object object = rootObject.getBuilder().newRing(5,"C");
		Assert.assertNotNull(object);
		Assert.assertTrue(object instanceof IChemObject);
		Assert.assertTrue(object instanceof IRing);
	}	

	@Test public void testNewRing_IAtomContainer() {
		Object object = rootObject.getBuilder().newRing(
			rootObject.getBuilder().newAtomContainer()
		);
		Assert.assertNotNull(object);
		Assert.assertTrue(object instanceof IChemObject);
		Assert.assertTrue(object instanceof IRing);
	}	

	@Test public void testNewRingSet() {
		Object object = rootObject.getBuilder().newRingSet();
		Assert.assertNotNull(object);
		Assert.assertTrue(object instanceof IRingSet);
	}

	@Test public void testNewAtomContainerSet() {
		Object object = rootObject.getBuilder().newAtomContainerSet();
		Assert.assertNotNull(object);
		Assert.assertTrue(object instanceof IChemObject);
		Assert.assertTrue(object instanceof IAtomContainerSet);
	}

	@Test public void testNewMoleculeSet() {
		Object object = rootObject.getBuilder().newMoleculeSet();
		Assert.assertNotNull(object);
		Assert.assertTrue(object instanceof IChemObject);
		Assert.assertTrue(object instanceof IMoleculeSet);
	}

	@Test public void testNewReactionSet() {
		Object object = rootObject.getBuilder().newReactionSet();
		Assert.assertNotNull(object);
		Assert.assertTrue(object instanceof IChemObject);
		Assert.assertTrue(object instanceof IReactionSet);
	}

	@Test public void testNewReactionScheme() {
		Object object = rootObject.getBuilder().newReactionScheme();
		Assert.assertNotNull(object);
		Assert.assertTrue(object instanceof IChemObject);
		Assert.assertTrue(object instanceof IReactionScheme);
	}

	@Test public void testNewSingleElectron() {
		Object object = rootObject.getBuilder().newSingleElectron();
		Assert.assertNotNull(object);
		Assert.assertTrue(object instanceof IChemObject);
		Assert.assertTrue(object instanceof ISingleElectron);
	}

	@Test public void testNewSingleElectron_IAtom() {
		Object object = rootObject.getBuilder().newSingleElectron(
			rootObject.getBuilder().newAtom()
		);
		Assert.assertNotNull(object);
		Assert.assertTrue(object instanceof IChemObject);
		Assert.assertTrue(object instanceof ISingleElectron);
	}
	
	@Test public void testNewStrand() {
		Object object = rootObject.getBuilder().newStrand();
		Assert.assertNotNull(object);
		Assert.assertTrue(object instanceof IChemObject);
		Assert.assertTrue(object instanceof IStrand);
	}
	
	@Test public void testNewFragmentAtom() {
		Object object = rootObject.getBuilder().newFragmentAtom();
		Assert.assertNotNull(object);
		Assert.assertTrue(object instanceof IChemObject);
		Assert.assertTrue(object instanceof IFragmentAtom);
	}

	@Test public void testNewMolecularFormula() {
		Object object = rootObject.getBuilder().newMolecularFormula();
		Assert.assertNotNull(object);
		Assert.assertTrue(object instanceof IMolecularFormula);
	}

	@Test public void testNewMolecularFormulaSet() {
		Object object = rootObject.getBuilder().newMolecularFormulaSet();
		Assert.assertNotNull(object);
		Assert.assertTrue(object instanceof IMolecularFormulaSet);
	}

	@Test public void testNewMolecularFormulaSet_IMolecularFormula() {
		Object object = rootObject.getBuilder().newMolecularFormulaSet(
			rootObject.getBuilder().newMolecularFormula()
		);
		Assert.assertNotNull(object);
		Assert.assertTrue(object instanceof IMolecularFormulaSet);
	}

	@Test public void testNewAdductFormula() {
		Object object = rootObject.getBuilder().newAdductFormula();
		Assert.assertNotNull(object);
		Assert.assertTrue(object instanceof IAdductFormula);
	}

	@Test public void testNewAdductFormula_IMolecularFormula() {
		Object object = rootObject.getBuilder().newAdductFormula(
			rootObject.getBuilder().newMolecularFormula()
		);
		Assert.assertNotNull(object);
		Assert.assertTrue(object instanceof IAdductFormula);
		
	}
	
}

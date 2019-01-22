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
 */
package org.openscience.cdk.reaction.type;


import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import org.junit.Assert;
import org.junit.Test;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.SingleElectron;
import org.openscience.cdk.atomtype.CDKAtomTypeMatcher;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.interfaces.IMolecule;
import org.openscience.cdk.interfaces.IMoleculeSet;
import org.openscience.cdk.interfaces.IReactionSet;
import org.openscience.cdk.isomorphism.UniversalIsomorphismTester;
import org.openscience.cdk.isomorphism.matchers.IQueryAtomContainer;
import org.openscience.cdk.isomorphism.matchers.QueryAtomContainerCreator;
import org.openscience.cdk.nonotify.NoNotificationChemObjectBuilder;
import org.openscience.cdk.reaction.IReactionProcess;
import org.openscience.cdk.reaction.ReactionProcessTest;
import org.openscience.cdk.reaction.type.parameters.IParameterReact;
import org.openscience.cdk.reaction.type.parameters.SetReactionCenter;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.cdk.tools.manipulator.ReactionManipulator;

/**
 * TestSuite that runs a test for the RadicalSiteHrGammaReaction.
 * Generalized Reaction: [A*]-(C)_3-C4[H] => A([H])-(C_3)-[C4*].
 *
 * @cdk.module test-reaction
 */
public class RadicalSiteHrGammaReactionTest extends ReactionProcessTest {
	
	private IChemObjectBuilder builder = NoNotificationChemObjectBuilder.getInstance();
	/**
	 *  The JUnit setup method
	 */
	public  RadicalSiteHrGammaReactionTest()  throws Exception {
			setReaction(RadicalSiteHrGammaReaction.class);
	 }
	 
	 /**
	  *  The JUnit setup method
	  */
	 @Test public void testRadicalSiteHrGammaReaction() throws Exception {
			IReactionProcess type = new RadicalSiteHrGammaReaction();
			Assert.assertNotNull(type);
	 }
	/**
	 * A unit test suite for JUnit. Reaction: [A*]-C1-C2-C3[H] => A([H])-C1-C2-[C3*]
	 * Automatic search of the center active.
	 *
	 * @return    The test suite
	 */
	@Test public void testInitiate_IMoleculeSet_IMoleculeSet() throws Exception {
		IReactionProcess type = new RadicalSiteHrGammaReaction();
		
		IMolecule molecule = getMolecule();
        
		IMoleculeSet setOfReactants = DefaultChemObjectBuilder.getInstance().newMoleculeSet();
		setOfReactants.addMolecule(molecule);

		/* initiate */
		makeSureAtomTypesAreRecognized(molecule);
		
        List<IParameterReact> paramList = new ArrayList<IParameterReact>();
	    IParameterReact param = new SetReactionCenter();
        param.setParameter(Boolean.FALSE);
        paramList.add(param);
        type.setParameterList(paramList);
        IReactionSet setOfReactions = type.initiate(setOfReactants, null);
        
        Assert.assertEquals(3, setOfReactions.getReactionCount());
        Assert.assertEquals(1, setOfReactions.getReaction(0).getProductCount());

        IMolecule product = setOfReactions.getReaction(0).getProducts().getMolecule(0);
        IMolecule molecule2 = getProduct1();
        
        IQueryAtomContainer queryAtom = QueryAtomContainerCreator.createSymbolAndChargeQueryContainer(product);
        Assert.assertTrue(UniversalIsomorphismTester.isIsomorph(molecule2,queryAtom));
        
	}
	/**
	 * create the compound
	 * 
	 * @return The IMolecule
	 * @throws Exception
	 */
	private IMolecule getMolecule() throws Exception {

		IMolecule molecule = builder.newMolecule();
		molecule.addAtom(builder.newAtom("C"));
		molecule.addAtom(builder.newAtom("C"));
		molecule.addBond(0, 1, IBond.Order.SINGLE);
		molecule.addAtom(builder.newAtom("C"));
		molecule.addBond(1, 2, IBond.Order.SINGLE);
		molecule.addAtom(builder.newAtom("C"));
		molecule.addBond(2, 3, IBond.Order.SINGLE);
		molecule.addAtom(builder.newAtom("C"));
		molecule.getAtom(4).setFormalCharge(1);
		molecule.addBond(3, 4, IBond.Order.SINGLE);
		molecule.addAtom(builder.newAtom("C"));
		molecule.addBond(4, 5, IBond.Order.SINGLE);
		addExplicitHydrogens(molecule);

		molecule.getAtom(4).setFormalCharge(0);
		molecule.addSingleElectron(new SingleElectron(molecule.getAtom(4)));
		
		AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(molecule);
		makeSureAtomTypesAreRecognized(molecule);
		return molecule;
	}
	/**
	 * create the compound
	 * 
	 * @return The IMolecule
	 * @throws Exception 
	 */
	private IMolecule getProduct1() throws Exception {
		IMolecule molecule2 = builder.newMolecule();
		molecule2.addAtom(builder.newAtom("C"));
		molecule2.getAtom(0).setFormalCharge(1);
		molecule2.addAtom(builder.newAtom("C"));
		molecule2.addBond(0, 1, IBond.Order.SINGLE);
		molecule2.addAtom(builder.newAtom("C"));
		molecule2.addBond(1, 2, IBond.Order.SINGLE);
		molecule2.addAtom(builder.newAtom("C"));
		molecule2.addBond(2, 3, IBond.Order.SINGLE);
		molecule2.addAtom(builder.newAtom("C"));
		molecule2.addBond(3, 4, IBond.Order.SINGLE);
		molecule2.addAtom(builder.newAtom("C"));
		molecule2.addBond(4, 5, IBond.Order.SINGLE);
		addExplicitHydrogens(molecule2);

		molecule2.getAtom(0).setFormalCharge(0);
        molecule2.addSingleElectron(new SingleElectron(molecule2.getAtom(0)));
        
		AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(molecule2);
		makeSureAtomTypesAreRecognized(molecule2);
		return molecule2;
	}
	/**
	 * A unit test suite for JUnit.
	 * 
	 * @return    The test suite
	 */
	@Test public void testCDKConstants_REACTIVE_CENTER() throws Exception {
		IReactionProcess type  = new RadicalSiteHrGammaReaction();
		IMoleculeSet setOfReactants = builder.newMoleculeSet();

		IMolecule molecule = getMolecule();
		
		/*manually put the reactive center*/
		molecule.getAtom(4).setFlag(CDKConstants.REACTIVE_CENTER,true);
		molecule.getAtom(0).setFlag(CDKConstants.REACTIVE_CENTER,true);
		molecule.getAtom(6).setFlag(CDKConstants.REACTIVE_CENTER,true);
		molecule.getBond(5).setFlag(CDKConstants.REACTIVE_CENTER,true);
		
		setOfReactants.addMolecule(molecule);
		List<IParameterReact> paramList = new ArrayList<IParameterReact>();
	    IParameterReact param = new SetReactionCenter();
        param.setParameter(Boolean.TRUE);
        paramList.add(param);
        type.setParameterList(paramList);
        
        /* initiate */
        IReactionSet setOfReactions = type.initiate(setOfReactants, null);

        IMolecule reactant = setOfReactions.getReaction(0).getReactants().getMolecule(0);
		Assert.assertTrue(molecule.getAtom(6).getFlag(CDKConstants.REACTIVE_CENTER));
		Assert.assertTrue(reactant.getAtom(6).getFlag(CDKConstants.REACTIVE_CENTER));
		Assert.assertTrue(molecule.getAtom(0).getFlag(CDKConstants.REACTIVE_CENTER));
		Assert.assertTrue(reactant.getAtom(0).getFlag(CDKConstants.REACTIVE_CENTER));
		Assert.assertTrue(molecule.getAtom(4).getFlag(CDKConstants.REACTIVE_CENTER));
		Assert.assertTrue(reactant.getAtom(4).getFlag(CDKConstants.REACTIVE_CENTER));
		Assert.assertTrue(molecule.getBond(5).getFlag(CDKConstants.REACTIVE_CENTER));
		Assert.assertTrue(reactant.getBond(5).getFlag(CDKConstants.REACTIVE_CENTER));
	}

	/**
	 * A unit test suite for JUnit.
	 *  
	 * @return    The test suite
	 */
	@Test public void testMapping() throws Exception {
		IReactionProcess type  = new RadicalSiteHrGammaReaction();
		IMoleculeSet setOfReactants = builder.newMoleculeSet();

		IMolecule molecule = getMolecule();
		
		setOfReactants.addMolecule(molecule);
		
		/*automatic search of the center active*/
        List<IParameterReact> paramList = new ArrayList<IParameterReact>();
	    IParameterReact param = new SetReactionCenter();
        param.setParameter(Boolean.FALSE);
        paramList.add(param);
        type.setParameterList(paramList);
        
        /* initiate */
		makeSureAtomTypesAreRecognized(molecule);
		
        IReactionSet setOfReactions = type.initiate(setOfReactants, null);
        
        IMolecule product = setOfReactions.getReaction(0).getProducts().getMolecule(0);

        Assert.assertEquals(4,setOfReactions.getReaction(0).getMappingCount());
        IAtom mappedProductA1 = (IAtom)ReactionManipulator.getMappedChemObject(setOfReactions.getReaction(0), molecule.getAtom(0));
        Assert.assertEquals(mappedProductA1, product.getAtom(0));
        IAtom mappedProductA2 = (IAtom)ReactionManipulator.getMappedChemObject(setOfReactions.getReaction(0), molecule.getAtom(6));
        Assert.assertEquals(mappedProductA2, product.getAtom(6));
        IAtom mappedProductA3 = (IAtom)ReactionManipulator.getMappedChemObject(setOfReactions.getReaction(0), molecule.getAtom(4));
        Assert.assertEquals(mappedProductA3, product.getAtom(4));
        IBond mappedProductB1 = (IBond)ReactionManipulator.getMappedChemObject(setOfReactions.getReaction(0), molecule.getBond(5));
        Assert.assertEquals(mappedProductB1, product.getBond(17));        
	}
	/**
	 * Test to recognize if a IMolecule matcher correctly identifies the CDKAtomTypes.
	 * 
	 * @param molecule          The IMolecule to analyze
	 * @throws CDKException
	 */
	private void makeSureAtomTypesAreRecognized(IMolecule molecule) throws Exception {

		Iterator<IAtom> atoms = molecule.atoms().iterator();
		CDKAtomTypeMatcher matcher = CDKAtomTypeMatcher.getInstance(molecule.getBuilder());
		while (atoms.hasNext()) {
				IAtom nextAtom = atoms.next();
				Assert.assertNotNull(
					"Missing atom type for: " + nextAtom, 
					matcher.findMatchingAtomType(molecule, nextAtom)
				);
		}
	}

}

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
import org.openscience.cdk.LonePair;
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
import org.openscience.cdk.isomorphism.matchers.QueryAtomContainer;
import org.openscience.cdk.isomorphism.matchers.QueryAtomContainerCreator;
import org.openscience.cdk.nonotify.NoNotificationChemObjectBuilder;
import org.openscience.cdk.reaction.IReactionProcess;
import org.openscience.cdk.reaction.ReactionProcessTest;
import org.openscience.cdk.reaction.type.parameters.IParameterReact;
import org.openscience.cdk.reaction.type.parameters.SetReactionCenter;
import org.openscience.cdk.tools.LonePairElectronChecker;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.cdk.tools.manipulator.ReactionManipulator;
/**
 * <p>IReactionProcess which participate in movement resonance. 
 * This reaction could be represented as |A-B=C => [A+]=B-[C-]. Due to 
 * the negative charge of the atom A, the double bond in position 2 is 
 * displaced.</p>
 * <pre>
 *  IMoleculeSet setOfReactants = DefaultChemObjectBuilder.getInstance().newMoleculeSet();
 *  setOfReactants.addMolecule(new Molecule());
 *  IReactionProcess type = new RearrangementLonePairReaction();
 *  HashMap<String,Object> params = new HashMap<String,Object>();
        params.put("hasActiveCenter",Boolean.FALSE);;
    type.setParameters(params);
 *  IReactionSet setOfReactions = type.initiate(setOfReactants, null);
 *  </pre>
 * 
 * <p>We have the possibility to localize the reactive center. Good method if you
 * want to localize the reaction in a fixed point</p>
 * <pre>atoms[0].setFlag(CDKConstants.REACTIVE_CENTER,true);</pre>
 * <p>Moreover you must put the parameter Boolean.TRUE</p>
 * <p>If the reactive center is not localized then the reaction process will
 * try to find automatically the possible reactive center.</p>
 * 
 * 
 * @author         Miguel Rojas
 * 
 * @cdk.created    2006-05-05
 * @cdk.module     test-reaction
 * @cdk.set        reaction-types
 * 
 **/
/**
 * TestSuite that runs a test for the RearrangementLonePairReactionTest.
 * Generalized Reaction: [A-]-B=C => A=B-[C-].
 *
 * @cdk.module test-reaction
 */
 
public class RearrangementLonePairReactionTest extends ReactionProcessTest {
	
	private IChemObjectBuilder builder = NoNotificationChemObjectBuilder.getInstance();
	/**
	 *  The JUnit setup method
	 */
	public  RearrangementLonePairReactionTest()  throws Exception {
			setReaction(RearrangementLonePairReaction.class);
	}
	 
	 /**
	  *  The JUnit setup method
	  */
	 @Test public void testRearrangementLonePairReaction() throws Exception {
			IReactionProcess type = new RearrangementLonePairReaction();
			Assert.assertNotNull(type);
	 }
	/**
	 * A unit test suite for JUnit. Reaction: O-C=C-C => [O+]=C-[C-]-C
	 * Automatic search of the center active.
	 *
	 * @cdk.inchi  InChI=1/C3H6O/c1-2-3-4/h2-4H,1H3
	 * 
	 * @return    The test suite
	 */
	@Test public void testInitiate_IMoleculeSet_IMoleculeSet() throws Exception {
		IReactionProcess type = new RearrangementLonePairReaction();
		
		IMolecule molecule = getMolecule1();
		
		IMoleculeSet setOfReactants = DefaultChemObjectBuilder.getInstance().newMoleculeSet();
		setOfReactants.addMolecule(molecule);

		/* initiate */
		
        List<IParameterReact> paramList = new ArrayList<IParameterReact>();
	    IParameterReact param = new SetReactionCenter();
        param.setParameter(Boolean.FALSE);
        paramList.add(param);
        type.setParameterList(paramList);
        IReactionSet setOfReactions = type.initiate(setOfReactants, null);
        
        Assert.assertEquals(1, setOfReactions.getReactionCount());
        Assert.assertEquals(1, setOfReactions.getReaction(0).getProductCount());

        IMolecule product = setOfReactions.getReaction(0).getProducts().getMolecule(0);
        Assert.assertEquals(-1, product.getAtom(2).getFormalCharge().intValue());
        Assert.assertEquals(0, product.getConnectedLonePairsCount(molecule.getAtom(1)));
        
        IMolecule molecule2 = getMolecule2();
        
        IQueryAtomContainer queryAtom = QueryAtomContainerCreator.createSymbolAndChargeQueryContainer(product);
        Assert.assertTrue(UniversalIsomorphismTester.isIsomorph(molecule2,queryAtom));
        
	}

	/**
	 * A unit test suite for JUnit. Reaction: O-C=C-C => [O+]=C-[C-]-C
	 * Manually put of the center active.
	 *
	 * @cdk.inchi  InChI=1/C3H6O/c1-2-3-4/h2-4H,1H3
	 * 
	 * @return    The test suite
	 */
	@Test public void testManuallyCentreActive() throws Exception {
		IReactionProcess type = new RearrangementLonePairReaction();
		
		IMoleculeSet setOfReactants = DefaultChemObjectBuilder.getInstance().newMoleculeSet();
		/*[C-]-C=C-C */
		IMolecule molecule = getMolecule1();
		setOfReactants.addMolecule(molecule);
		
		/*manually put the center active*/
		molecule.getAtom(0).setFlag(CDKConstants.REACTIVE_CENTER,true);
		molecule.getBond(0).setFlag(CDKConstants.REACTIVE_CENTER,true);
		molecule.getAtom(1).setFlag(CDKConstants.REACTIVE_CENTER,true);
		molecule.getBond(1).setFlag(CDKConstants.REACTIVE_CENTER,true);
		molecule.getAtom(2).setFlag(CDKConstants.REACTIVE_CENTER,true);
		molecule.getBond(2).setFlag(CDKConstants.REACTIVE_CENTER,true);
		molecule.getAtom(3).setFlag(CDKConstants.REACTIVE_CENTER,true);
		
        List<IParameterReact> paramList = new ArrayList<IParameterReact>();
	    IParameterReact param = new SetReactionCenter();
        param.setParameter(Boolean.TRUE);
        paramList.add(param);
        type.setParameterList(paramList);
        
		/* initiate */
		
        IReactionSet setOfReactions = type.initiate(setOfReactants, null);
        
        Assert.assertEquals(1, setOfReactions.getReactionCount());
        Assert.assertEquals(1, setOfReactions.getReaction(0).getProductCount());

        IMolecule product = setOfReactions.getReaction(0).getProducts().getMolecule(0);
        
        /*C=C-[C-]-C*/
        IMolecule molecule2 = getMolecule2();
        
        IQueryAtomContainer queryAtom = QueryAtomContainerCreator.createSymbolAndChargeQueryContainer(product);
        Assert.assertTrue(UniversalIsomorphismTester.isIsomorph(molecule2,queryAtom));
        
	}

	/**
	 * A unit test suite for JUnit.
	 * 
	 * @return    The test suite
	 */
	@Test public void testCDKConstants_REACTIVE_CENTER() throws Exception {
		IReactionProcess type  = new RearrangementLonePairReaction();
		IMoleculeSet setOfReactants = builder.newMoleculeSet();

		IMolecule molecule = getMolecule1();
		
		/*manually put the reactive center*/
		molecule.getAtom(0).setFlag(CDKConstants.REACTIVE_CENTER,true);
		molecule.getAtom(1).setFlag(CDKConstants.REACTIVE_CENTER,true);
		molecule.getAtom(2).setFlag(CDKConstants.REACTIVE_CENTER,true);
		molecule.getBond(0).setFlag(CDKConstants.REACTIVE_CENTER,true);
		molecule.getBond(1).setFlag(CDKConstants.REACTIVE_CENTER,true);
		
		setOfReactants.addMolecule(molecule);
		List<IParameterReact> paramList = new ArrayList<IParameterReact>();
	    IParameterReact param = new SetReactionCenter();
        param.setParameter(Boolean.TRUE);
        paramList.add(param);
        type.setParameterList(paramList);
        
        /* initiate */
        IReactionSet setOfReactions = type.initiate(setOfReactants, null);

        IMolecule reactant = setOfReactions.getReaction(0).getReactants().getMolecule(0);
		Assert.assertTrue(molecule.getAtom(0).getFlag(CDKConstants.REACTIVE_CENTER));
		Assert.assertTrue(reactant.getAtom(0).getFlag(CDKConstants.REACTIVE_CENTER));
		Assert.assertTrue(molecule.getAtom(1).getFlag(CDKConstants.REACTIVE_CENTER));
		Assert.assertTrue(reactant.getAtom(1).getFlag(CDKConstants.REACTIVE_CENTER));
		Assert.assertTrue(molecule.getAtom(2).getFlag(CDKConstants.REACTIVE_CENTER));
		Assert.assertTrue(reactant.getAtom(2).getFlag(CDKConstants.REACTIVE_CENTER));
		Assert.assertTrue(molecule.getBond(0).getFlag(CDKConstants.REACTIVE_CENTER));
		Assert.assertTrue(reactant.getBond(0).getFlag(CDKConstants.REACTIVE_CENTER));
		Assert.assertTrue(molecule.getBond(1).getFlag(CDKConstants.REACTIVE_CENTER));
		Assert.assertTrue(reactant.getBond(1).getFlag(CDKConstants.REACTIVE_CENTER));
	}

	/**
	 * A unit test suite for JUnit.
	 *  
	 * @return    The test suite
	 */
	@Test public void testMapping() throws Exception {
		IReactionProcess type = new RearrangementLonePairReaction();
		
		IMoleculeSet setOfReactants = DefaultChemObjectBuilder.getInstance().newMoleculeSet();
		/*[C-]-C=C-C*/
		IMolecule molecule = getMolecule1();
		molecule.addLonePair(new LonePair(molecule.getAtom(0)));
		setOfReactants.addMolecule(molecule);
		
		/*automatic search of the center active*/
        List<IParameterReact> paramList = new ArrayList<IParameterReact>();
	    IParameterReact param = new SetReactionCenter();
        param.setParameter(Boolean.FALSE);
        paramList.add(param);
        type.setParameterList(paramList);
		/* initiate */
		
        IReactionSet setOfReactions = type.initiate(setOfReactants, null);
        
        IMolecule product = setOfReactions.getReaction(0).getProducts().getMolecule(0);

        Assert.assertEquals(5,setOfReactions.getReaction(0).getMappingCount());
        IAtom mappedProductA1 = (IAtom)ReactionManipulator.getMappedChemObject(setOfReactions.getReaction(0), molecule.getAtom(0));
        Assert.assertEquals(mappedProductA1, product.getAtom(0));
        mappedProductA1 = (IAtom)ReactionManipulator.getMappedChemObject(setOfReactions.getReaction(0), molecule.getAtom(1));
        Assert.assertEquals(mappedProductA1, product.getAtom(1));
        mappedProductA1 = (IAtom)ReactionManipulator.getMappedChemObject(setOfReactions.getReaction(0), molecule.getAtom(2));
        Assert.assertEquals(mappedProductA1, product.getAtom(2));
        
	}
	/**
	 * Test to recognize if this IMolecule_1 matches correctly into the CDKAtomTypes.
	 */
	@Test public void testAtomTypesMolecule1() throws Exception{
		IMolecule moleculeTest = getMolecule1();
		makeSureAtomTypesAreRecognized(moleculeTest);
		
	}

	/**
	 * Test to recognize if this IMolecule_2 matches correctly into the CDKAtomTypes.
	 */
	@Test public void testAtomTypesMolecule2() throws Exception{
		IMolecule moleculeTest = getMolecule2();
		makeSureAtomTypesAreRecognized(moleculeTest);
		
	}
	/**
	 * get the molecule 1: O-C=C-C
	 * 
	 * @cdk.inchi  InChI=1/C3H6O/c1-2-3-4/h2-4H,1H3
	 * 
	 * @return The IMolecule
	 */
	private IMolecule getMolecule1()throws Exception {
		IMolecule molecule = builder.newMolecule();
		molecule.addAtom(builder.newAtom("O"));
		molecule.addLonePair(new LonePair(molecule.getAtom(0)));
		molecule.addAtom(builder.newAtom("C"));
		molecule.addBond(0, 1, IBond.Order.SINGLE);
		molecule.addAtom(builder.newAtom("C"));
		molecule.addBond(1, 2, IBond.Order.DOUBLE);
		molecule.addAtom(builder.newAtom("C"));
		molecule.addBond(2, 3, IBond.Order.SINGLE);
		
		addExplicitHydrogens(molecule);
		AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(molecule);

        LonePairElectronChecker lpcheck = new LonePairElectronChecker();
        lpcheck.saturate(molecule);
        
		return molecule;
	}
	/**
	 * get the molecule 2: [O+]=C-[C-]-C
	 * 
	 * @return The IMolecule
	 */
	private IMolecule getMolecule2()throws Exception {
		IMolecule molecule = builder.newMolecule();
		molecule.addAtom(builder.newAtom("O"));
		molecule.getAtom(0).setFormalCharge(+1);
		molecule.addAtom(builder.newAtom("C"));
		molecule.addBond(0, 1, IBond.Order.DOUBLE);
		molecule.addAtom(builder.newAtom("C"));
		molecule.getAtom(2).setFormalCharge(-1);
		molecule.addLonePair(new LonePair(molecule.getAtom(2)));
		molecule.addBond(1, 2, IBond.Order.SINGLE);
		molecule.addAtom(builder.newAtom("C"));
		molecule.addBond(2, 3, IBond.Order.SINGLE);
		
		addExplicitHydrogens(molecule);
		AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(molecule);

        LonePairElectronChecker lpcheck = new LonePairElectronChecker();
        lpcheck.saturate(molecule);
        
        return molecule;
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
	/**
	 * A unit test suite for JUnit: Resonance Fluorobenzene  Fc1ccccc1 <=> ...
	 *
	 * InChI=1/C6H5F/c7-6-4-2-1-3-5-6/h1-5H
	 *
	 * @return    The test suite
	 */
	@Test public void testFluorobenzene() throws Exception {

		 IMolecule molecule = builder.newMolecule();
		 molecule.addAtom(builder.newAtom("F"));
		 molecule.addAtom(builder.newAtom("C"));
		 molecule.addBond(0, 1, IBond.Order.SINGLE);
		 molecule.addAtom(builder.newAtom("C"));
		 molecule.addBond(1, 2, IBond.Order.DOUBLE);
		 molecule.addAtom(builder.newAtom("C"));
		 molecule.addBond(2, 3, IBond.Order.SINGLE);
		 molecule.addAtom(builder.newAtom("C"));
		 molecule.addBond(3, 4, IBond.Order.DOUBLE);
		 molecule.addAtom(builder.newAtom("C"));
		 molecule.addBond(4, 5, IBond.Order.SINGLE);
		 molecule.addAtom(builder.newAtom("C"));
		 molecule.addBond(5, 6, IBond.Order.DOUBLE);
		 molecule.addBond(6, 1, IBond.Order.SINGLE);
		
		addExplicitHydrogens(molecule);
		AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(molecule);
        LonePairElectronChecker lpcheck = new LonePairElectronChecker();
        lpcheck.saturate(molecule);
		
		IReactionProcess type = new RearrangementLonePairReaction();
		
		IMoleculeSet setOfReactants = DefaultChemObjectBuilder.getInstance().newMoleculeSet();
		setOfReactants.addMolecule(molecule);
		/*automatic search of the center active*/
        List<IParameterReact> paramList = new ArrayList<IParameterReact>();
	    IParameterReact param = new SetReactionCenter();
        param.setParameter(Boolean.FALSE);
        paramList.add(param);
        type.setParameterList(paramList);
		/* initiate */
        
        IReactionSet setOfReactions = type.initiate(setOfReactants, null);
        Assert.assertEquals(1, setOfReactions.getReactionCount());
        Assert.assertEquals(1, setOfReactions.getReaction(0).getProductCount());
        IMolecule product1 = setOfReactions.getReaction(0).getProducts().getMolecule(0);
        
        IMolecule molecule1 = builder.newMolecule();
	 	molecule1.addAtom(builder.newAtom("F"));
	 	molecule1.getAtom(0).setFormalCharge(1);
		molecule1.addAtom(builder.newAtom("C"));
		molecule1.addBond(0, 1, IBond.Order.DOUBLE);
		molecule1.addAtom(builder.newAtom("C"));
	 	molecule1.getAtom(2).setFormalCharge(-1);
		molecule1.addBond(1, 2, IBond.Order.SINGLE);
		molecule1.addAtom(builder.newAtom("C"));
		molecule1.addBond(2, 3, IBond.Order.SINGLE);
		molecule1.addAtom(builder.newAtom("C"));
		molecule1.addBond(3, 4, IBond.Order.DOUBLE);
		molecule1.addAtom(builder.newAtom("C"));
		molecule1.addBond(4, 5, IBond.Order.SINGLE);
		molecule1.addAtom(builder.newAtom("C"));
		molecule1.addBond(5, 6, IBond.Order.DOUBLE);
		molecule1.addBond(6, 1, IBond.Order.SINGLE);
		addExplicitHydrogens(molecule1);
		AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(molecule1);
		lpcheck.saturate(molecule1);
		
		QueryAtomContainer qAC = QueryAtomContainerCreator.createSymbolAndChargeQueryContainer(molecule1);
	    Assert.assertTrue(UniversalIsomorphismTester.isIsomorph(product1,qAC));
        
	}
}

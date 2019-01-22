/* $Revision$ $Author$ $Date$
 *
 * Copyright (C) 2007  Egon Willighagen <egonw@users.sf.net>
 *
 * Contact: cdk-devel@lists.sourceforge.net
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
package org.openscience.cdk.aromaticity;

import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.annotations.TestClass;
import org.openscience.cdk.annotations.TestMethod;
import org.openscience.cdk.atomtype.CDKAtomTypeMatcher;
import org.openscience.cdk.config.AtomTypeFactory;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.graph.ConnectivityChecker;
import org.openscience.cdk.graph.SpanningTree;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomType;
import org.openscience.cdk.interfaces.IAtomType.Hybridization;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IRingSet;
import org.openscience.cdk.ringsearch.AllRingsFinder;
import org.openscience.cdk.ringsearch.SSSRFinder;

import java.util.List;

/**
 * This aromaticity detector detects the aromaticity based on the H&uuml;ckel
 * 4n+2 pi-electrons rule applied to isolated ring systems. It assumes
 * CDK atom types to be perceived with the {@link CDKAtomTypeMatcher} or with
 * any compatible class. For example:
 * <pre>
 * Molecule molecule = MoleculeFactory.makePyridineOxide();
 * AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(molecule);
 * CDKHueckelAromaticityDetector.detectAromaticity(molecule);
 * </pre>
 *
 * @author         egonw
 * @cdk.module     standard
 * @cdk.githash
 * @cdk.created    2007-10-05
 * 
 * @see org.openscience.cdk.CDKConstants
 */
@TestClass("org.openscience.cdk.aromaticity.CDKHueckelAromaticityDetectorTest")
public class CDKHueckelAromaticityDetector {

	private static AtomTypeFactory factory = null;

    @TestMethod("testDetectAromaticity_IAtomContainer")
    public static boolean detectAromaticity(IAtomContainer atomContainer) throws CDKException {
		SpanningTree spanningTree = new SpanningTree(atomContainer);
		IAtomContainer ringSystems = spanningTree.getCyclicFragmentsContainer();
		if (ringSystems.getAtomCount() == 0) {
			// If there are no rings, then there cannot be any aromaticity
			return false;
		}
		// FIXME: should not really mark them here
        for (IAtom atom : ringSystems.atoms()) atom.setFlag(CDKConstants.ISINRING, true);
        for (IBond bond : ringSystems.bonds()) bond.setFlag(CDKConstants.ISINRING, true);

		boolean foundSomeAromaticity = false;
        for (IAtomContainer isolatedSystem : ConnectivityChecker.partitionIntoMolecules(ringSystems).atomContainers()) {
			IRingSet singleRings = new SSSRFinder(isolatedSystem).findSSSR();
			int maxRingSize = 20;
			boolean atLeastOneRingIsSprouted = false;
			// test single rings in SSSR
            for (IAtomContainer singleRing : singleRings.atomContainers()) {
				if (singleRing.getAtomCount() > maxRingSize) maxRingSize = singleRing.getAtomCount();
				if (isRingSystemSproutedWithNonRingDoubleBonds(atomContainer, singleRing)) {
					// OK, this ring is not aromatic
					atLeastOneRingIsSprouted = true;
				} else {
					// possibly aromatic
					foundSomeAromaticity |= isRingAllSP2AndHueckelValid(singleRing);
				}
			}
			// OK, what about the one larger ring (if no aromaticity found in SSSR)?
			if (!foundSomeAromaticity && !atLeastOneRingIsSprouted &&
				singleRings.getAtomContainerCount() <= 3) {
				// every ring system consisting of more than two rings is too difficult
                IRingSet allRings = new AllRingsFinder().findAllRingsInIsolatedRingSystem(isolatedSystem);
                for (IAtomContainer ring : allRings.atomContainers()) {
					// there should be exactly three rings, of which only one has a size larger
					// than the two previous ones
					if (ring.getAtomCount() <= maxRingSize) {
						// possibly aromatic
						foundSomeAromaticity |= isRingAllSP2AndHueckelValid(ring);
					}
				}
			}
		}
		
		return foundSomeAromaticity;
	}
	
	/**
	 * Tests if the method contains only SP2 atoms, and if the electron count
	 * matches the H&uuml;ckel 4n+2 rule.
	 */
	private static boolean isRingAllSP2AndHueckelValid(IAtomContainer singleRing) throws CDKException {
		boolean foundSomeAromaticity = false;
		// OK, all atoms must be sp2, or sp3 with a lone pair
		boolean allAtomsSuitable = true;
		int electronCount = 0;
		for (IAtom ringAtom : singleRing.atoms()) {
			if (ringAtom.getHybridization() != CDKConstants.UNSET &&
				(ringAtom.getHybridization() == Hybridization.SP2) ||
                    ringAtom.getHybridization() == Hybridization.PLANAR3) {
				// for example, a carbon
				// note: the double bond is in the ring, that has been tested earlier
				// FIXME: this does assume bond orders to be resolved too, when detecting
				// sprouting double bonds
				if ("N.planar3".equals(ringAtom.getAtomTypeName())) {
					electronCount += 2;
				} else if ("N.minus.planar3".equals(ringAtom.getAtomTypeName())) {
					electronCount += 2;
				} else if ("S.2".equals(ringAtom.getAtomTypeName())) {
                    electronCount += 2;
				} else if ("S.planar3".equals(ringAtom.getAtomTypeName())) {
                    electronCount += 2;
                } else if ("C.minus.planar".equals(ringAtom.getAtomTypeName())) {
                    electronCount += 2;
                } else if ("O.planar3".equals(ringAtom.getAtomTypeName())) {
                    electronCount += 2;
                } else {
                    if (factory == null) {
						factory = AtomTypeFactory.getInstance(
							"org/openscience/cdk/dict/data/cdk-atom-types.owl",
							ringAtom.getBuilder()
						);
					}
					IAtomType type = factory.getAtomType(ringAtom.getAtomTypeName());
					Object property = type.getProperty(CDKConstants.PI_BOND_COUNT);
					if (property != null && property instanceof Integer) {
						electronCount += (Integer) property;
					}
				}
			} else if (ringAtom.getHybridization() != null &&
					   ringAtom.getHybridization() == Hybridization.SP3 &&
					   getLonePairCount(ringAtom) > 0) {
				// for example, a nitrogen or oxygen
				electronCount += 2;
			} else {
				// OK, not aromatic, need to break the while loop now
				allAtomsSuitable = false;
			}
		}
		boolean fourNplusTwoRule = (electronCount % 4 == 2) && (electronCount > 2);
		if (allAtomsSuitable && fourNplusTwoRule) {
			foundSomeAromaticity = true;
			markRingAtomsAndBondsAromatic(singleRing);
		}
		return foundSomeAromaticity;
	}

	/**
	 * Determines if the isolatedRingSystem has attached double bonds, which are not part of the ring system itself,
	 * and not part of any other ring system.
	 */
	private static boolean isRingSystemSproutedWithNonRingDoubleBonds(IAtomContainer fullContainer, IAtomContainer isolatedRingSystem) {
		for (IAtom atom : isolatedRingSystem.atoms()) {
			List<IBond> neighborBonds = fullContainer.getConnectedBondsList(atom);
			for (IBond neighborBond : neighborBonds) {
				if (!neighborBond.getFlag(CDKConstants.ISINRING) &&
					neighborBond.getOrder() == CDKConstants.BONDORDER_DOUBLE ||
					neighborBond.getOrder() == CDKConstants.BONDORDER_TRIPLE) {
					return true;
				}
			}
		}
		return false;
	}
	
	private static int getLonePairCount(IAtom atom) {
		Integer count = (Integer)atom.getProperty(CDKConstants.LONE_PAIR_COUNT);
		if (count == null) {
			return 0;
		} else {
			return count;
		}
	}
	
	private static void markRingAtomsAndBondsAromatic(IAtomContainer container) {
	    for (IAtom atom : container.atoms()) atom.setFlag(CDKConstants.ISAROMATIC, true);
		  for (IBond bond : container.bonds()) bond.setFlag(CDKConstants.ISAROMATIC, true);
	}
}


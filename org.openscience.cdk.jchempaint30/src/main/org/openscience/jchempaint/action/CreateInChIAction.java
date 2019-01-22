/*
 *  $RCSfile$
 *  $Author: egonw $
 *  $Date: 2007-01-04 17:26:00 +0000 (Thu, 04 Jan 2007) $
 *  $Revision: 7634 $
 *
 *  Copyright (C) 1997-2008 Christoph Steinbeck, Stefan Kuhn
 *
 *  Contact: cdk-jchempaint@lists.sourceforge.net
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public License
 *  as published by the Free Software Foundation; either version 2.1
 *  of the License, or (at your option) any later version.
 *  All we ask is that proper credit is given for our work, which includes
 *  - but is not limited to - adding the above copyright notice to the beginning
 *  of your source code files, and to any copyright notice that you may distribute
 *  with programs based on this work.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 */
package org.openscience.jchempaint.action;

import java.awt.event.ActionEvent;

import javax.swing.JFrame;

import net.sf.jniinchi.INCHI_RET;

import org.openscience.cdk.ChemModel;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.inchi.InChIGenerator;
import org.openscience.cdk.inchi.InChIGeneratorFactory;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.jchempaint.GT;
import org.openscience.jchempaint.dialog.TextViewDialog;

/**
 * Creates an InChI from the current model
 *
 * @cdk.module jchempaint
 * @author     Sam Adams
 */
public class CreateInChIAction extends JCPAction
{

	private static final long serialVersionUID = -4886982931009753347L;
	
	TextViewDialog dialog = null;
	JFrame frame = null;

	public void actionPerformed(ActionEvent e)
	{
		logger.debug("Trying to create InChI: ", type);
		
		if (dialog == null)
		{
			dialog = new TextViewDialog(frame, "InChI", null, false, 40, 2);
		}
        
        InChIGeneratorFactory factory = null;
        try {
            factory = InChIGeneratorFactory.getInstance();
        } catch (CDKException cdke) {
            String message = GT._("Error loading InChI library:")+" " + cdke.getMessage();
            logger.error(message);
            logger.debug(cdke);
            dialog.setMessage(GT._("Error"), message);
        }
        
        if (factory != null) {
            ChemModel model = (ChemModel) jcpPanel.getChemModel();
            if(model.getReactionSet()!=null && model.getReactionSet().getReactionCount()>0){
            	dialog.setMessage(GT._("Problem"), GT._("You have reactions in JCP. Reactions cannot be shown as InChI!"));
            }else{
            	StringBuffer dialogText=new StringBuffer();
            	int i=0;
            	String eol = System.getProperty("line.separator");
				for(IAtomContainer container : model.getMoleculeSet().molecules())
				{
	                if (model.getMoleculeSet().getAtomContainerCount() > 1) {
	                    dialogText.append(GT._("Structure")+" #" + (i+1) + eol);
	                }
	                try {
	                    InChIGenerator inchiGen = factory.getInChIGenerator(container);
	                    INCHI_RET ret = inchiGen.getReturnStatus();
	                    String inchi = inchiGen.getInchi();
	                    String auxinfo = inchiGen.getAuxInfo();
	                    String message = inchiGen.getMessage();
	                    if (ret == INCHI_RET.OKAY) {
	                        dialogText.append(inchi + eol + auxinfo + eol + eol);
	                    } else if (ret == INCHI_RET.WARNING) {
	                        dialogText.append(inchi + eol + auxinfo + eol + "Warning: " + message + eol + eol);
	                    } else {
	                        dialogText.append(GT._("InChI generation failed")+" (" + ret.toString() + ")" + eol + message + eol + eol);
	                    }
	                } catch (CDKException cdke) {
	                    dialogText.append(GT._("InChI generation failed")+": " + cdke.getMessage() + eol + eol);
	                }
	                i++;
	            }
	            dialog.setMessage(GT._("Generated InChI")+":", dialogText.toString());
            }
	            
        }
        
        
		dialog.setVisible(true);
	}
}


/*
 *  $RCSfile$
 *  $Author: egonw $
 *  $Date: 2007-01-04 17:26:00 +0000 (Thu, 04 Jan 2007) $
 *  $Revision: 7634 $
 *
 *  Copyright (C) 1997-2008 Stefan Kuhn, Christoph Steinbeck
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

import org.openscience.jchempaint.GT;
import org.openscience.jchempaint.dialog.HelpDialog;

/**
 * Pops up the help.
 *
 */
public class HelpAction extends JCPAction
{

	private static final long serialVersionUID = -9213900779679488824L;

	public void actionPerformed(ActionEvent e)
	{
		if (type.equals("tutorial"))
		{
			new HelpDialog(null, "org/openscience/jchempaint/resources/userhelp_jcp/contain/tutorial.html", GT._("JChemPaint Help")).setVisible(true);
		} else if (type.equals("feedback"))
        {
            new HelpDialog(null, "org/openscience/jchempaint/resources/userhelp_jcp/contain/feedback.html", GT._("JChemPaint Help")).setVisible(true);
        } else if (type.equals("license"))
        {
            new HelpDialog(null, "org/openscience/jchempaint/resources/userhelp_jcp/license.html", GT._("JChemPaint License")).setVisible(true);
        } else
        {
            new HelpDialog(null, "org/openscience/jchempaint/resources/userhelp_jcp/jcp.html", GT._("JChemPaint Help")).setVisible(true);
    	}
	}
}


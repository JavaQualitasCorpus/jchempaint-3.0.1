/* Copyright (C) 2009  Egon Willighagen <egonw@users.lists.sf>
 *
 * Contact: cdk-devel@lists.sourceforge.net
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 * All I ask is that proper credit is given for my work, which includes
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
package org.openscience.jchempaint.renderer;

import javax.vecmath.Point2d;

/**
 * Interface that all 2D renderers implement.
 *
 * @cdk.module render
 */
public interface IRenderer {

    /**
     * Returns the drawing model, giving access to drawing parameters.
     * 
     * @return the rendering model
     */
	public RendererModel getRenderer2DModel();

	/**
	 * Converts screen coordinates into model (or world) coordinates.
	 *
	 * @param screenXTo the screen's x coordinate
	 * @param screenYTo the screen's y coordinate
	 * @return          the matching model coordinates
	 *
	 * @see #toScreenCoordinates(double, double)
	 */
	public Point2d toModelCoordinates(double screenXTo, double screenYTo);

    /**
     * Converts model (or world) coordinates into screen coordinates.
     *
     * @param screenXTo the model's x coordinate
     * @param screenYTo the model's y coordinate
     * @return          the matching screen coordinates
     *
     * @see #toModelCoordinates(double, double)
     */
	public Point2d toScreenCoordinates(double screenXTo, double screenYTo);
	
	/**
	 * Set a new zoom factor.
	 *
	 * @param zoomFactor the new zoom factor
	 */
	public void setZoom(double zoomFactor);

    /**
     * Set a new drawing center in screen coordinates.
     *
     * @param zoomFactor the new new drawing center
     */
	public void shiftDrawCenter(double screenX, double screenY);

	/**
	 * Sets the mouse cursor shown on the renderPanel.
	 * 
	 * @param cursor One of the constants from java.awt.Cursor.
	 */
	public void setCursor(int cursor);
	/**
	 * Tells the mouse cursor shown on the renderPanel.
	 * 
	 * @return One of the constants from java.awt.Cursor.
	 */
	public int getCursor();
}

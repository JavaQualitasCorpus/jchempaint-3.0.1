/*
 *  $RCSfile$
 *  $Author: egonw $
 *  $Date: 2007-01-04 17:26:00 +0000 (Thu, 04 Jan 2007) $
 *  $Revision: 7634 $
 *
 *  Copyright (C) 1997-2008 Stefan Kuhn
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
package org.openscience.jchempaint;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Container;
import java.awt.Cursor;
import java.awt.Image;
import java.awt.event.ActionEvent;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.awt.event.WindowListener;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.EventObject;
import java.util.List;

import javax.swing.JButton;
import javax.swing.JComponent;
import javax.swing.JFrame;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JToolBar;
import javax.swing.SwingConstants;
import javax.swing.filechooser.FileFilter;
import javax.swing.undo.UndoManager;

import org.openscience.cdk.Atom;
import org.openscience.cdk.Bond;
import org.openscience.cdk.ChemModel;
import org.openscience.cdk.PseudoAtom;
import org.openscience.cdk.config.IsotopeFactory;
import org.openscience.cdk.event.ICDKChangeListener;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IChemModel;
import org.openscience.cdk.interfaces.IIsotope;
import org.openscience.cdk.tools.manipulator.ChemModelManipulator;
import org.openscience.jchempaint.action.SaveAction;
import org.openscience.jchempaint.applet.JChemPaintAbstractApplet;
import org.openscience.jchempaint.applet.JChemPaintEditorApplet;
import org.openscience.jchempaint.controller.AddAtomModule;
import org.openscience.jchempaint.controller.ControllerHub;
import org.openscience.jchempaint.controller.IChangeModeListener;
import org.openscience.jchempaint.controller.IChemModelEventRelayHandler;
import org.openscience.jchempaint.controller.IControllerModule;
import org.openscience.jchempaint.controller.MoveModule;
import org.openscience.jchempaint.renderer.RendererModel;
import org.openscience.jchempaint.renderer.selection.AbstractSelection;

public class JChemPaintPanel extends AbstractJChemPaintPanel implements
        IChemModelEventRelayHandler, ICDKChangeListener, KeyListener, IChangeModeListener {

    private static final long serialVersionUID = -8932765332441119177L;
    private JComponent lastActionButton;
    private JComponent lastSecondaryButton;
    private File currentWorkDirectory;
    private File lastOpenedFile;
    private FileFilter currentOpenFileFilter;
    private File isAlreadyAFile;
    private boolean isModified = false;
    private FileFilter currentSaveFileFilter;
    public static List<JChemPaintPanel> instances = new ArrayList<JChemPaintPanel>();
    private boolean showInsertTextField = true;
    private JPanel topContainer = null;
    private JPanel centerContainer = null;
    private boolean showToolBar = true;
    private boolean showMenuBar = true;
    private JMenuBar menu;
    private JToolBar uppertoolbar;
    private JToolBar lefttoolbar;
    private JToolBar lowertoolbar;
    private JToolBar righttoolbar;
    protected JMenuItem undoMenu;
    protected JMenuItem redoMenu;
    protected JMenu atomMenu;
    protected JMenu bondMenu;
    private boolean debug=false;
    private String lastSelectId;

	/**
     * Builds a JCPPanel with a certain model and a certain gui
     *
     * @param chemModel   The model
     * @param gui         The gui configuration string
     * @param debug       Should we be in debug mode?
     * @param applet      If this panel is to be in an applet, pass the applet here, else null.
     */
    public JChemPaintPanel(IChemModel chemModel, String gui, boolean debug, JChemPaintAbstractApplet applet) {
        GT.setLanguage(JCPPropertyHandler.getInstance().getJCPProperties().getProperty("General.language"));
        this.guistring = gui;
        this.debug = debug;
        this.setLayout(new BorderLayout());
        topContainer = new JPanel(new BorderLayout());
        topContainer.setLayout(new BorderLayout());
        this.add(topContainer, BorderLayout.NORTH);
        try {
			renderPanel = new RenderPanel(chemModel, getWidth(), getHeight(), false, debug, false, applet);
		} catch (IOException e) {
			announceError(e);
		}
        renderPanel.getHub().addChangeModeListener(this);
        renderPanel.setName("renderpanel");
        centerContainer=new JPanel();
        centerContainer.setLayout(new BorderLayout());
        centerContainer.add(new JScrollPane(renderPanel), BorderLayout.CENTER);
        this.add(centerContainer);

        customizeView();
        updateUndoRedoControls();
        SwingPopupModule inputAdapter = new SwingPopupModule(renderPanel,
                renderPanel.getHub());
        setupPopupMenus(inputAdapter);
        renderPanel.getHub().registerGeneralControllerModule(inputAdapter);
        renderPanel.getHub().setEventHandler(this);
        renderPanel.getRenderer().getRenderer2DModel().addCDKChangeListener(
                this);
        instances.add(this);
        //we set this to true always, the user should have no option to switch it off
        renderPanel.getHub().getController2DModel().setAutoUpdateImplicitHydrogens(true);
        this.addKeyListener(this);
        renderPanel.addMouseListener(new MouseAdapter(){
            public void mouseExited(MouseEvent e) {
                //this avoids ghost phantom rings if the user leaves the panel
                JChemPaintPanel.this.get2DHub().clearPhantoms();
                JChemPaintPanel.this.get2DHub().updateView();
            }            
        });
    }

    /**
     * Gets the top level container (JFrame, Applet) of this panel.
     * 
     * @return The top level container.
     */
    public Container getTopLevelContainer() {
        return this.getParent().getParent().getParent().getParent();
    }

    /**
     * If this panel is in a JFrame, sets the title of the JFrame.
     * 
     * @param title The title to set.
     */
    public void setTitle(String title) {
        Container topLevelContainer = this.getTopLevelContainer();
        if (topLevelContainer instanceof JFrame) {
            ((JFrame) topLevelContainer).setTitle(title);
        }
    }

    /**
     * Installs popup menus for this panel.
     * 
     * @param inputAdapter The SwingPopupModule to use for the popup menus.
     */
    public void setupPopupMenus(SwingPopupModule inputAdapter) {
        if (inputAdapter.getPopupMenu(PseudoAtom.class) == null) {
            inputAdapter.setPopupMenu(PseudoAtom.class,
                    new JChemPaintPopupMenu(this, "pseudo", this.guistring));
        }
        if (inputAdapter.getPopupMenu(Atom.class) == null) {
            inputAdapter.setPopupMenu(Atom.class, new JChemPaintPopupMenu(this,
                    "atom", this.guistring));
        }
        if (inputAdapter.getPopupMenu(Bond.class) == null) {
            inputAdapter.setPopupMenu(Bond.class, new JChemPaintPopupMenu(this,
                    "bond", this.guistring));
        }
        if (inputAdapter.getPopupMenu(ChemModel.class) == null) {
            inputAdapter.setPopupMenu(ChemModel.class, new JChemPaintPopupMenu(
                    this, "chemmodel", this.guistring));
        }
        /*if (inputAdapter.getPopupMenu(Reaction.class) == null) {
            inputAdapter.setPopupMenu(Reaction.class, new JChemPaintPopupMenu(
                    this, "reaction", this.guistring));
        }*/
    }

    /**
     * Called to force a re-centring of the displayed structure.
     *
     * @param isNewChemModel
     */
    public void setIsNewChemModel(boolean isNewChemModel) {
        this.renderPanel.setIsNewChemModel(isNewChemModel);
    }

    /**
     * Helps in keeping the current action button highlighted
     *
     * @return The last action button used
     */
    public JComponent getLastActionButton() {
        return lastActionButton;
    }

    /**
     * Allows setting of the is modified stage (e. g. after save)
     *
     * @param isModified
     *            is modified
     */
    public void setModified(boolean isModified) {
        this.isModified = isModified;
        Container c = this.getTopLevelContainer();
        if (c instanceof JFrame) {
            String id = renderPanel.getChemModel().getID();
            if (isModified)
                ((JFrame) c).setTitle(id + "*");
            else
                ((JFrame) c).setTitle(id);
        }
    }

    /**
     * Helps in keeping the current action button highlighted - needs to be set
     * if a new action button is choosen
     *
     * @param actionButton
     *            The new action button
     */
    public void setLastActionButton(JComponent actionButton) {
        lastActionButton = actionButton;
    }

    /**
     * Gets the currentWorkDirectory attribute of the JChemPaintPanel object
     *
     *@return The currentWorkDirectory value
     */
    public File getCurrentWorkDirectory() {
        return currentWorkDirectory;
    }

    /**
     * Sets the currentWorkDirectory attribute of the JChemPaintPanel object
     *
     *@param cwd
     *            The new currentWorkDirectory value
     */
    public void setCurrentWorkDirectory(File cwd) {
        this.currentWorkDirectory = cwd;
    }

    /**
     * Gets the lastOpenedFile attribute of the JChemPaintPanel object
     *
     *@return The lastOpenedFile value
     */
    public File getLastOpenedFile() {
        return lastOpenedFile;
    }

    /**
     * Sets the lastOpenedFile attribute of the JChemPaintPanel object
     *
     *@param lof
     *            The new lastOpenedFile value
     */
    public void setLastOpenedFile(File lof) {
        this.lastOpenedFile = lof;
    }

    /**
     * Gets the currentOpenFileFilter attribute of the JChemPaintPanel object
     *
     *@return The currentOpenFileFilter value
     */
    public FileFilter getCurrentOpenFileFilter() {
        return currentOpenFileFilter;
    }

    /**
     * Sets the currentOpenFileFilter attribute of the JChemPaintPanel object
     *
     *@param ff
     *            The new currentOpenFileFilter value
     */
    public void setCurrentOpenFileFilter(FileFilter ff) {
        this.currentOpenFileFilter = ff;
    }

    /**
     * Gets the currentSaveFileFilter attribute of the JChemPaintPanel object
     *
     *@return The currentSaveFileFilter value
     */
    public FileFilter getCurrentSaveFileFilter() {
        return currentSaveFileFilter;
    }

    /**
     * Sets the currentSaveFileFilter attribute of the JChemPaintPanel object
     *
     *@param ff
     *            The new currentSaveFileFilter value
     */
    public void setCurrentSaveFileFilter(FileFilter ff) {
        this.currentSaveFileFilter = ff;
    }

    /**
     * Tells if a menu is shown
     *
     *@return The showMenu value
     */
    public boolean getShowMenuBar() {
        return showMenuBar;
    }

    /**
     * Sets if a menu is shown
     *
     *@param showMenuBar
     *            The new showMenuBar value
     */
    public void setShowMenuBar(boolean showMenuBar) {
        this.showMenuBar = showMenuBar;
        customizeView();
    }

    /**
     * Shows and hides menus, statusbar, toolbars according to settings.
     */
    public void customizeView() {
        if (showMenuBar) {
            if (menu == null) {
                menu = new JChemPaintMenuBar(this, this.guistring);
            }
            topContainer.add(menu, BorderLayout.NORTH);
        } else {
            topContainer.remove(menu);
        }
        if (showStatusBar) {
            if (statusBar == null) {
                statusBar = new JCPStatusBar();
            }
            add(statusBar, BorderLayout.SOUTH);
        } else {
            remove(statusBar);
        }
        if (showToolBar) {
            if (uppertoolbar == null) {
                uppertoolbar = JCPToolBar.getToolbar(this, "uppertoolbar", SwingConstants.HORIZONTAL);
            }
            centerContainer.add(uppertoolbar, BorderLayout.NORTH);
            if (lefttoolbar == null) {
            	lefttoolbar = JCPToolBar.getToolbar(this, "lefttoolbar", SwingConstants.VERTICAL);
            }
            centerContainer.add(lefttoolbar, BorderLayout.WEST);
            if (righttoolbar == null) {
            	righttoolbar = JCPToolBar.getToolbar(this, "righttoolbar", SwingConstants.VERTICAL);
            }
            centerContainer.add(righttoolbar, BorderLayout.EAST);
            if (lowertoolbar == null) {
            	lowertoolbar = JCPToolBar.getToolbar(this, "lowertoolbar", SwingConstants.HORIZONTAL);
            }
            centerContainer.add(lowertoolbar, BorderLayout.SOUTH);
        } else {
        	centerContainer.remove(uppertoolbar);
        	centerContainer.remove(lowertoolbar);
        	centerContainer.remove(lefttoolbar);
        	centerContainer.remove(righttoolbar);
        }
        if (showInsertTextField) {
            if (insertTextPanel == null)
                insertTextPanel = new InsertTextPanel(this, null);
            topContainer.add(insertTextPanel, BorderLayout.SOUTH);
        } else {
            topContainer.remove(insertTextPanel);
        }
        revalidate();
    }

    /**
     * Tells if a status bar is shown
     *
     *@return The showStatusBar value
     */
    public boolean getShowStatusBar() {
        return showStatusBar;
    }

    /**
     * Sets the value of showToolbar.
     *
     *@param showToolBar
     *            The value to assign showToolbar.
     */
    public void setShowToolBar(boolean showToolBar) {
        setShowToolBar(showToolBar);
    }


    /**
     * Returns the value of showToolbar.
     *
     *@return The showToolbar value
     */
    public boolean getShowToolBar() {
        return showToolBar;
    }

    /**
     * Sets if statusbar should be shown
     *
     *@param showStatusBar
     *            The value to assign showStatusBar.
     */
    public void setShowStatusBar(boolean showStatusBar) {
        this.showStatusBar = showStatusBar;
        customizeView();
    }

    /**
     * Sets the file currently used for saving this Panel.
     *
     *@param value
     *            The new isAlreadyAFile value
     */
    public void setIsAlreadyAFile(File value) {
        isAlreadyAFile = value;
    }

    /**
     * Returns the file currently used for saving this Panel, null if not yet
     * saved
     *
     *@return The currently used file
     */
    public File isAlreadyAFile() {
        return isAlreadyAFile;
    }

    /**
     * Gets the current gui configuration string of this panel.
     * 
     * @return The current gui configuration string of this panel.
     */
    public String getGuistring() {
        return guistring;
    }

    /**
     * Set to indicate whether the insert text field should be used.
     *
     * @param showInsertTextField
     *            true is the text entry widget is to be shown
     */
    public void setShowInsertTextField(boolean showInsertTextField) {
        this.showInsertTextField = showInsertTextField;
        customizeView();
    }

    /**
     * Tells if the enter text field is currently shown or not.
     *
     * @return text field shown or not
     */
    public boolean getShowInsertTextField() {
        return showInsertTextField;
    }

    /**
     * Gets the SVG of the chemical entities in this panel.
     * 
     * @return The SVG of the chemical entities in this panel.
     */
    public String getSVGString() {
        return this.renderPanel.toSVG();
    }

    /**
     * Takes an image snapshot of this panel.
     * 
     * @return The snapshot.
     */
    public Image takeSnapshot() {
        return this.renderPanel.takeSnapshot();
    }

    /**
     * Shows a warning if the JCPPanel has unsaved content and does save, if the
     * user wants to do it.
     *
     * @return
     *         OptionPane.YES_OPTION/OptionPane.NO_OPTION/OptionPane.CANCEL_OPTION
     */
    public int showWarning() {
        if (isModified && !guistring.equals(JChemPaintEditorApplet.GUI_APPLET)) { 
            int answer = JOptionPane.showConfirmDialog(this, renderPanel
                    .getChemModel().getID()
                    + " " + GT._("has unsaved data. Do you want to save it?"),
                    GT._("Unsaved data"), JOptionPane.YES_NO_CANCEL_OPTION,
                    JOptionPane.WARNING_MESSAGE);
            if (answer == JOptionPane.YES_OPTION) {
                SaveAction saveaction = new SaveAction(this, false);
                saveaction.actionPerformed(new ActionEvent(
                        this, 12, ""));
                if(saveaction.getWasCancelled())
                    answer = JOptionPane.CANCEL_OPTION;
            }
            return answer;
        } else if (guistring.equals(JChemPaintEditorApplet.GUI_APPLET)) {
            return JOptionPane.YES_OPTION;
        } else {
            return JOptionPane.YES_OPTION;
        }
    }

    /**
     * Class for closing jcp
     *
     *@author shk3
     *@cdk.created November 23, 2008
     */
    public final static class AppCloser extends WindowAdapter {

        /**
         * closing Event. Shows a warning if this window has unsaved data and
         * terminates jvm, if last window.
         *
         * @param e
         *            Description of the Parameter
         */
        public void windowClosing(WindowEvent e) {
            int clear = ((JChemPaintPanel) ((JFrame) e.getSource())
                    .getContentPane().getComponents()[0]).showWarning();
            if (JOptionPane.CANCEL_OPTION != clear) {
                for (int i = 0; i < instances.size(); i++) {
                    if (instances.get(i).getTopLevelContainer() == (JFrame) e
                            .getSource()) {
                        instances.remove(i);
                        break;
                    }
                }
                ((JFrame) e.getSource()).setVisible(false);
                ((JFrame) e.getSource()).dispose();
                if (instances.size() == 0) {// TODO &&
                                            // !((JChemPaintPanel)rootFrame.getContentPane().getComponent(0)).isEmbedded())
                                            // {
                    System.exit(0);
                }
            }
        }
    }

    /**
     * Closes all currently opened JCP instances.
     */
    public static void closeAllInstances() {
        int instancesNumber = instances.size();
        for (int i = instancesNumber - 1; i >= 0; i--) {
            JFrame frame = (JFrame) instances.get(i).getTopLevelContainer();
            WindowListener[] wls = (WindowListener[]) (frame
                    .getListeners(WindowListener.class));
            wls[0].windowClosing(new WindowEvent(frame,
                    WindowEvent.WINDOW_CLOSING));
        }
    }

    /* (non-Javadoc)
     * @see org.openscience.jchempaint.controller.IChemModelEventRelayHandler#coordinatesChanged()
     */
    public void coordinatesChanged() {
        setModified(true);
        updateStatusBar();
    }

    /* (non-Javadoc)
     * @see org.openscience.jchempaint.controller.IChemModelEventRelayHandler#selectionChanged()
     */
    public void selectionChanged() {
        updateStatusBar();
        if(this.getRenderPanel().getRenderer().getRenderer2DModel().getSelection()!=null 
        		&& this.getRenderPanel().getRenderer().getRenderer2DModel().getSelection().getConnectedAtomContainer()!=null 
        		&& this.getRenderPanel().getRenderer().getRenderer2DModel().getSelection().getConnectedAtomContainer().getAtomCount()>0)
            enOrDisableMenus(atomMenu,true);
        else
            enOrDisableMenus(atomMenu,false);
        if(this.getRenderPanel().getRenderer().getRenderer2DModel().getSelection()!=null 
        		&& this.getRenderPanel().getRenderer().getRenderer2DModel().getSelection().getConnectedAtomContainer()!=null 
        		&& this.getRenderPanel().getRenderer().getRenderer2DModel().getSelection().getConnectedAtomContainer().getBondCount()>0)
            enOrDisableMenus(bondMenu,true);
        else
            enOrDisableMenus(bondMenu,false);
    }
    
    /**
     * Enables or disables all JMenuItems in a JMenu recursivly.
     * 
     * @param root  The JMenu to search in.
     * @param b     Enable or disable.
     */
    protected void enOrDisableMenus(JMenu root, boolean b) {
        for(int i=0;i<root.getItemCount();i++){
            if(root.getItem(i) instanceof JMenu){
                this.enOrDisableMenus((JMenu)root.getItem(i), b);
            }else if(root.getItem(i) instanceof JMenuItem){
                ((JMenuItem)root.getItem(i)).setEnabled(b);
            }
        }
    }

    /* (non-Javadoc)
     * @see org.openscience.jchempaint.controller.IChemModelEventRelayHandler#structureChanged()
     */
    public void structureChanged() {
        setModified(true);
        updateStatusBar();
        //if something changed in the structure, selection should be cleared
        //this is behaviour like eg in word processors, if you type, selection goes away
        this.getRenderPanel().getRenderer().getRenderer2DModel().setSelection(AbstractSelection.EMPTY_SELECTION);
        updateUndoRedoControls();
        this.get2DHub().updateView();
    }

    /* (non-Javadoc)
     * @see org.openscience.jchempaint.controller.IChemModelEventRelayHandler#structurePropertiesChanged()
     */
    public void structurePropertiesChanged() {
        setModified(true);
        updateStatusBar();
        //if something changed in the structure, selection should be cleared
        //this is behaviour like eg in word processors, if you type, selection goes away
        this.getRenderPanel().getRenderer().getRenderer2DModel().setSelection(AbstractSelection.EMPTY_SELECTION);
    }

    /**
     * Enables/disables the undo/redo button depending on if something can be undone/redone.
     */
    public void updateUndoRedoControls() {
        UndoManager undoManager = renderPanel.getUndoManager();
        JButton redoButton=buttons.get("redo");
        JButton undoButton=buttons.get("undo");
        if (undoManager.canRedo()) {
            redoButton.setEnabled(true);
            redoMenu.setEnabled(true);
            redoButton.setToolTipText(GT._("Redo")+": "+undoManager.getRedoPresentationName());
        } else {
            redoButton.setEnabled(false);
            redoMenu.setEnabled(false);
            redoButton.setToolTipText(GT._("No redo possible"));
        }

        if (undoManager.canUndo()) {
            undoButton.setEnabled(true);
            undoMenu.setEnabled(true);
            undoButton.setToolTipText(GT._("Undo")+": "+undoManager.getUndoPresentationName());
        } else {
            undoButton.setEnabled(false);
            undoMenu.setEnabled(false);
            undoButton.setToolTipText(GT._("No undo possible"));
        }
    }

    /* (non-Javadoc)
     * @see org.openscience.cdk.event.ICDKChangeListener#stateChanged(java.util.EventObject)
     */
    public void stateChanged(EventObject event) {
    	updateUndoRedoControls();
    }

    /* (non-Javadoc)
     * @see java.awt.event.KeyListener#keyPressed(java.awt.event.KeyEvent)
     */
    public void keyPressed(KeyEvent arg0) {
    }

    /* (non-Javadoc)
     * @see java.awt.event.KeyListener#keyReleased(java.awt.event.KeyEvent)
     */
    public void keyReleased(KeyEvent arg0) {
        RendererModel model = renderPanel.getRenderer().getRenderer2DModel();
        ControllerHub relay = renderPanel.getHub();
        if (model.getHighlightedAtom() != null) {
            try {
                IAtom closestAtom = model.getHighlightedAtom();
                char x = arg0.getKeyChar();
                if (Character.isLowerCase(x))
                    x = Character.toUpperCase(x);
                IsotopeFactory ifa;
                ifa = IsotopeFactory.getInstance(closestAtom.getBuilder());
                IIsotope iso = ifa.getMajorIsotope(Character.toString(x));
                if (iso != null) {
                    relay.setSymbol(closestAtom, Character.toString(x));
                }
                this.get2DHub().updateView();
            } catch (IOException e) {
                announceError(e);
            }
        }
    }

    /* (non-Javadoc)
     * @see java.awt.event.KeyListener#keyTyped(java.awt.event.KeyEvent)
     */
    public void keyTyped(KeyEvent arg0) {
    }

    /* (non-Javadoc)
     * @see org.openscience.jchempaint.controller.IChemModelEventRelayHandler#zoomChanged()
     */
    public void zoomChanged() {
        this.updateStatusBar();
    }

    /**
     * Tells if debug output is desired or not.
     *
     * @return debug output or not.
     */
    public boolean isDebug() {
		return debug;
	}

	/* (non-Javadoc)
	 * @see org.openscience.cdk.controller.ChangeModeListener#modeChanged(org.openscience.cdk.controller.IControllerModule)
	 */
	public void modeChanged(IControllerModule newActiveModule) {
	    //we set the old button to inactive colour
        if (this.getLastActionButton() != null)
            this.getLastActionButton().setBackground(JCPToolBar.BUTTON_INACTIVE_COLOR);
        if (this.lastSecondaryButton != null)
            this.lastSecondaryButton.setBackground(JCPToolBar.BUTTON_INACTIVE_COLOR);
        String actionid = newActiveModule.getID();
        //this is because move mode does not have a button
        if(actionid.equals("move"))
            actionid=lastSelectId;
        //we remember the last activated move mode so that we can switch back to it after move
        if(newActiveModule.getID().equals("select") || newActiveModule.getID().equals("lasso"))
            lastSelectId = newActiveModule.getID();
        //we set the new button to active colour
        JButton newActionButton=buttons.get(actionid);
        if(newActionButton!=null){
            this.setLastActionButton(newActionButton);
            newActionButton.setBackground(Color.GRAY);
        }
        if(JCPToolBar.getToolbarResourceString("lefttoolbar", getGuistring()).indexOf(newActiveModule.getID())>-1){
            if(this.buttons.get(this.get2DHub().getController2DModel().getDrawElement())!=null){
                this.buttons.get(this.get2DHub().getController2DModel().getDrawElement()).setBackground(Color.GRAY);
                lastSecondaryButton = this.buttons.get(this.get2DHub().getController2DModel().getDrawElement());
            }else if(buttons.get("periodictable")!=null){
                buttons.get("periodictable").setBackground(Color.GRAY);
                lastSecondaryButton = buttons.get("periodictable");
            }
        }
        if(JCPToolBar.getToolbarResourceString("lowertoolbar", getGuistring()).indexOf(newActiveModule.getID())>-1){
            //the newActiveModule should always be an AddAtomModule, but we still check
            if(newActiveModule instanceof AddAtomModule){
                if(((AddAtomModule)newActiveModule).getStereoForNewBond().equals(IBond.Stereo.NONE)){
                    this.buttons.get("bond").setBackground(Color.GRAY);
                    lastSecondaryButton = this.buttons.get("bond");
                }else if(((AddAtomModule)newActiveModule).getStereoForNewBond().equals(IBond.Stereo.UP)){
                    this.buttons.get("up_bond").setBackground(Color.GRAY);
                    lastSecondaryButton = this.buttons.get("up_bond");
                }else if(((AddAtomModule)newActiveModule).getStereoForNewBond().equals(IBond.Stereo.DOWN)){
                    this.buttons.get("down_bond").setBackground(Color.GRAY);
                    lastSecondaryButton = this.buttons.get("down_bond");
                }else if(((AddAtomModule)newActiveModule).getStereoForNewBond().equals(IBond.Stereo.E_OR_Z)){
                    this.buttons.get("undefined_bond").setBackground(Color.GRAY);
                    lastSecondaryButton = this.buttons.get("undefined_bond");
                }else if(((AddAtomModule)newActiveModule).getStereoForNewBond().equals(IBond.Stereo.UP_OR_DOWN)){
                    this.buttons.get("undefined_stereo_bond").setBackground(Color.GRAY);
                    lastSecondaryButton = this.buttons.get("undefined_stereo_bond");
                }
            }else{
                this.buttons.get("bond").setBackground(Color.GRAY);
                lastSecondaryButton = this.buttons.get("bond");
            }
        }
        if(!(newActiveModule instanceof MoveModule)){
            this.renderPanel.setCursor(new Cursor(Cursor.DEFAULT_CURSOR));
            this.get2DHub().updateView();
        }
        this.updateStatusBar();
	}
	
    /**
     * Gets all atomcontainers of a chemodel in one AtomContainer.
     * 
     * @param chemModel The chemodel
     * @return The result.
     */
    public static IAtomContainer getAllAtomContainersInOne(IChemModel chemModel){
		List<IAtomContainer> acs=ChemModelManipulator.getAllAtomContainers(chemModel);
		IAtomContainer allinone=chemModel.getBuilder().newAtomContainer();
		for(int i=0;i<acs.size();i++){
			allinone.add(acs.get(i));
		}
		return allinone;
    }

    /**
     * Sets the lastSecondaryButton attribute. Only to be used once from JCPToolBar.
     * 
     * @param lastSecondaryButton The lastSecondaryButton.
     */
    public void setLastSecondaryButton(JComponent lastSecondaryButton) {
        this.lastSecondaryButton = lastSecondaryButton;
    }
}

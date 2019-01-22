package org.openscience.cdk.structgen.deterministic;

import java.util.HashMap;
import java.util.Map;

import org.openscience.cdk.Molecule;
import org.openscience.cdk.graph.invariant.MorganNumbersTools;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

public class Fragment extends Molecule {
    private String smiles=null;
    private Map<Long, Map<IAtom,Integer>> attachmentPoints=null;
    long[] morganNumbers = null;
        
    
    public Fragment(){
        super();
    }
    
    public Fragment(Fragment fragment1, Fragment fragment2, long attachpoint1, long attachpoint2){
        super();
        for(int i=0;i<fragment1.atomCount;i++){
            this.addAtom(fragment1.getAtom(i));
        }
        for(int i=0;i<fragment2.atomCount;i++){
            this.addAtom(fragment2.getAtom(i));
        }
        for(int i=0;i<fragment1.bondCount;i++){
            this.addBond(fragment1.getBond(i));
        }
        for(int i=0;i<fragment2.bondCount;i++){
            this.addBond(fragment2.getBond(i));
        }
        this.addBond(this.getBuilder().newBond(fragment1.getAtomByMorgan(attachpoint1), fragment2.getAtomByMorgan(attachpoint2)));
    }
    
    public Fragment(String symbol) {
        super();
        this.addAtom(this.getBuilder().newAtom(symbol));
    }

    private IAtom getAtomByMorgan(long attachpoint1) {
        getAttachmentPoints();// to make sure things are properly initialized
        for(int i=0;i<morganNumbers.length;i++){
            if(morganNumbers[i]==attachpoint1)
                return this.getAtom(i);
        }
        return null;
    }

    public String getSmiles(){
        if(smiles==null){
            smiles = new SmilesGenerator().createSMILES(this);
        }
        return smiles;
    }
    
    public Map<Long, Map<IAtom,Integer>> getAttachmentPoints(){
        if(attachmentPoints==null){
            morganNumbers = MorganNumbersTools.getMorganNumbers(this);
            attachmentPoints = new HashMap<Long,Map<IAtom,Integer>>();
            for(int i=0;i<this.getAtomCount();i++){
                if(this.getAtom(i).getSymbol().equals("C") || this.getAtom(i).getSymbol().equals("H")){
                    int sum = (int) AtomContainerManipulator.getBondOrderSum(this,this.getAtom(i));
                    int desired = 1;
                    if(this.getAtom(i).getSymbol().equals("C"))
                        desired = 4;
                    for(int k=sum;k<desired;k++){
                        if(attachmentPoints.get(morganNumbers[i])==null){
                            Map<IAtom, Integer> map = new HashMap<IAtom, Integer>();
                            attachmentPoints.put(morganNumbers[i], map);
                        }
                        int count = 0;
                        if(attachmentPoints.get(morganNumbers[i]).get(this.getAtom(i))!=null)
                            count = attachmentPoints.get(morganNumbers[i]).get(this.getAtom(i));
                        attachmentPoints.get(morganNumbers[i]).put(this.getAtom(i), count+1);
                    }
                }
            }
        }
        return attachmentPoints;
    }

}

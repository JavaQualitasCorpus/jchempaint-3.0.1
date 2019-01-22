package org.openscience.cdk.structgen.deterministic;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import org.openscience.cdk.Bond;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IIsotope;
import org.openscience.cdk.interfaces.IMolecularFormula;
import org.openscience.cdk.interfaces.IMolecule;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.SaturationChecker;
import org.openscience.cdk.tools.manipulator.MolecularFormulaManipulator;

public class DeterministicGenerator {

    private int level=0;

    public List<String> generateIsomers(String formulaString) throws CloneNotSupportedException, CDKException, IOException, ClassNotFoundException{
        IMolecularFormula formula = DefaultChemObjectBuilder.getInstance().newMolecularFormula();
        return generateIsomers(MolecularFormulaManipulator.getMolecularFormula(formulaString, formula));
    }
    
    
    public List<String> generateIsomers(IMolecularFormula formula) throws CloneNotSupportedException, CDKException, IOException, ClassNotFoundException{
        List<String> solutions=new ArrayList<String>();
        Map<Fragment, Integer> input=new HashMap<Fragment, Integer>();
        for(IIsotope isotope : formula.isotopes()){
            input.put(new Fragment(isotope.getSymbol()),formula.getIsotopeCount(isotope));
        }
        makeNextGeneration(input, solutions);
        return solutions;
    }

    private void makeNextGeneration(Map<Fragment, Integer> input, List<String> solutions) throws CloneNotSupportedException, CDKException, IOException, ClassNotFoundException{
        level++;
        System.err.println("Runde: "+level);
        Iterator<Fragment> it=input.keySet().iterator();
        while(it.hasNext()){
            Fragment x=it.next();
            System.err.println(x.getSmiles()+" "+input.get(x));
        }
        System.err.println(solutions.size());
        if(input.keySet().size()==1 && input.get(input.keySet().iterator().next())==1){
            solutions.add(input.keySet().iterator().next().getSmiles());
        }else{
            //über alle aps iteririeren, mit allen aps joinen
            //jeweils rekursiv aufrufen mit "Rückbau
            List<Map<String, Integer>> result=new ArrayList<Map<String,Integer>>();
            
            String x=input.keySet().iterator().next();
            input.put(x, input.get(x)-1);
            if(input.get(x)==0)
                input.remove(x);
            Iterator<String> it2=((Map<String, Integer>)((HashMap<String, Integer>)input).clone()).keySet().iterator();
            while(it2.hasNext()){
                String c=it2.next();
                System.err.println(c);
                IMolecule acc = new SmilesParser(DefaultChemObjectBuilder.getInstance()).parseSmiles(x);
                for(int i=0;i<acc.getAtomCount();i++){
                    IMolecule mol1=(IMolecule)acc.clone();
                    IMolecule toadd = new SmilesParser(DefaultChemObjectBuilder.getInstance()).parseSmiles(c);
                    for(int k=0;k<toadd.getAtomCount();k++){
                        IMolecule mol2=(IMolecule)toadd.clone();
                        if(isUnsaturated(mol1.getAtom(i), mol1) && isUnsaturated(mol2.getAtom(k),mol2)){
                            IMolecule molnew=mol1.getBuilder().newMolecule();
                            molnew.add(mol1);
                            molnew.add(mol2);
                            molnew.addBond(new Bond(mol1.getAtom(i), mol2.getAtom(k),IBond.Order.SINGLE));
                            String newsmiles = addMoleculeToClasses(molnew, c, x, input);
                            /*for(int firstatom=0;firstatom<molnew.getAtomCount();firstatom++){
                                for(int secondatom=0;secondatom<molnew.getAtomCount();secondatom++){
                                    if(secondatom>=firstatom){*/
                                        makeNextGeneration(input, solutions);

                                        /*try{
                                        }catch(IllegalArgumentException ex){
                                            ex.printStackTrace();
                                            System.out.println("looks like a loop");
                                        }
                                    
                                    }
                                }
                            }*/
                            input.put(newsmiles, input.get(newsmiles)-1);
                            if(input.get(newsmiles)==0)
                                input.remove(newsmiles);
                            if(input.containsKey(c))
                                input.put(c, input.get(c).intValue()+1);
                            else
                                input.put(c, 1);
                            /*if(input.containsKey(x))
                                input.put(x, input.get(x).intValue()+1);
                            else
                                input.put(x, 1);*/                          
                        }
                    }
                }
            }
            if(input.containsKey(x))
                input.put(x, input.get(x).intValue()+1);
            else
                input.put(x, 1);
        }           
        level--;
    }
    
    private boolean isUnsaturated(IAtom atom, IAtomContainer ac) throws CDKException, IOException, ClassNotFoundException {
        if(new SaturationChecker().isSaturated(atom, ac))
            return false;
        else
            return true;
    }

}

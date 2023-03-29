import os
import sys
from PyAutoDock.AutoDockScorer import AutoDock41Scorer
from MolKit import Read
import MolKit.molecule
import MolKit.protein
from AutoDockTools.MoleculePreparation import AD4ReceptorPreparation
from AutoDockTools.MoleculePreparation import AD4LigandPreparation
from AutoDockTools.GridParameters import GridParameters, grid_parameter_list4
from AutoDockTools.GridParameters import GridParameter4FileMaker
from AutoDockTools.atomTypeTools import AutoDock4_AtomTyper
from AutoDockTools.DockingParameters import DockingParameters, DockingParameter4FileMaker, genetic_algorithm_list, \
                genetic_algorithm_local_search_list4, local_search_list4,\
                simulated_annealing_list4
from PyAutoDock.InternalEnergy import InternalEnergy
from PyAutoDock.MolecularSystem import MolecularSystem

if len(sys.argv) == 3:
	pdbqt = sys.argv[1]
	ligqt = sys.argv[2]
	
else:
	print ("AD_score.py <prot> <lig>")
	sys.exit()
	
def autodock_scoring(receptor, ligand):
    receptorfilename =  receptor
    ligandfilename =  ligand
    write_file_mode = False
    parameter_library_filename = None
    exclude_torsFreeEnergy = False
    verbose = None
    ad_scorer = AutoDock41Scorer(exclude_torsFreeEnergy=exclude_torsFreeEnergy)
    supported_types = ad_scorer.supported_types
    receptor = Read(receptorfilename)[0]
    receptor.buildBondsByDistance()
    ligand = Read(ligandfilename)[0]
    ligand.buildBondsByDistance()

    ms = MolecularSystem()
    ms.add_entities(receptor.allAtoms)
    ms.add_entities(ligand.allAtoms)
    ad_scorer.set_molecular_system(ms)
    #get the scores, score per term:
    [estat, hb, vdw ,dsolv] = ad_scorer.get_score_per_term()
    torsEnrg = ligand.TORSDOF * ad_scorer.tors_weight
    score = estat +hb +vdw 
    #print score,estat,hb,vdw,dsolv
    output_score = "%8.4f %6.4f %6.4f %6.4f %6.4f %6.4f" %(score, estat, hb, vdw, dsolv, torsEnrg)
    return score
    #print output_score.strip().split(" ")[0]
    #return float(output_score.strip().split(" ")[0])	
    #return "-".join(map(str, [ score,estat,hb,vdw,dsolv,torsEnrg]))
try:
	ad_val = autodock_scoring(pdbqt, ligqt)
except:
	ad_val = 1000
	pass

print(ad_val)





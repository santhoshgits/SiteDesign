import warnings
warnings.simplefilter("ignore")


import time
import numpy as np
import os
import shutil
import sys
import random
import copy
#from multiprocess import Pool
from multiprocessing import Pool
from CustomClasses.SiteGenIterate import SiteGenIter
from CustomClasses.Dock import Dock
from CustomClasses.GiveCoord import GiveCoord 
#from Fitness4 import Fitness
from CustomClasses.Fitness6 import Fitness
from CustomClasses.Mutation3 import Mutation
from CustomClasses.CrossOver1 import CrossOver
#from CrossOver2 import CrossOver


if len(sys.argv) == 6:
	lig = sys.argv[1]
	nprocs = int(sys.argv[2])
	outFold = sys.argv[3]
	minim = int(sys.argv[4])
	maxim = int(sys.argv[5])
else:
	print('\n\n')
	print("INSTRUCTION TO EXECUTE THE CODE")
	print ("SiteDesignGenerator.py <lig.pdb> <no.of.procs> <Folder> <Min> <Max>")
	print('\n\n')
	sys.exit()
	
	

	
dire = os.getcwd()

site_gen_iter = SiteGenIter()
print ('Loaded SiteGenIter Class')
give_coord = GiveCoord(outFold)
print ('Loaded GiveCoord Class')
dock = Dock()
print ('Loaded AutoDock Class')
fitness = Fitness(outFold)
print ('Loaded Fitness Class')
cross_over = CrossOver(outFold)
print ('Loaded Cross-Over Class')
mutation = Mutation(outFold)
print ('Loaded Mutation Class')

# export PATH=/home/santhosh/mgltools_x86_64Linux2_1.5.6/bin:$PATH

# export PATH=/media/santhosh/7add6b54-492d-474d-b478-871833f96d09/home/santhosh/Downloads/mgltools_x86_64Linux2_1.5.6/bin:$PATH

os.chdir(dire+'/'+outFold)

os.system('python2.5 prepare_ligand4.py -l lig.pdb -o lig.pdbqt')
if os.path.isfile(dire+'/'+outFold+"/temp/lig.pdbqt"):
	os.system('rm '+dire+'/'+outFold+"/temp/lig.pdbqt")
shutil.move(dire+'/'+outFold+"/lig.pdbqt", dire+'/'+outFold+"/temp/")
shutil.copy(dire+'/'+'CustomClasses'+'/prepare_receptor4.py', dire+'/'+outFold+'/temp')
shutil.copy(dire+'/'+'CustomClasses'+'/AD_score.py', dire+'/'+outFold+'/temp')

os.chdir(dire)

def automaton_dock(Total):
	nos, val = Total
	if os.path.isfile(dire+'/'+outFold+"/temp/"+str(nos)+".pdbqt"):
		os.system('rm '+dire+'/'+outFold+"/temp/"+str(nos)+".pdbqt")
	if os.path.isfile(dire+'/'+outFold+"/temp/"+str(nos)+".pdb"):
		os.system('rm '+dire+'/'+outFold+"/temp/"+str(nos)+".pdb")	
	return dock.file_read(val, nos, outFold)


'''
def automaton_cross_over(Sites, Nos):
	SitesTotal = copy.deepcopy(Sites)
	for _ in range(Nos):
		site1 = random.choice(Sites)
		site2 = random.choice(Sites)

		site1 = give_coord.file_read(site1).strip().split('\n')
		site2 = give_coord.file_read(site2).strip().split('\n')
		#print site1
		site_new = cross_over.file_read(site1, site2)

		SitesTotal.append(site_new)
	return SitesTotal	
'''

def automaton_cross_over(Total):
	site1, site2 = Total
	site1 = give_coord.file_read(site1).strip().split('\n')
	site2 = give_coord.file_read(site2).strip().split('\n')
	site_new = cross_over.file_read(site1, site2, outFold)
	return site_new
	
def automaton_mutation(Total):
	Sites = Total
	Nos = 1
	NewSites = []
	for i in Sites:
		NewSites.append(i)
	aline = give_coord.file_read(NewSites).strip().split('\n')
	#print(aline)
	new_site = mutation.file_read(aline, 25, outFold)
	return new_site


	
'''	
def automaton_mutation(Sites, Nos):
	NewSites = []
	for i in Sites:
		NewSites.append(i)
	for _ in range(Nos):
		site = random.choice(Sites)
		aline = give_coord.file_read(site).strip().split('\n')
		new_site = mutation.file_read(aline, 25)
		NewSites.append(new_site)
	return NewSites	
'''




def sort_fitness(Nos, Nos1, fit_arr):
	arr = []
	for line in fit_arr:
		line = line.strip()
		l = line.split('\t')
		arr.append(l)
	#print(arr)	
	ans = fitness.file_read(arr, outFold)
	SitesNew = []
	count_dic = {}
	for i in ans[:Nos]:
		i1 = i.split('\t')
		pdb = int(i1[0].split('.')[0])
		count_dic[pdb] = 0
		#print pdb
		SitesNew.append(arr[pdb][0].split('_'))
		
	SitesRandom = []
	for i in range(0, len(arr)):
		if i not in count_dic:
			SitesRandom.append( arr[i][0].split('_') )
		
	for _ in range(Nos1):
		if SitesRandom:
			SitesNew.append(random.choice(SitesRandom))
	return SitesNew	
		
		
	

ga_run = 10



res_point_line = open(dire+'/'+outFold+'/residue_position.txt', 'r').readlines()	
ligand_aline = open(dire+'/'+outFold+'/'+lig, 'r').readlines()	


p = Pool(nprocs)



out = open(dire+'/'+outFold+'/PocketDesignOutput.txt', 'w')
for ga_count in range(ga_run):
	print (ga_count, 'GA_count')


	automaton_initial = site_gen_iter.file_read(res_point_line, ligand_aline, minim, maxim, 20, nprocs, outFold)
	
	initial_sites = []
	#print(automaton_initial)
	for i in automaton_initial:
		initial_sites.append(i[0])

	total_sites = []
	for i in initial_sites:
		#print (i)
		total_sites.append(give_coord.file_read(i))	


	automaton_dock_array = []
	for i in range(len(total_sites)):
		automaton_dock_array.append( [ i, total_sites[i] ] )
	
	dock_val = p.map(automaton_dock, automaton_dock_array)
	dock_collate = []
	print('Dock Score Completed')
	#print(len(dock_val))
	out1 = open(dire+'/'+outFold+'/dock_val.txt', 'w')
	for i in range(len(dock_val)):
		#print(dock_val[i])
		dock_collate.append('_'.join(initial_sites[i])+'\t'+str(dock_val[i])+'\n')
		out1.write('_'.join(initial_sites[i])+'\t'+str(dock_val[i])+'\n')
	out1.close()


	#sys.exit()
	NewSites = sort_fitness(10, 5, dock_collate)	
	#print(NewSites)
	# population gen starts from here
	while_count = 0
	
	automaton_extend = site_gen_iter.file_read(res_point_line, ligand_aline, minim, maxim, 25, nprocs, outFold)
	automaton_extend = [ i[0] for i in automaton_extend ]
	aut_count = 0
	
	while while_count < 5:
		while_count += 1
		print (ga_count,' - ',while_count)
		co_inp = []
		for _ in range(50):
			co_inp.append([ random.choice(NewSites), random.choice(NewSites) ])
		#co_inp.extend(NewSites)	
		
		cross_over_sites = p.map(automaton_cross_over, co_inp)
		cross_over_sites = [ i for i in cross_over_sites ]
		cross_over_sites.extend(NewSites)
		
		cross_over_sites.extend(automaton_extend[aut_count:aut_count+5])
		
		total_sites = []
		for i in cross_over_sites:
			total_sites.append(give_coord.file_read(i))
		automaton_dock_array = []
		for i in range(len(total_sites)):
			automaton_dock_array.append( [ i, total_sites[i] ] )
		dock_val = p.map(automaton_dock, automaton_dock_array)
		
		dock_collate = []
		for i in range(len(dock_val)):
			dock_collate.append('_'.join(cross_over_sites[i])+'\t'+str(dock_val[i])+'\n')
		NewSites = sort_fitness(10, 5, dock_collate)

		cross_over_top = []
		mut_inp = []
		for i in NewSites:
			cross_over_top.append(i)
			#mut_inp.append(i)
		
		
		for _ in range(50):
			mut_inp.append( random.choice(NewSites) )
		#print(mut_inp,'----')	
		print ('entering mutation')
		mutation_sites = p.map(automaton_mutation, mut_inp)
		mutation_sites.extend(NewSites)
		
		mutation_sites.extend(automaton_extend[aut_count:aut_count+5])
		
		total_sites = []
		for i in mutation_sites:
			total_sites.append(give_coord.file_read(i))
		
		automaton_dock_array = []
		for i in range(len(total_sites)):
			automaton_dock_array.append( [ i, total_sites[i] ] )
		dock_val = p.map(automaton_dock, automaton_dock_array)
		
		dock_collate = []
		for i in range(len(dock_val)):
			dock_collate.append('_'.join(mutation_sites[i])+'\t'+str(dock_val[i])+'\n')
			
			
		NewSites = sort_fitness(10, 5, dock_collate)
		
		
		for i in NewSites:
			out.write(str(ga_count)+'\t'+str(while_count)+'\t'+'_'.join(i)+'\n')
		
		#for i in cross_over_top:
		#	NewSites.append(i)
		print (len(NewSites))	
		print ('\n')
		
		aut_count += 5
	
	

	#time.sleep(11)


out.close()



time.sleep(3)

print ('Generating Sites in PDB')
os.system('python3.8 '+dire+'/CustomClasses/GiveSites.py '+outFold)
print ('PDB files are generated. folder "site" contains the *.pdb files\n')



out = open(dire+'/'+outFold+'/Pairs.txt', 'w')

for i in os.listdir(dire+'/'+outFold+'/site'):
	for j in os.listdir(dire+'/'+outFold+'/site'):
		out.write(i+'\t'+j+'\n')
out.close()

time.sleep(1)

os.system('python '+dire+'/CustomClasses/CreateVectors.py site SiteVec '+outFold)

time.sleep(1)

print('Comparing Generated Sites among them for grouping')
os.system('python '+dire+'/CustomClasses/PocketApproxAlignParallel.py SiteVec Pairs.txt SiteSim.txt '+str(nprocs)+' '+outFold)
time.sleep(1)


print('Aligning Designed Sites Against Known PDB sites')
os.system('python '+dire+'/CustomClasses/NRScanSelf.py '+outFold+' '+str(nprocs)+' '+outFold+' '+lig[:-4])

#print('Aligning Designed Sites Against Non-Redundant PDBs')
#os.system('python '+dire+'/CustomClasses/NRScan.py '+outFold+' '+str(nprocs)+' '+outFold)







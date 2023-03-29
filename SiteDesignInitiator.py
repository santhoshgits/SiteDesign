import warnings
warnings.simplefilter("ignore")

import os
import sys
from collections import Counter, defaultdict
import shutil
import time
import random
import copy
import numpy as np


from multiprocessing import Pool
from CustomClasses.MCSS4 import MCSS
from CustomClasses.AlignSites1 import AlignSites
from CustomClasses.Rename_site import RenameSite
from CustomClasses.Residue_Outlier import ResidueOutlier
from CustomClasses.Rename_site_whole import RenameSiteWhole
from CustomClasses.SortLigHits2 import SortLigHits
from CustomClasses.ResidueFrequency import ResidueFrequency


#from CustomClasses.ClusterResidues import Clustering
#from CustomClasses.ClusterPeptide1 import ClusterPeptide

from CustomClasses.RenameSecond import RenameSecond
from CustomClasses.ResidueStretch import ResidueStretch
from CustomClasses.ResPoints1 import ResPoints

from CustomClasses.StrideMapper import StrideMapper
from CustomClasses.GeometryAndCoordinate import GeometryCoordinate

MCSS = MCSS()
align_sites = AlignSites()
rename_site = RenameSite()
residue_outlier = ResidueOutlier()
rename_site_whole = RenameSiteWhole()
sort_lig_hits = SortLigHits()
residue_frequency = ResidueFrequency()
#cluster = Clustering()	
#cluster_energy = ClusterEnergy()


residue_stretch = ResidueStretch()
rename_second = RenameSecond()	
stride_mapper = StrideMapper()
geomCoord = GeometryCoordinate()


if len(sys.argv) == 4:
	lig = sys.argv[1]
	nprocs = int(sys.argv[2])
	outFold = sys.argv[3]
else:
	print('\n\n')
	print("INSTRUCTION TO EXECUTE THE CODE")
	print('For Example: To run SiteInitiator.py for ATP ligand, the syntax would be like the below.')
	print('python3.8 PocketDesignInitiator.py ATP.pdb 4 ATPOutput')
	print('Here: 4 indicates total no. of processor that will be utilised for the computation and ATPOutput is the folder onto which all outputs will be pasted')
	print('NOTE: Please remember the folder name: ATPOutput. This folder needs to be specified as the argument for the subsequent python codes')
	print('\n')
	
	sys.exit()
	
dire = os.getcwd()

path_out_fold = dire+'/'+outFold
if not os.path.exists(path_out_fold):
	os.mkdir(path_out_fold)
else:
	shutil.rmtree(path_out_fold)
	os.mkdir(path_out_fold)	
	
shutil.copy(dire+'/CustomDatas/ResWt', dire+'/'+outFold)


path_out_temp = dire+'/'+outFold+'/temp'

if not os.path.exists(path_out_temp):
	os.mkdir(path_out_temp)
	shutil.copy(dire+"/checkmol", path_out_temp)
else:
	shutil.rmtree(path_out_temp)
	os.mkdir(path_out_temp)	
	shutil.copy(dire+"/checkmol", path_out_temp)


path_out_matched_hits = dire+'/'+outFold+'/matched_hits'
if not os.path.exists(path_out_matched_hits):
	os.mkdir(path_out_matched_hits)
else:
	shutil.rmtree(path_out_matched_hits)
	os.mkdir(path_out_matched_hits)		


path_out_rename = dire+'/'+outFold+'/matched_hits_renamed'
if not os.path.exists(path_out_rename):
	os.mkdir(path_out_rename)
else:
	shutil.rmtree(path_out_rename)
	os.mkdir(path_out_rename)		


lig_dir = dire+'/CustomDatas/nrsite'


lig_aline = open(lig, 'r').readlines()
lig_aline = [ line for line in lig_aline if line[:6] == 'HETATM' and line[73:].strip() != 'H']


out = open(dire+'/'+outFold+'/lig.pdb', 'w')
for line in lig_aline:
	out.write(line)
out.close()	


shutil.copy(dire+'/'+lig, dire+'/'+outFold)



res_dict = {'GLY':1 ,'ALA':2 ,'LEU':2 ,'ILE':2 ,'TRP':5 ,'SER':2 ,'THR':2 ,'TYR':4 ,'PHE':4 ,\
					 'PRO':2 ,'ASP':4,'GLU':4 ,'HIS':4 ,'CYS':3 ,'MET':3 ,'VAL':2 ,'ASN':3 ,\
					 'GLN':3 ,'ARG':7 ,'LYS':7 }




def fn_lig_association():
	aline = open(dire+'/CustomDatas/NrLigandCluster', 'r').readlines()
	dic = defaultdict(list)
	ligand_cluster_dic = defaultdict(list)
	for line in aline:
		line = line.strip()
		#line = line.replace('\t', ' ')
		l = line.split(' ')
		#print l
		'''
		l = sorted(set(l))
		for i in l:
			#dic[i].append(j)
			for j in l:
				dic[i].append(j)
		'''
		for i in l[1:]:
			dic[i[0]].append(i)		
	return dic			
	
lig_dic_id = fn_lig_association()		


def mcss(line):
	line, count = line.split(' ')
	aline = open(lig_dir+"/"+line, 'r').readlines()
	aline = [ line for line in aline if line[:6] == 'HETATM' and line[73:].strip() != 'H']
	#print aline
	#print lig_aline
	#print count
	matches = None
	try:
		matches = MCSS.file_read(aline, lig_aline, count, outFold)
	except:
		pass

	return matches

lig_dic_track = {}
lig_arr = os.listdir(lig_dir)	

os.system('python3.8 '+dire+'/CustomClasses/MCSSNumba.py '+lig+' '+outFold+' '+str(nprocs))


p = Pool(nprocs)



time.sleep(1)
aline = open(dire+'/'+outFold+'/BestSite.txt','r').readlines()
BestSiteArr = []
for line in aline:
	line = line.strip()
	l = line.split(' ')
	if int(l[1]) > 0:
		BestSiteArr.append(l[0])
	
BestSiteArr = sorted(set(BestSiteArr))
L = len(BestSiteArr)
Step = 100
out = open(dire+'/'+outFold+'/ligand_mcss_hits.txt', 'w')
for i in range(0, L, Step*nprocs):
	count = 0
	arr = []
	print('Total No. of Alignable Ligand is {:0.0f}. Site that have been processed so far is {:0.0f}. Percent completed : {:3.3f}%'.format(L, i, (i/float(L))*100 ), end = '\r')
	for j in BestSiteArr[i:i+Step*nprocs]:
		arr.append(j+" "+str(count))
		count += 1
	match1 = p.map(mcss, arr)	
	for i1 in range(len(arr)):
		if match1[i1]:
			if match1[i1] == 'SANT':
				#print match1[i1], arr[i1]
				lig = arr[i1].split(' ')[0].split('_')[1]

			else:
				lig = arr[i1].split(' ')[0].split('_')[1]
				for j in range(len(match1[i1])):
					out.write(arr[i1].split(' ')[0]+" "+match1[i1][j]+"\n")

out.close()



time.sleep(3)


print ('\nstarting python class SortLigHits')
sort_lig_hits.file_read(lig_aline, outFold)
print ('class SortLigHits completed successfully\n')
time.sleep(2)


print ('starting python class AlignSites')
align_sites.file_read(lig_aline, outFold)
print ('class AlignSites completed successfully\n')
time.sleep(2)

print ('starting python class Rename_site')
rename_site.file_read(outFold)	
print ('class Rename_site completed successfully\n')
time.sleep(2)


print('Running Stride Mapper')
stride_mapper.file_read(outFold)
time.sleep(2)



print ('starting python class ResidueFrequency')
residue_frequency.file_read(outFold)
print ('class ResidueFrequencycompleted successfully\n')
time.sleep(2)

#print 'done'
#time.sleep(11)
out = open(dire+'/'+outFold+"/collated.pdb", 'w')
for i in os.listdir(dire+'/'+outFold+"/matched_hits_renamed"):
	aline = open(dire+'/'+outFold+"/matched_hits_renamed/"+i, 'r').readlines()
	out.write("".join(aline))
out.close()		
time.sleep(2)

print ('starting python class Residue_Outlier')
residue_outlier.file_read("collated.pdb", lig_aline, outFold)	
print ('class Residue_Outlier completed successfully\n')
time.sleep(2)

print ('starting python class  Rename_site_whole')
rename_site_whole.file_read(outFold)
print ('class Rename_site_whole completed successfully\n')
time.sleep(2)


print ('starting python class ClusterResidues')
#cluster.file_read("final_residue_pool.pdb")
os.system('python3.8 '+dire+'/CustomClasses/ClusterResiduesPar.py final_residue_pool.pdb '+str(nprocs)+' '+outFold)
print ('class ClusterResidues completed successfully\n')
time.sleep(2)


#shutil.copy(dire+'/residue_cluster', dire+'/residue_cluster1')

print ('starting python class ClusterEnergy1')
os.system('python3.8 '+dire+'/CustomClasses/ClusterEnergy2.py lig.pdb '+str(nprocs)+' '+outFold)
print ('class ClusterEnergy1 completed successfully\n')
time.sleep(2)
print ('----cluster done')


shutil.copy(dire+"/CustomClasses/prepare_receptor4.py", dire+'/'+outFold+"/temp/")
print ('----- receptor prepared')
shutil.copy(dire+"/CustomClasses/AD_score.py", dire+'/'+outFold+"/temp/")
print ('------ad score calculated')
shutil.copy(dire+'/'+outFold+'/residue_cluster1.txt', dire+'/'+outFold+'/residue_cluster2.txt')
time.sleep(2)

'''
print ('starting python class ProtLigAngle1')
prot_lig_angle.file_read(outFold)
print ('class ProtLigAngle1 completed succesfully\n')
time.sleep(2)

print ('starting python class ProtProtAngle1')
prot_angle.file_read(outFold)
print ('class ProtProtAngle1 completed sucessfully\n')
time.sleep(2)
'''

print ('starting python class RenameSecond')
rename_second.file_read(outFold)
print ('class RenameSecond completed successfully')


print ('starting python class ResidueStretch')
residue_stretch.file_read(lig_aline, outFold)
print ('class ResidueStretch completed successfully')
#print 'done'

time.sleep(2)


# ----- kind of new module recently written




time.sleep(2)
#print (dire)
out = open(dire+'/'+outFold+'/final.pdb', 'w')
for i in os.listdir(dire+'/'+outFold+'/matched_hits_renamed'):
	aline = open(dire+'/'+outFold+'/matched_hits_renamed/'+i, 'r').readlines()
	for line in aline:
		line = line.strip()
		if line[:4] == "ATOM" and line[17:20] in res_dict:
			out.write(line+'\n')
out.close()


time.sleep(2)




aline = open(dire+'/'+outFold+'/final.pdb', 'r').readlines()
# 11:15 17:26
dic = {}
out = open(dire+'/'+outFold+'/final2.pdb', 'w')
for line in aline:
	var = line[11:16]+line[17:26]
	if var not in dic:
		dic[var] = 0
		out.write(line)
out.close()	


time.sleep(2)


print('Running Geometry and Coordinate')
geomCoord.file_read(outFold)

time.sleep(2)


res_points = ResPoints(lig_aline, outFold)
print ('Starting Res Points')
res_points.file_read(outFold)
print ('ResPoints class completed\n')

time.sleep(2)
print ('Starting K-means clustering')
os.system('python3.8 '+dire+'/CustomClasses/KMCluster.py '+outFold)
print ('K-means clustering completed\n')
time.sleep(2)

print ('Starting Residue Association Matrix')
#os.system('python3.8 '+dire+'/CustomClasses/ResidueAssociationMatrix1.py '+outFold)
os.system('python3.8 '+dire+'/CustomClasses/ResidueAssociationMatrix2.py '+outFold)
print ('Resdiue Matrix build completed\n')

time.sleep(2)





	

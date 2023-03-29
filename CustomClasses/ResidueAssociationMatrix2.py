import time
import copy
import numpy as np
import os
import sys
from collections import Counter, defaultdict
import math
from numpy import save
from numba import jit


dire = os.getcwd()

'''
dic = {'GLY':'G','ALA':'A','VAL':'V','LEU':'L','ILE':'I','THR':'T','SER':'S','TYR':'Y','TRP':'W',\
	  'PRO':'P','PHE':'F','ASN':'N','GLN':'Q','ASP':'D','GLU':'E','ARG':'R','LYS':'K','HIS':'H',\
	   'CYS':'C','MET':'M'
	  }
'''

if len(sys.argv) == 2:
	Folder = sys.argv[1]
else:
	print('python3.8 ResidueAssociationMatrix2.py <Folder>')
	sys.exit()
		

ResArr = ['GLY', 'ALA', 'VAL', 'LEU', 'ILE', 'THR', 'SER', 'TYR', 'TRP',
'PRO', 'PHE', 'ASN', 'GLN', 'ASP', 'GLU', 'ARG', 'LYS', 'HIS', 'CYS', 'MET']



def GroupResidues(aline, query_lig, target_lig, het_mat):
	
	lig_dic = {}
	for i in range(len(query_lig)):
		lig_dic[target_lig[i]] = query_lig[i]
	dic = defaultdict(list)
	
	het_arr = []
	#print het_mat
	for i in aline:
		i1 = i.strip()
		if i1[:4] == 'ATOM' and i1[13:16].strip() == 'CA':
			dic[i1[17:26]].append([ i1[28:38].strip(), i1[38:46].strip(), i1[46:54].strip() ])
		elif i1[:6] == 'HETATM':
			het_arr.append([ i1[6:11].strip(), i1[28:38].strip(), i1[38:46].strip(), i1[46:54].strip() ])
	
	prot_arr = []
	for i,j in dic.items():
		prot_arr.append([i,np.asarray(j, dtype='float').mean(axis=0)])
		
	ResMatrix = []	
	ResMatrix_dic = {}
	
	het_res_dist_dic = {}
	query_count = {}
	c = 0
	for i in het_mat:
		query_count[i] = c
		het_res_dist_dic[i] = np.zeros((20,10), dtype='int')
		c += 1

	for i in het_arr:
		#print(i, lig_dic[i[0]])
		
		x, y, z = list(map(float,i[1:]))
		het_dic = {}
		for i1 in ResArr:
			het_dic[i1] = 0
			
		mat = []
		for i1 in range(len(ResArr)):
			Ans = [0]*10
			index = 0
			for j in prot_arr:
				if ResArr[i1] == j[0][:3]:
					x1, y1, z1 = j[1]
					x_ans = pow(( x - x1),2)
					y_ans = pow(( y - y1),2)
					z_ans = pow(( z - z1),2)
					ans = math.sqrt(x_ans + y_ans + z_ans)
					if index < 9:
						Ans[index] = int(ans)
					index += 1
			#print(Ans)		
			mat.append(Ans)	
					
		het_res_dist_dic[lig_dic[i[0]]] = mat
			
	#print(het_res_dist_dic.keys())		
	matrix = []
	for i in het_mat:
		matrix.append(het_res_dist_dic[i])
	matrix = np.asarray(matrix, dtype='int')			
		
	return matrix
	

aline = open(dire+'/'+Folder+'/align_site_mapper.txt', 'r').readlines()
align = copy.deepcopy(aline)
bline = open(dire+'/'+Folder+'/lig.pdb','r').readlines()
het_mat = [ i1[6:11].strip() for i1 in bline if i1[:6] == 'HETATM' ]


MatrixTotal = []
for line in aline:
	line = line.strip()
	l = line.split(' ')
	#if l[0] != '1KPI_SAH_A_900.pdb':
	#	continue
	#print(l)
	site_aline = open(dire+'/'+Folder+'/matched_hits/'+l[-1], 'r').readlines()
	ResMatrix_dic = GroupResidues(site_aline, l[1].split('-'), l[2].split('-'), copy.deepcopy(het_mat))
	MatrixTotal.append(ResMatrix_dic)
	
MatrixTotal = np.asarray(MatrixTotal, dtype='int')	
save(dire+'/'+Folder+'/Association2.npy',MatrixTotal)
	
print(MatrixTotal.shape)






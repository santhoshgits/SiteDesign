import time
import os
import math
import numpy as np
from collections import Counter, defaultdict
import re

dire = os.getcwd()

class SiteGenFit():

	def blosum(self):
		residue_dict_single = {'G':'GLY','A':'ALA','V':'VAL','L':'LEU','I':'ILE','T':'THR','S':'SER','Y':'TYR','W':'TRP',\
	  'P':'PRO','F':'PHE','N':'ASN','Q':'GLN','D':'ASP','E':'GLU','R':'ARG','K':'LYS','H':'HIS',\
	   'C':'CYS','M':'MET'
	  }
		aline = []
		aline.append("A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *")
		aline.append("A  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1  0 -4")
		aline.append("R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1  0 -1 -4")
		aline.append("N -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3  3  0 -1 -4")
		aline.append("D -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3  4  1 -1 -4")
		aline.append("C  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -3 -2 -4")
		aline.append("Q -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2  0  3 -1 -4")
		aline.append("E -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4")
		aline.append("G  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 -1 -2 -1 -4")
		aline.append("H -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3  0  0 -1 -4")
		aline.append("I -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3 -3 -3 -1 -4")
		aline.append("L -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1 -4 -3 -1 -4")
		aline.append("K -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2  0  1 -1 -4")
		aline.append("M -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1 -3 -1 -1 -4")
		aline.append("F -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1 -3 -3 -1 -4")
		aline.append("P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2 -2 -1 -2 -4")
		aline.append("S  1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2  0  0  0 -4")
		aline.append("T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0 -1 -1  0 -4")
		aline.append("W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3 -4 -3 -2 -4")
		aline.append("Y -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1 -3 -2 -1 -4")
		aline.append("V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4 -3 -2 -1 -4")
		res_info = aline[0].strip()
		res_info = re.sub(" {1,}"," ", res_info)
		res_info_split = res_info.split(" ")[:-4]
		ans = []
		for line in aline[1:]:
			line = line.strip()
			line = re.sub(" {1,}"," ",line)
			residue = line.split(" ")[0]
			min_max = []
			for i in line.split(" ")[1:-4]:
				min_max.append(int(i))
			minimum = min(min_max)
			#print minimum
			#if minimum < 0:
			new_min = minimum
			#new_min = abs(minimum)
			new_min_max = []
			for i in min_max:
				new_min_max.append(i+new_min)
			maxima = float(max(new_min_max))

			for j,k in zip(res_info_split, new_min_max):
				ans.append(residue_dict_single[residue]+" "+residue_dict_single[j]+" "+str(k))

		dic_temp = {}
		residue_pairs_dictionary = {}
		for i in ans:
			i0 = i.split(" ")[0]
			i1 = i.split(" ")[1]
			val = int(i.split(" ")[2])
			if i0 not in dic_temp:
				dic_temp[i0] = 0
				residue_pairs_dictionary[i0] = {}
			residue_pairs_dictionary[i0][i1] = val		
		return residue_pairs_dictionary
	
	
	
	def get_cn(self, res, res_arr):
		arr = []
		x, y, z = np.asarray(res_coord[res], dtype='float').mean(axis=0)
		for i in res_arr:
			x1, y1, z1 = np.asarray(res_coord[i], dtype='float').mean(axis=0)
			x_ans = pow(( x - x1 ),2)
			y_ans = pow(( y - y1 ),2)
			z_ans = pow(( z - z1 ),2)
			ans = math.sqrt(x_ans + y_ans + z_ans)
			if ans < 1.5:
				arr.append(i)
		sum1 = float(sum([ NewGridVal['_'.join(map(str, map(int, map(float, res_coord[i][1]))))] for i in arr ] ))
		sum2 = float(sum([ NewGridVal['_'.join(map(str, map(int, map(float, res_coord[i][1]))))] for i in res_arr ]) )
		diff = sum2 - float(sum1)		
		return (sum([ residue_pairs_dictionary[res[:3]][i[:3]] for i in arr ]) * (len(res_arr)/float(len(arr)))), diff
		
		
	def get_ca(self, res_coord1):
		final, final1 = [], []
		for i in res_coord1.items():
			x, y, z = map(float, i[1][1])
			arr = []
			for j in res_coord.items():
				x1, y1, z1 = map(float, j[1][1])
				x_ans = pow(( x - x1 ),2)
				y_ans = pow(( y - y1 ),2)
				z_ans = pow(( z - z1 ),2)
				ans = math.sqrt(x_ans + y_ans + z_ans)
				if ans < 1:
					arr.append(j[0])
			#print len(arr)		
			arr = sorted(set(arr))	
			#print sum([ NewGridVal['_'.join(map(str, map(int, map(float, res_coord[i][1]))))] for i in arr ])		
			ans, ans1 = self.get_cn(i[0], arr)
			#print ans, ans1
			final1.append(ans1)
			final.append(ans)

		return sum(final), sum(final1)
		
		
	def file_read(self, res_arr, Folder):
		global residue_pairs_dictionary, res_coord, NewGridVal
		
		'''
		res_count = {'GLY':1 ,'ALA':2 ,'LEU':2 ,'ILE':2 ,'TRP':9 ,'SER':3 ,'THR':3 ,'TYR':9 ,'PHE':9 ,\
						 'PRO':3 ,'ASP':6,'GLU':6 ,'HIS':4 ,'CYS':4 ,'MET':4 ,'VAL':2 ,'ASN':3 ,\
						 'GLN':3 ,'ARG':12 ,'LYS':10 }	
		'''
		res_count = {}
		aline = open(dire+'/'+Folder+'/ResWt','r').readlines()
		for line in aline:
			line = line.strip()
			l = line.split(' ')
			res_count[l[0]] = int(l[1])		 
						 
		aline = open(dire+'/'+Folder+'/final2.pdb','r').readlines()
		res_coord = defaultdict(list)
		res_line = defaultdict(list)
		TotalGrids = []
		GridRes = defaultdict(list)
		for line in aline:
			line = line.strip()
			res_coord[line[17:26]].append([line[28:38].strip(), line[38:46].strip(), line[46:54].strip()])	
			res_line[line[17:26]].append(line)
			if line[:4] == 'ATOM' and line[13:16].strip() == 'CA':
				coord = [line[28:38].strip(), line[38:46].strip(), line[46:54].strip()]
				coord = map(float, coord)
				coord = map(int, coord)
				coord = "_".join(map(str, coord))
				TotalGrids.append(coord)
				#for _ in range(res_count[line[17:20]]):
				GridRes[coord].append(line[17:20])
				
		NewGridVal = {}		
		for i in Counter(TotalGrids).items():
			val = i[1] + sum([ res_count[j] for j in  GridRes[i[0]] ])
			NewGridVal[i[0]] = val	
			
		'''	
		res_count = {'GLY':1 ,'ALA':2 ,'LEU':2 ,'ILE':2 ,'TRP':5 ,'SER':2 ,'THR':2 ,'TYR':4 ,'PHE':4 ,\
			 'PRO':2 ,'ASP':4,'GLU':4 ,'HIS':4 ,'CYS':3 ,'MET':3 ,'VAL':2 ,'ASN':3 ,\
			 'GLN':3 ,'ARG':7 ,'LYS':7 }
		'''

		aline = open(dire+'/'+Folder+'/reside_rename_map.txt', 'r').readlines()			
		rename = {}
		for line in aline:
			line = line.strip()
			l = line.split('\t')
			rename[l[1]] = l[0]
		aline = open(dire+'/'+Folder+'/residue_cluster1.txt', 'r').readlines()	
		res_energy = {}
		for line in aline:
			line = line.strip()
			l = line.split('-')
			if len(l) >= 1:
				for i in l:
					if i in rename:
						res_energy[rename[i]] = 0		
		residue_pairs_dictionary = self.blosum()
		
		arr = []
		for i in res_arr:
			for j in res_line[i]:
				arr.append(j)
		
		res_coord1 = defaultdict(list)	
		for line in arr:
			line = line.strip()
			res_coord1[line[17:26]].append([line[28:38].strip(), line[38:46].strip(), line[46:54].strip()])	
		ans, ans1 = self.get_ca(res_coord1)
		return ans/float(len(res_coord1)), ans1/len(res_coord1)
		
		

				

'''
site_gen_fit = SiteGenFit()
aline = open(dire+'/temp1/'+'0a.pdb', 'r').readlines()
arr = sorted(set([ i[17:26] for i in aline if i[:4] == 'ATOM']))
for _ in range(10):
	print site_gen_fit.file_read(arr)
'''





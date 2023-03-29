import os
from collections import defaultdict, Counter
import time
import random
import numpy as np
import math
import random
import shutil
#from BuildStretcha import BuildStretch
from CustomClasses.BuildStretch16 import BuildStretch
from numba import jit



dire = os.getcwd()


@jit(nopython=True)
def compare(C1, C2):

	#print(C1, C2)
	T = [[0,1], [0,2], [0,3], [1,0], [1,2], [1,3], [2,0], [2,1], [2,3], [3,0], [3,1], [3,2]]
	
	cn_ans = float(1000)
	for n in [[0, 2],[2, 0]]:
		#print(C1[n[0]])
		ans = math.sqrt( (pow(( C1[n[0]][0] - C2[n[1]][0] ),2)) + (pow(( C1[n[0]][1] - C2[n[1]][1] ),2)) + (pow(( C1[n[0]][2] - C2[n[1]][2] ),2)) )  
		cn_ans = min(cn_ans, ans)
		if ans < 1.19:
			#print('tru1')
			return True
	
	#print(cn_ans)
	for i in range(len(C1)):
		for j in range(len(C2)):
			ans = math.sqrt( (pow(( C1[i][0] - C2[j][0] ),2)) + (pow(( C1[i][1] - C2[j][1] ),2)) + (pow(( C1[i][2] - C2[j][2] ),2))   )
			#print(i,j, ans)
			
			if [i,j] not in T:
				#print('not in T')
				if ans < 2.45:
					#print('true3')
					return True
		
			
			if [i,j] in T:

				if [i,j] not in [[0,2],[2,0]]:
					if ans < 2.11: 
						#print('true4')
						return True
					if cn_ans < 3.0:
						if ans-.1 < cn_ans:
							#print('true5', ans, cn_ans)
							return True
					
				else:
					if ans < 1.19:
						return True
					if 2.5 > ans > 1.8:
						return True
					
	return False



class Mutation():

	def __init__(self, Folder):
		res_point_line = open(dire+'/'+Folder+'/residue_position.txt', 'r').readlines()	
		ligand_aline = open(dire+'/'+Folder+'/lig.pdb', 'r').readlines()	
		self.build_stretch = BuildStretch(res_point_line, ligand_aline, Folder)
		self.res_id = []
		self.res_coord, self.res_line = defaultdict(list), defaultdict(list)
		self.res_dict = {'GLY':'G' ,'ALA':'A' ,'LEU':'L' ,'ILE':'I' ,'TRP':'W' ,'SER':'S' ,'THR':'T' ,'TYR':'Y' ,'PHE':'F' ,\
						 'PRO':'P' ,'ASP':'D','GLU':'E' ,'HIS':'H' ,'CYS':'C' ,'MET':'M' ,'VAL':'V' ,'ASN':'N' ,\
						 'GLN':'Q' ,'ARG':'R' ,'LYS':'K' }
						 
		'''				 
		self.res_count = {'GLY':1 ,'ALA':2 ,'LEU':2 ,'ILE':2 ,'TRP':9 ,'SER':3 ,'THR':3 ,'TYR':9 ,'PHE':9 ,\
						 'PRO':3 ,'ASP':6,'GLU':6 ,'HIS':6 ,'CYS':4 ,'MET':4 ,'VAL':2 ,'ASN':3 ,\
						 'GLN':3 ,'ARG':12 ,'LYS':9 }
		'''
		
		self.res_count = {}			 
		aline = open(dire+'/'+Folder+'/ResWt','r').readlines()	 
		for line in aline:
			line = line.strip()
			l = line.split(' ')
			self.res_count[l[0]] = int(l[1])			 	
						 
		aline = open(dire+'/'+Folder+'/final2.pdb', 'r').readlines()
		TotalGrids = []
		GridRes = defaultdict(list)
		for line in aline:
			line = line.strip()
			if line[:4] == 'ATOM' and line[13:16].strip() == 'CA':
				coord = [line[28:38].strip(), line[38:46].strip(), line[46:54].strip()]
				coord = map(float, coord)
				coord = map(int, coord)
				coord = "_".join(map(str, coord))
				GridRes[coord].append(line[17:20])
				TotalGrids.append(coord)
			
		NewGridVal = {}			
		for i in Counter(TotalGrids).items():
			val = i[1] * sum([ self.res_count[j] for j in  GridRes[i[0]] ])
			NewGridVal[i[0]] = val		
				
				
		aline = open(dire+'/'+Folder+'/residue_position.txt', 'r').readlines()		
		dic = defaultdict(list)
		for line in aline:
			line = line.strip()
			l = line.split('\t')
			r1 = l[1].split('-')
			r2 = l[3].split(' ')
			ln1 = len(r1) 
			#print(r1, [ i[:3] for i in r1 ])
			val = sum([ self.res_count[i[:3]] for i in r1 ])*(len(r1)**2)
			val1 = sum([ NewGridVal[i] for i in r2 if i in NewGridVal  ])
			for i in zip(r1, r2):
				dic[i[1]].append([i[0], val1])
			
			
		
		self.GridToResidues = defaultdict(list)	
		for i,j in dic.items():
			dic1 = {}	
			for i1 in j:
				if i1[0] not in dic1:
					dic1[i1[0]] = 0
				else:
					val = dic1[i1[0]]
					val += i1[1]
					dic1[i1[0]] = val
			for i1 in dic1.items():
				self.GridToResidues[i].append([i1[0],i1[1]])	
					
				
		aline = open(dire+'/'+Folder+'/final2.pdb','r').readlines()		
		self.res_coord_dict = defaultdict(list)
		self.res_line = defaultdict(list)

		res_id = []
		for line in aline:
			line = line.strip()
			#if line[17:26] not in ResidueNotConsider:
			if line[:4] == 'ATOM':
				self.res_coord_dict[line[17:26]].append([  float(line[27:38].strip()), float(line[38:46].strip()), float(line[46:54].strip()) ])
				#self.res_line[line[17:26]].append(line)
				res_id.append(line[17:26])
		self.res_id = sorted(set(res_id))		
				
		self.res_coord = {}
		for i,j in self.res_coord_dict.items():
			self.res_coord[i] = np.asarray(j, dtype='float')	
				
		
		del dic
		del dic1
		del NewGridVal
		del TotalGrids
		del GridRes
		
	
	def check_clash(self, res_taken, res_inp):	
		for i in res_taken:
			if compare( self.res_coord[i], self.res_coord[res_inp] ):
				return True, [res_inp, i]
		return False, None
		
		
	'''	
	def check_clash(self, res_taken, res_inp):	
		#print res_inp
		#print res_taken
		res_coord = np.asarray(self.res_coord_dict[res_inp], dtype='float')
		dic_pep = {"0 2":0, '2 0':0}
		#dic = {"0 2":0, '2 0':0, '0 3':0, '3 0':0, '0 1':0, '1 0':0, '1 2':0, '2 1':0}
		dic = {}
		
		res1_n_c = res_coord[0]
		res1_n_c = np.append([res1_n_c], [res_coord[2]], axis=0)
		c_n = [[0,1], [1,0]]
		
		nos = [0,1,2,3]
		for i in nos:
			for j in nos:
				if i != j:
					val = str(i)+" "+str(j)
					dic[val] = 0
		dic_ca = {}
		for i in range(len(res_coord)):
			if i not in dic_ca:
				x, y, z = res_coord[i]
				for j in res_taken:
					#print j
					#print res_coord_dic[j], j
					res_coord1 = np.asarray(self.res_coord_dict[j], dtype='float')
					res2_n_c = res_coord1[0]
					res2_n_c = np.append([res2_n_c], [res_coord1[2]], axis=0)
					
					cn_ans = []
					for n in c_n:
						xn, yn, zn = res1_n_c[n[0]]
						xc, yc, zc = res2_n_c[n[1]]
						x_ans = pow(( xn - xc ),2)
						y_ans = pow(( yn - yc ),2)
						z_ans = pow(( zn - zc ),2)
						ans = math.sqrt(x_ans + y_ans + z_ans)
						cn_ans.append(ans)
					#print cn_ans	
					cn_ans = min(cn_ans)
					#print cn_ans
					if cn_ans < 1.15:
						return True
				
					for j1 in range(len(res_coord1)):
						if j1 not in dic_ca:
							val = str(i)+" "+str(j1)
							#if j1 not in dic and i not in dic:
							x1, y1, z1 = res_coord1[j1]
							x_ans = pow(( x - x1),2)
							y_ans = pow(( y - y1),2)
							z_ans = pow(( z - z1),2)
							ans = math.sqrt(x_ans + y_ans + z_ans)

							if val not in dic:
								if ans < 2.45:
									return True

							if val in dic:
								if val not in dic_pep:			
									if ans < 2:
										return True
									if cn_ans < 3.2:
										if ans < cn_ans:
											return True
								else:
									if ans < 1.15:
										return True
									if 2.5 > ans > 1.8:
										return True
		return False
	'''
	
	
	
	def mutate_residues(self, res_id_taken, res_id):
		grid = '_'.join(map(str,map(int,map(float, self.res_coord_dict[res_id_taken][1]))))
		#print grid
		grid_res1 = [ random.choice(self.GridToResidues[grid]) for _ in range(10) ]
		random.shuffle(grid_res1)
		grid_res, dic =[], {}
		for i in grid_res1:
			if i[0] not in dic:
				dic[i[0]] = 0
				grid_res.append(i)
		grid_res = sorted(grid_res, key = lambda x:int(x[1]), reverse=True)
		for i in grid_res:
			#print i,'--'
			if not self.check_clash(res_id, i[0]):
				return i[0], True
				#break
		return None, False
	
	
	
	
	def file_read(self, aline, rate, Folder):
		
		ligand_aline = open(dire+'/'+Folder+'/lig.pdb', 'r').readlines()
		
		res_id = sorted(set( [ line[17:26] for line in aline if line[:4] == 'ATOM' ] ))
		rates = int((float(rate)/100)*len(res_id))
		random.shuffle(res_id)
		taken_res = []
		#print rates
		for i in range(rates):
			taken_res.append(res_id[i])
		res_new = [ i for i in res_id if i not in taken_res ]
		
		for i in taken_res:
			#print i
			i1, check = self.mutate_residues(i, res_new)
			if check:
				res_new.append(i1)
		if len(res_new) < len(res_id):
			res_new = self.build_stretch.file_read(1, ligand_aline, res_new, len(res_id))
	
		return res_new
		
	
'''				
dire = os.getcwd()
aline = open(dire+'/temp1/0a.pdb', 'r').readlines()
arr = sorted(set([ line[17:26] for line in aline if line[13:16].strip() == 'CA' ]))
		
				
mutation = Mutation()			
arr1 = mutation.file_read(arr, 25)				
print (arr)
print (arr1)

dire = os.getcwd()
aline = open('final2.pdb', 'r').readlines()
dic = defaultdict(list)
for line in aline:
	line = line.strip()
	dic[line[17:26]].append(line)
	
out = open('a.pdb', 'w')	
for i in arr:
	for line in dic[i]:
		out.write(line+'\n')
out.close()
		
out = open('b.pdb', 'w')	
for i in arr1:
	for line in dic[i]:
		out.write(line+'\n')
out.close()				
'''						 








import math
import numpy as np
from collections import defaultdict
import time
import random
#from BuildStretcha import BuildStretch
from CustomClasses.BuildStretch16 import BuildStretch
import os
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
	
	
	
class CrossOver:

	
	def __init__(self, Folder):
		'''
		self.res_count = {'GLY':1 ,'ALA':2 ,'LEU':2 ,'ILE':2 ,'TRP':6 ,'SER':3 ,'THR':3 ,'TYR':4 ,'PHE':4 ,\
					 'PRO':2 ,'ASP':4,'GLU':4 ,'HIS':3 ,'CYS':3 ,'MET':3 ,'VAL':2 ,'ASN':4 ,\
					 'GLN':4 ,'ARG':6 ,'LYS':5 }
		'''
		
		self.res_count = {}			 
		aline = open(dire+'/'+Folder+'/ResWt','r').readlines()	 
		for line in aline:
			line = line.strip()
			l = line.split(' ')
			self.res_count[l[0]] = int(l[1])
					 
		res_point_line = open(dire+'/'+Folder+'/residue_position.txt', 'r').readlines()	
		ligand_aline = open(dire+'/'+Folder+'/lig.pdb', 'r').readlines()	
		self.build_stretch = BuildStretch(res_point_line, ligand_aline, Folder)			 
		
		
		aline = open(dire+'/'+Folder+'/final2.pdb', 'r').readlines()
		self.res_id = []
		res_coord_dic, self.res_line = defaultdict(list), defaultdict(list)
		for line in aline:
			line = line.strip()
			if line[:4] == 'ATOM':
				res_coord_dic[line[17:26]].append([line[28:38].strip(), line[38:46].strip(), line[46:54].strip()])	
				#self.res_line[line[17:26]].append(line)	
				self.res_id.append(line[17:26])
				
		self.res_coord = {}
		for i,j in res_coord_dic.items():
			self.res_coord[i] = np.asarray(j, dtype='float')
		
		
		
		'''
		aline = open('final2.pdb', 'r').readlines()
		self.res_line =  defaultdict(list)
		for line in aline:
			line = line.strip()	
			self.res_line[line[17:26]].append(line)	
		'''
	
	
	
	def coord_rep(self, site1, site2):
		'''
		site1, site2 = [], []
		for i in site1a:
			for j in self.res_line[i]:
				site1.append(j)
		for i in site2a:
			for j in self.res_line[i]:
				site2.append(j)		
		'''		
		
		res_coord_dic = defaultdict(list)
		res_id = []
		dic = {}
		for line in site1:
			line = line.strip()
			res_id.append(line[17:26])
			dic[line[17:26]] = 0
			res_coord_dic[line[17:26]].append([line[28:38].strip(), line[38:46].strip(), line[46:54].strip()])
		ln = len(sorted(set(res_id)))	
		for line in site2:
			line = line.strip()
			if line[17:26] not in dic:
				#print line[17:26]
				res_id.append(line[17:26])
				res_coord_dic[line[17:26]].append([line[28:38].strip(), line[38:46].strip(), line[46:54].strip()])	
		res_id = sorted(set(res_id))	
		return res_id, res_coord_dic, ln
	
	
	def check_clash(self, res_taken, res_inp):	
		for i in res_taken:
			if compare( self.res_coord[i], self.res_coord[res_inp] ):
				return True, [res_inp, i]
		return False, None
	
	
	'''
	def check_clash(self, res_taken, res_inp, res_coord_dic):	
		#print res_inp
		#print res_taken
		res_coord = np.asarray(res_coord_dic[res_inp], dtype='float')
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
					res_coord1 = np.asarray(res_coord_dic[j], dtype='float')
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
					if cn_ans < 1.19:
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
									if ans < 2.11:
										return True
									if cn_ans < 3.0:
										if ans-.1 < cn_ans:
											return True
								else:
									if ans < 1.19:
										return True
									if 2.5 > ans > 1.8:
										return True

		return False
	'''
	
	
	
	def GetOverlap(self, res_taken, res_inp, res_coord_dic):
		
		res_coord = np.asarray(res_coord_dic[res_inp], dtype='float')
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
		res_overlap = []
		for i in range(len(res_coord)):
			if i not in dic_ca:
				x, y, z = res_coord[i]
				for j in res_taken:
					#print j
					res_coord1 = np.asarray(res_coord_dic[j], dtype='float')
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
						res_overlap.append([res_inp, j])
				
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
									res_overlap.append([res_inp, j])

							if val in dic:
								if val not in dic_pep:			
									if ans < 2.11:
										res_overlap.append([res_inp, j])
									if cn_ans < 3.0:
										if ans-.1 < cn_ans:
											res_overlap.append([res_inp, j])
								else:
									if ans < 1.19:
										res_overlap.append([res_inp, j])
									if 2.5 > ans > 1.8:
										res_overlap.append([res_inp, j])
		if res_overlap:								
			r1, r2 = res_overlap[0]
			#print r1,r2
			arr = []
			for _ in range(0, self.res_count[r1[:3]]):
				arr.append(r1)
			for _ in range(0, self.res_count[r2[:3]]):
				arr.append(r2)
			return random.choice(arr)	
		else:
			return res_inp
		
	
	
	
	
	def cross_sites(self, res_id, res_coord_dic, rng):
		random.shuffle(res_id)
		arr = []
		for i in res_id:
			if len(arr) == rng:
				break
			if not arr:
				arr.append(i)
			else:
				#time.sleep(11)
				if not self.check_clash(arr, i):
					arr.append(i)
		return arr	
		
	
	def file_read(self, site1, site2, Folder):
		#res_id1, res_coord_dic1 = self.coord_rep(site1)
		#res_id2, res_coord_dic2 = self.coord_rep(site2)
		
		res_id, res_coord_dic, rng = self.coord_rep(site1, site2)
		#self.region_overlap(res_id1, res_coord_dic1, res_id2, res_coord_dic2)
		random.shuffle(res_id)
		#print(res_id)
		res_total = []
		dic_track = {}
		for i in res_id:
			arr = []
			for j in res_id:
				if i != j:
					if self.check_clash([j], i):
						arr.append(i)
						arr.append(j)
			arr.append(i)			
			arr = sorted(set(arr))
			arr1 = []
			for k in arr:	
				#print(k,'--')
				for _ in range(0, self.res_count[k[:3]]):
					arr1.append(k)
			#if arr1:	
			val = random.choice(arr1)
			if val not in dic_track:
				res_total.append(val)
				dic_track[val] = 0
		for i in res_id:
			if i not in dic_track:
				res_total.append(i)
				dic_track[i] = 0
		
		#self.GetOverlap(res_taken, res_inp)
		arr = self.cross_sites(res_total, res_coord_dic, rng)
		ligand_aline = open(dire+'/'+Folder+'/lig.pdb', 'r').readlines()
		if len(arr) < rng:
			#print ('co culprit', rng, len(arr))
			arr = self.build_stretch.file_read(1, ligand_aline, arr, rng)
		
		return arr
	
'''
cross_over = CrossOver()
dire = os.getcwd()
site2 = open(dire+'/temp1/0a.pdb', 'r').readlines()
#site1 = sorted(set([ line[17:26] for line in aline if line[13:16].strip() == 'CA' ]))

site1 = open(dire+'/temp1/1a.pdb', 'r').readlines()
#site2 = sorted(set([ line[17:26] for line in aline if line[13:16].strip() == 'CA' ]))
site_arr = cross_over.file_read(site1, site2)
print (site_arr)
'''

'''
dic = defaultdict(list)
aline = open('final.pdb', 'r').readlines()
for line in aline:
	dic[line[17:26]].append(line.strip())

out = open('b0.pdb', 'w')
for i in site_arr:
	print i
	for j in dic[i]:
		out.write(j+'\n')
out.close()
'''


	






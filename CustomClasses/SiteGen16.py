import time
import numpy as np
import os
from collections import Counter, defaultdict
import random
import copy
import math
import sys
import re
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
	
	

class SiteGen():
	def __init__(self, res_point_line, ligand_aline, Folder):
		self.lig_atom_id = sorted(set([i[6:11].strip() for i in ligand_aline if i[:6] =='HETATM']))
		self.res_dict = {'GLY':'G' ,'ALA':'A' ,'LEU':'L' ,'ILE':'I' ,'TRP':'W' ,'SER':'S' ,'THR':'T' ,'TYR':'Y' ,'PHE':'F' ,\
						 'PRO':'P' ,'ASP':'D','GLU':'E' ,'HIS':'H' ,'CYS':'C' ,'MET':'M' ,'VAL':'V' ,'ASN':'N' ,\
						 'GLN':'Q' ,'ARG':'R' ,'LYS':'K' }
		
		
		self.res_count1 = {}
		self.res_count_stretch ={}	 
		aline = open(dire+'/'+Folder+'/ResWt','r').readlines()	 
		for line in aline:
			line = line.strip()
			l = line.split(' ')
			self.res_count1[l[0]] = int(l[1])
			self.res_count_stretch[l[0]] = int(l[1])
		self.res_count = self.res_count1
		self.res_lig_points = defaultdict(list)
		
	
	
		aline = open(dire+'/'+Folder+'/final2.pdb', 'r').readlines()
		self.res_id = []
		res_coord_dic, self.res_line = defaultdict(list), defaultdict(list)
		for line in aline:
			line = line.strip()
			if line[17:20] in self.res_dict:
				res_coord_dic[line[17:26]].append([line[28:38].strip(), line[38:46].strip(), line[46:54].strip()])	
				#self.res_line[line[17:26]].append(line)	
				self.res_id.append(line[17:26])
				
		self.res_coord = {}
		for i,j in res_coord_dic.items():
			self.res_coord[i] = np.asarray(j, dtype='float')
				
				
				
		self.res_id = sorted(set(self.res_id))
		
		aline = ligand_aline
		self.het_coord = defaultdict(list)
		for line in aline:
			if line[:6] == 'HETATM':
				self.het_coord[line[6:11].strip()].append(float(line[27:38].strip()))
				self.het_coord[line[6:11].strip()].append(float(line[38:46].strip()))
				self.het_coord[line[6:11].strip()].append(float(line[46:54].strip()))
		
		
		aline = open(dire+'/'+Folder+'/final2.pdb','r').readlines()
		TotalGrids = []
		GridRes = defaultdict(list)
		for line in aline:
			line = line.strip()
			if line[:4] == 'ATOM' and line[13:16].strip() == 'CA':
				coord = [line[28:38].strip(), line[38:46].strip(), line[46:54].strip()]
				coord = map(float, coord)
				coord = map(int, coord)
				coord = "_".join(map(str, coord))
				TotalGrids.append(coord)
				GridRes[coord].append(line[17:20])

		
		NewGridVal = {}	
		for i in Counter(TotalGrids).items():
			val = sum([ self.res_count[j] for j in  GridRes[i[0]] ]) # modified
			NewGridVal[i[0]] = val

		self.GridToRes = defaultdict(list)
		for i in Counter(GridRes).items():
			for j in Counter(i[1]).items():
				self.GridToRes[i[0]].append([j[0], j[1]])
			
	
		self.NewGridVal = NewGridVal 		
		self.res_lig_points = defaultdict(list)
		self.residue_gap_fill = [] # avoid early stopping
		for res_ids in self.lig_atom_id:
			for line in res_point_line:
				line = line.strip()
				l = line.split('\t')
				if l[0] == res_ids:
					#print l
					val = sum([ NewGridVal[i] for i in l[3].split(' ') if i in NewGridVal ]) #modified
					self.res_lig_points[res_ids].append([l, val])
					#print res_ids, l, val
					self.residue_gap_fill.append([ l[0], l[1],val])
					
		aline = open(dire+'/'+Folder+'/reside_rename_map.txt', 'r').readlines()			
		rename = {}
		for line in aline:
			line = line.strip()
			l = line.split('\t')
			rename[l[1]] = l[0]
		
		aline = open(dire+'/'+Folder+'/residue_cluster1.txt', 'r').readlines()	
		self.res_energy = {}
		for line in aline:
			line = line.strip()
			l = line.split('-')
			if len(l) >= 1:
				for i in l:
					#if i in rename:
					self.res_energy[rename[i]] = 0
					
		
		self.associated_grid_residues = defaultdict(list)			
		for line in res_point_line: # can incorprate length criteria
			line = line.strip()
			l = line.split('\t')
			r1 = l[1].split('-')
			r2 = l[3].split(' ')
			ln1 = len(r1) 
			#print ln1,l
			#time.sleep(1)
			for i in range(len(r1)):
				if r2[i] in NewGridVal:
					res = NewGridVal[r2[i]]
					self.associated_grid_residues[r2[i]].append( [ '_'.join(r1), ln1*res ] ) # modified		
						
	
	def check_clash(self, res_taken, res_inp):	
		#print (res_inp,'res inp')
		#print (res_taken)
		res_coord = np.asarray(self.res_coord[res_inp], dtype='float')
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
					#print j,' is j', res_taken
					res_coord1 = np.asarray(self.res_coord[j], dtype='float')
					#print j, len(res_coord1), res_coord1
					
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
						return True, [res_inp, j]
				
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
									return True, [res_inp, j]

							if val in dic:
								if val not in dic_pep:			
									if ans < 2.11: # can change it to 2.0
										return True, [res_inp, j]
									if cn_ans < 3.0:
										if ans-.1 < cn_ans:
											return True, [res_inp, j]
								else:
									if ans < 1.19:
										return True, [res_inp, j]
									if 2.5 > ans > 1.8:
										return True, [res_inp, j]

		return False, None
		
	
	
	def FillGaps(self, res_taken):
		#print 'Enter Fill'
		L = len(res_taken)
		count = 0
		dic = {}
		arr = []
		while count < 10000:
			count += 1
			val = random.choice(self.residue_gap_fill)
			arr.append(val)
			
		arr1 = sorted(arr, key = lambda x:int(x[2]), reverse=True)
		
		arr = arr1[:1000]
		#for _ in range(300):
		#	arr.append(random.choice(arr1))
		random.shuffle(arr)
		
		arr1 = []
		for i in arr:
			if i[0] not in dic:
				dic[i[0]] = 0
			if i[0] in dic:
				count = dic[i[0]]
				count += 1
				dic[i[0]] = count
			if dic[i[0]] < 100:
				arr1.append(i)
		#print len(arr), len(arr1)
		
		
		ca_coord = np.asarray([ self.res_coord[i][1] for i in res_taken ], dtype='float')
		cn_coord = np.asarray([ self.res_coord[i][0] for i in res_taken ], dtype='float')
		cc_coord = np.asarray([ self.res_coord[i][2] for i in res_taken ], dtype='float')
		cnt_coord = np.asarray([ np.asarray(self.res_coord[i], dtype='float').mean(axis=0) for i in res_taken ], dtype='float')
		pep_arr = []
		for i in range(len(res_taken)):
			x, y, z = cc_coord[i]
			for j in range(len(res_taken)):
				if i != j:
					x1, y1, z1 = cn_coord[j]
					x_ans = pow(( x - x1),2)
					y_ans = pow(( y - y1),2)
					z_ans = pow(( z - z1),2)
					ans = math.sqrt(x_ans + y_ans + z_ans)
					if ans < 2.4:
						pep_arr.append(res_taken[i])
						pep_arr.append(res_taken[j])
		pep_arr = sorted(set(pep_arr))
		
		
		
		for res in arr1[:100]:
			#break
			for i in res[1].split('-'):	
				if i in self.res_energy:
					count = 0
					for j in range(len(res_taken)):
						ClashCheck, ClashPairs = self.check_clash([res_taken[j]], i)
						if ClashCheck:
							count += 1
					if count == 1:
						add_list = []
						for j in range(len(res_taken)):
							ClashCheck_temp, ClashPairs_temp = self.check_clash([res_taken[j]], i)
							if ClashCheck_temp:
								ClashCheck, ClashPairs = ClashCheck_temp, ClashPairs_temp
						
						#if ClashPairs[] not in pep_arr:
						if ClashCheck and ClashPairs[1] not in pep_arr :
							#print ClashPairs, i, pep_arr
							r1 = '_'.join(map(str,map(int,map(float, self.res_coord[i][1]))))
							r2 = '_'.join(map(str,map(int,map(float, self.res_coord[ ClashPairs[1] ][1]))))
							#print r1, r2
							#print self.NewGridVal[r1], self.NewGridVal[r2]
							val1 = [ self.res_count[j[0]]*j[1]*self.NewGridVal[r1] for j in self.GridToRes[r1] if j[0] == i[:3] ][0]
							val2 = [ self.res_count[j[0]]*j[1]*self.NewGridVal[r2] for j in self.GridToRes[r2] if j[0] == ClashPairs[1][:3] ][0]
							
							if val2 > val1:
								if random.randint(1,10) > 2:
									add_list.append(ClashPairs[1])
								else:
									add_list.append(i)	
							else:
								if random.randint(1,10) < 5:
									add_list.append(i)
										
							res_taken2 = add_list + res_taken
							#print add_list
							res_taken1 = []
							for j in res_taken2:
								ClashCheck, ClashPairs = self.check_clash(res_taken1, i)
								if not ClashCheck:
									res_taken1.append(i)	
							#print len(res_taken1), len(res_taken) , res_taken
							if len(res_taken1) > len(res_taken):	
								res_taken = res_taken1
			
		
		count = 0
		while count < 300:
			count += 1
			res = random.choice(arr1)
			
			for i in res[1].split('-'):	
				#print i,'==='
				ClashCheck, ClashPairs = self.check_clash(res_taken, i)
				add_list = []
				if not ClashCheck:
					res_taken.append(i)
					if len(res_taken) > L+2:
						count = 1000
						#print L, len(res_taken)
						break
				
					
				if len(res_taken) > L+2:
					break
									
			
		return res_taken	
		
		
	
	def FindNewRes(self, res_taken):
		grids = [ '_'.join(map(str,map(int,map(float, self.res_coord[i][1])))) for i in res_taken ]
		grids = sorted(set(grids))
		random.shuffle(grids)
		#print '\n'
		
		ca_coord = np.asarray([ self.res_coord[i][1] for i in res_taken ], dtype='float')
		cn_coord = np.asarray([ self.res_coord[i][0] for i in res_taken ], dtype='float')
		cc_coord = np.asarray([ self.res_coord[i][2] for i in res_taken ], dtype='float')
		cnt_coord = np.asarray([ np.asarray(self.res_coord[i], dtype='float').mean(axis=0) for i in res_taken ], dtype='float')
		pep_arr = []
		for i in range(len(res_taken)):
			x, y, z = cc_coord[i]
			for j in range(len(res_taken)):
				if i != j:
					x1, y1, z1 = cn_coord[j]
					x_ans = pow(( x - x1),2)
					y_ans = pow(( y - y1),2)
					z_ans = pow(( z - z1),2)
					ans = math.sqrt(x_ans + y_ans + z_ans)
					if ans < 2.4:
						pep_arr.append(res_taken[i])
						pep_arr.append(res_taken[j])
		pep_arr = sorted(set(pep_arr))
		
		
		picked_res = []
		
		L = len(res_taken)
		check = False
		picked_count = 0
		for i in grids:
			check = False
			
			#print i
			#PickedGrid = [ random.choice(self.associated_grid_residues[i]) for _ in range(100) ]
			PickedGrid = []
			for _ in range(1000):
				val = random.choice(self.associated_grid_residues[i])
				if val not in PickedGrid:
					PickedGrid.append(val)

			PickedGrid1 = sorted(PickedGrid, key = lambda x:int(x[1]), reverse=True)
			PickedGrid = PickedGrid1[:100]
			#for _ in range(300):
			#	PickedGrid.append(random.choice(PickedGrid1))
			random.shuffle(PickedGrid)	
							
			#print PickedGrid
			for j in PickedGrid:
				#print j
				if check:
					break
				count = 0
				for j1 in j[0].split('_'):
					ClashCheck, ClashPairs = self.check_clash(res_taken, j1)
					if ClashCheck:
						count += 1
				#if count == 1:
				#print j
				if count <= 1:
					
					#print j,'got it'
					picked_res.append(j[0].split('_'))
					picked_count += 1
					if picked_count > 10:
						check = True
						break
					
			#print '\n\n'
			#time.sleep(1)
	
		
		#picked_res = ['ARG A2000', 'CYS A 257', 'SER A1639', 'GLY A3531']
		#print res_taken
		#print picked_res
		random.shuffle(picked_res)
		L = len(res_taken)
		if picked_res:
			#return res_taken
			#print 'if', L
			for ii in picked_res:
				add_list = []
				for i in ii:
					ClashCheck, ClashPairs = self.check_clash(res_taken, i)
					if not ClashCheck:
						res_taken.append(i)
					else:
						if ClashPairs[1] not in pep_arr:
							#print ClashPairs,'---',i
							#print self.res_coord[i],'\n\n\n'
							#print self.res_coord[ClashPairs[1]]
							r1 = '_'.join(map(str,map(int,map(float, self.res_coord[i][1]))))
							r2 = '_'.join(map(str,map(int,map(float, self.res_coord[ ClashPairs[1] ][1]))))
							#print r1, r2
							#print self.NewGridVal[r1], self.NewGridVal[r2]
							val1 = [ self.res_count[j[0]]*j[1]*self.NewGridVal[r1] for j in self.GridToRes[r1] if j[0] == i[:3] ][0]
							val2 = [ self.res_count[j[0]]*j[1]*self.NewGridVal[r2] for j in self.GridToRes[r2] if j[0] == ClashPairs[1][:3] ][0]
							if val2 > val1:
								add_list.append(ClashPairs[1])
							else:
								if random.randint(1,10) < 5:
									add_list.append(i)
								
				res_taken2 = add_list + res_taken
				#print add_list
				res_taken1 = []
				for i in res_taken2:
					ClashCheck, ClashPairs = self.check_clash(res_taken1, i)
					if not ClashCheck:
						res_taken1.append(i)
				#print len(res_taken1), add_list		
				
				if len(res_taken1) > len(res_taken):
					#print 'hi'
					res_taken = res_taken1
					return res_taken
				
				#cc = 10
				#while cc = 0:
				#	cc -= 1
				if len(res_taken1) > L+1:
					#print 'react 1', len(res_taken1), L 
					return res_taken1
					
				
				
			#print 'out'
			
			if len(res_taken1) <= len(res_taken):
				#print 'fill'
				res_taken1 = self.FillGaps(res_taken1)	
				return res_taken1
					
				
								
		else:		
			#print '  ELSE  '
			res_taken = self.FillGaps(res_taken)	
			return res_taken
		
		#return res_taken		
			
			
			
	def NoMapRes1(self, res_taken):	
		
		random.shuffle(res_taken)
		#print 'NomapRes1'
		ca_coord = np.asarray([ self.res_coord[i][1] for i in res_taken ], dtype='float')
		cn_coord = np.asarray([ self.res_coord[i][0] for i in res_taken ], dtype='float')
		cc_coord = np.asarray([ self.res_coord[i][2] for i in res_taken ], dtype='float')
		cnt_coord = np.asarray([ np.asarray(self.res_coord[i], dtype='float').mean(axis=0) for i in res_taken ], dtype='float')
		pep_arr = []
		for i in range(len(res_taken)):
			x, y, z = cc_coord[i]
			for j in range(len(res_taken)):
				if i != j:
					x1, y1, z1 = cn_coord[j]
					x_ans = pow(( x - x1),2)
					y_ans = pow(( y - y1),2)
					z_ans = pow(( z - z1),2)
					ans = math.sqrt(x_ans + y_ans + z_ans)
					if ans < 2.4:
						pep_arr.append(res_taken[i])
						pep_arr.append(res_taken[j])
		pep_arr = sorted(set(pep_arr))
		#ln = int(len(pep_arr)/2)
		ln = len(pep_arr)
		pep_arr = pep_arr[:ln]
		res_sparse = []
		p1 = len(pep_arr)
		for i in range(len(res_taken)):
			minim = []
			if res_taken[i] not in pep_arr:
				x, y, z = ca_coord[i]
				for j in range(len(res_taken)):
					if i != j and res_taken[j] not in pep_arr:
						x1, y1, z1 = cn_coord[j]
						x_ans = pow(( x - x1),2)
						y_ans = pow(( y - y1),2)
						z_ans = pow(( z - z1),2)
						ans = math.sqrt(x_ans + y_ans + z_ans)
						minim.append(ans)
						if ans < 4.0:
							#print 'cnt ca ', ans
							res_sparse.append(res_taken[i])
							res_sparse.append(res_taken[j])
		
		for i in range(len(res_taken)):
			minim = []
			if res_taken[i] not in pep_arr:
				x, y, z = cn_coord[i]
				for j in range(len(res_taken)):
					if i != j and res_taken[j] not in pep_arr:
						#print res_taken[i], res_taken[j], '++++'
						x1, y1, z1 = cn_coord[j]
						x_ans = pow(( x - x1),2)
						y_ans = pow(( y - y1),2)
						z_ans = pow(( z - z1),2)
						ans = math.sqrt(x_ans + y_ans + z_ans)
						minim.append(ans)
						if ans < 4.0: # keep 4
							#print 'cnt cn', ans
							res_sparse.append(res_taken[i])
							res_sparse.append(res_taken[j])
		
		res_sparse = sorted(set(res_sparse))	
		#res_sparse = []
		res_sparse = res_sparse[1:]
		res_sparse = [ i for i in res_sparse if i not in pep_arr ]
		ln = len(res_taken)

		
		
		if not res_sparse:
			
			#print 'if'
			res_taken = self.FindNewRes(res_taken)
			return res_taken
		else:
			#print 'else-----'
			res_taken = [ i for i in res_taken if i not in res_sparse ]
			res_taken = self.FindNewRes(res_taken)
			return res_taken
		
			
		#res_taken = self.FindNewRes(res_taken)
		return res_taken
			
			
	
		
			
		
	def GiveStretch(self, res_stretch, lig, res_taken):
		#res_stretch = 'VAL A1409/-8_3_-3/6.4+GLY A3146/-5_0_-3/12.1+ASP A2651/-8_7_-7/7.5+ASP A2650/-2_0_-6/5.9+ALA A1214/-5_2_-8/5.3'
		res_taken = copy.deepcopy(res_taken)
		random.shuffle(self.res_id)
		x1, y1, z1 = self.het_coord[lig]
		res_total_grid = []
		
		
		ResCurrent = [ i.split('/')[0] for i in res_stretch.split('+') ]
		
		for i in ResCurrent:
			ClashCheck, ClashPairs = self.check_clash(res_taken, i)
			if not ClashCheck:
				res_taken.append(i)	
			
		#self.ClashOrderResidue
		
		#print(res_taken)
		L = len(res_taken)
		res_taken = self.NoMapRes1(res_taken)
		
		if len(res_taken) < L:
			res_taken = self.NoMapRes1(res_taken)
			
		return res_taken
								
	
	
	
	
	
	def PickResStretch(self):	
		lig_taken = random.choice(self.lig_atom_id)
		#lig_taken = '8192' # '8211'
		#print(lig_taken)
		#print(random.choice(self.res_lig_points[lig_taken]))
		
		L = len(self.res_lig_points[lig_taken]) // 10
		if L <= 10:
			L = 10
		arr = []
		for _ in range(L):
			arr.append(random.choice(self.res_lig_points[lig_taken]))
		arr = sorted(arr, key = lambda x:x[1], reverse=True)[:50]
		
		#for i in arr[:5]:
		#	print i
		arr = arr[:10]
		random.shuffle(arr)
		
		'''
		arr2 = []
		for i in arr:
			for j in i[0][3].split(' '):
				val = self.associated_grid_residues[j]
				val = sorted(val, key = lambda x:int(x[1]), reverse=True)
				print val[:10]
				print val[-10:]
				time.sleep(11)
		'''
		
		res = arr[0]
		#print res
		r1 = res[0][1].split('-')
		r2 = res[0][3].split(' ')
		r3 = res[0][4].split(' ')	
		arr = []
		for j in range(len(r1)):
			arr.append(r1[j]+'/'+r2[j]+'/'+r3[j])
		#print '+'.join(arr), lig_taken,'got'
		#time.sleep(11)
		return  '+'.join(arr), lig_taken
		
		
		
	
	def file_read(self, nos, ligand_aline, minim, maxim):
		global NoMapReslen 
		for no in range(nos):
			prev, cur = 0, 0
			count = 0
			stop_count = 0
			res_taken = []
			check = True
			lig_track_dic = {}
			rng = random.randint(minim, maxim)
			#rng = 30
			#res_taken = ['GLU A1501', 'TYR A 112', 'CYS A  33', 'PHE A 446', 'ALA A 419', 'MET A 165', 'GLY A2976', 'PRO A  60', 'GLY A2629', 'GLY A2977', 'HIS A 704', 'ALA A 538', 'LYS A 363', 'GLY A3470', 'MET A 724', 'GLN A 612', 'PHE A1307', 'LEU A 129', 'GLY A  55', 'ASP A2806', 'THR A  99', 'GLY A3471', 'LEU A 679', 'GLY A2630']
			while check:
				count += 1
				if not res_taken:
					res_stretch, lig = self.PickResStretch()
					lig_track_dic[lig] = 0
					res_taken = self.GiveStretch(res_stretch, lig, res_taken)
					#print res_taken
				else:
					#print 'else'
					res_stretch, lig = self.PickResStretch()
					if not res_stretch:
						count -= 1
						continue
					#print res_stretch	
					res_taken = self.GiveStretch(res_stretch, lig, res_taken)
					if len(res_taken) >= rng:
						res_taken = res_taken[:rng]
						check = False
			
				#print res_taken, len(res_taken), rng
				prev = len(res_taken)
				
				if abs(cur - prev) <= 0:
					#print cur, prev, cur-prev, stop_count
					stop_count += 1
					if stop_count > 10: # > 7
						#print('STOP')
						check = False
				else:
					stop_count = 0
				cur = prev		
				#time.sleep(1)
				#print '\n'
			return res_taken
			
'''			
			
Folder = 'SAMSingleOut'

res_point_line = open(dire+'/'+Folder+'/residue_position', 'r').readlines()	
ligand_aline = open(dire+'/'+Folder+'/lig.pdb', 'r').readlines()	
site_gen = SiteGen(res_point_line, ligand_aline, Folder)
		


dic = defaultdict(list)
dire = os.getcwd()
aline = open(dire+'/'+Folder+'/final2.pdb', 'r').readlines()
for line in aline:
	line = line.strip()
	dic[line[17:26]].append(line)

nos = 50
for i in range(nos):	
	t = time.time()
	sites = site_gen.file_read(1, ligand_aline, 20, 26)	
	#sites = site_gen.file_read(1, ligand_aline, 10, 15)	
	out = open(dire+'/'+Folder+'/temp2a/'+str(i)+'a.pdb', 'w')
	
	print (i)
	for j in sites:
		#print j
		for k in dic[j]:
			out.write(k+'\n')
	s = time.time()
	print (s-t,'\n')		
	out.close()	

'''


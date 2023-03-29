import os
import time
from collections import Counter, defaultdict
import math
import numpy as np
from numba import jit
import random
from sklearn.preprocessing import MinMaxScaler
import numba

from CustomClasses.RankSite import RankSite

dire = os.getcwd()


@jit(nopython=True, cache=True)
def Calc(MatrixVal, TotalMatrix):	
	#return 1
	L = len(MatrixVal)
	TotalCount = 0		
			
	ln = len(TotalMatrix)			
	
	for i in TotalMatrix:
		#print(i.shape)
		#print(MatrixVal.shape)
		
		count = 0
		count1 = 0
		
		CM = np.zeros(20, dtype='float')#[0]*20 # Count Matrix
		for ligNo in range(L):
			for j in range(20):
				A1 = i[ligNo][j][np.where(i[ligNo][j] > 0)]
				A2 = MatrixVal[ligNo][j][np.where(MatrixVal[ligNo][j] > 0)]	
				#print(A1, A2)			
				local = sum([ 1 for k in A1 if k in A2 ])
				#local = (local * LM[ligNo])/ln
				#local = local * ResWt[j]
				CM[j] = CM[j]+local
		

		TotalCount += CM.sum()
		
	TotalCount = float(TotalCount)
		
	return TotalCount//10




class Fitness():
	
	def __init__(self, Folder):
		self.res_coord = defaultdict(list)
		aline = open(dire+'/'+Folder+'/final2.pdb', 'r').readlines()
		for line in aline:
			line = line.strip()
			if line[:4] == 'ATOM':
				#self.res_coord[line[17:26]].append([line[28:38].strip(), line[38:46].strip(), line[46:54].strip()])	
				self.res_coord[line[17:26]].append(line)	

		self.KMeanClusterRes = defaultdict(list)
		aline = open(dire+'/'+Folder+'/clusters.txt', 'r').readlines()
		for line in aline:
			line = line.strip()
			if len(line.split('\t')) == 3:
				#print (line)
				linea, val, val1 = line.strip().split('\t')
				val = int(val)
				val1 = int(val1)
				for i in linea.split('-'):
					self.KMeanClusterRes[i].append(val) # earlier val1 was used
		
		self.stride_line = open(dire+'/'+Folder+'/StrideMap2.txt', 'r').readlines()
		
		aline = open(dire+'/'+Folder+'/StrideMap1.txt', 'r').readlines()
		self.ResidueGeometry = {}
		for line in aline:
			line = line.strip()
			l = line.split('\t')
			self.ResidueGeometry[l[0]] = l[1:]
		
	
	def GeometryFn(self, aline):
		return 1
		Residue = {}
		for line in aline:
			line = line.strip()
			if line[:4] == 'ATOM' and line[13:16].strip() == 'CA':
				val = [ line[28:38].strip(), line[38:46].strip(), line[46:54].strip()  ]
				val = list(map(float, val))
				val = list(map(int, val))
				val = list(map(str, val))
				val = '_'.join(val)
				if line[17:26] in self.ResidueGeometry:
					Residue[line[17:20]+' '+val] = line[17:26]
				#Residue.append(line[17:20]+' '+val)
		
		
		Score = []
		for i in Residue.keys():
			#print(i)
			for line in self.stride_line:
				line = line.strip()
				l = line.split('\t')
				count = 0
				if len(l) > 2:
					#print(line)
					try:
						coord = l[1].split('/')
						res = l[2].split('/')
						ss = l[3].split('/')
						phi = list(map(float,l[4].split('/')))
						psi = list(map(float,l[5].split('/')))
						#print(l)
						for j in range(len(res)):
							if res[j]+' '+coord[j] in Residue:
								if ss[j] == self.ResidueGeometry[Residue[res[j]+' '+coord[j]]][0]:
									angle_phi = 180 - abs(abs(phi[j] - float(self.ResidueGeometry[Residue[res[j]+' '+coord[j]]][1])) - 180)
									angle_psi = 180 - abs(abs(psi[j] - float(self.ResidueGeometry[Residue[res[j]+' '+coord[j]]][2])) - 180)
									if angle_phi < 10 and angle_psi < 10:
										#print(ss[j], phi[j], psi[j], ResidueGeometry[Residue[res[j]+' '+coord[j]]], angle_phi, angle_psi)
										count += 1
					except:
						pass					
				Score.append(count)
		Score = np.asarray(Score, dtype='int')		
		val = sum([ i**i for i in Score ])
		#print(val, max(Score))
		return val


	
	'''		
		
	def GroupResidues(self, aline, lig_aline):
		ResArr = ['GLY', 'ALA', 'VAL', 'LEU', 'ILE', 'THR', 'SER', 'TYR', 'TRP',
			'PRO', 'PHE', 'ASN', 'GLN', 'ASP', 'GLU', 'ARG', 'LYS', 'HIS', 'CYS', 'MET']
		het_arr = []
		dic = defaultdict(list)
		for i in aline:
			i1 = i.strip()
			if i1[:4] == 'ATOM' and i1[13:16].strip() == 'CA':
				dic[i1[17:26]].append([ i1[28:38].strip(), i1[38:46].strip(), i1[46:54].strip() ])
		for i in lig_aline:
			i1 = i.strip()		
			if i1[:6] == 'HETATM':
				het_arr.append([ i1[6:11].strip(), i1[28:38].strip(), i1[38:46].strip(), i1[46:54].strip() ])
		
		prot_arr = []
		for i,j in dic.items():
			prot_arr.append([i,np.asarray(j, dtype='float').mean(axis=0)])
			
		Matrix = []
		for i in het_arr:
			#print i
			het_res_arr = []
			#print i, query_count[lig_dic[i[0]]]
			x, y, z = list(map(float,i[1:]))
			het_dic = {}
			for i1 in ResArr:
				het_dic[i1] = 0
			
			for j in prot_arr:
				x1, y1, z1 = j[1]
				#print j
				x_ans = pow(( x - x1),2)
				y_ans = pow(( y - y1),2)
				z_ans = pow(( z - z1),2)
				ans = math.sqrt(x_ans + y_ans + z_ans)
				#print ans
				het_res_arr.append([i[0], j[0][:3], ans])
			het_res_arr = sorted(het_res_arr, key = lambda x:float(x[2]))	
		
			dic_temp = {}	
			for i in het_res_arr:
				if i[1] not in dic_temp:
					dic_temp[i[1]] = 0
					het_dic[i[1]] = float(i[2])
					
			mat = []		
			for j in ResArr:
				mat.append(het_dic[j])
			#print mat	
			Matrix.append(mat)
		Matrix = np.asarray(Matrix, dtype='float')
		Matrix = [ np.round(i, decimals=1) for i in Matrix ]
		Matrix = np.asarray(Matrix, dtype='float')
		return Matrix
	'''
	
	def GroupResidues(self, arr, lig_aline):
		ResArr = ['GLY', 'ALA', 'VAL', 'LEU', 'ILE', 'THR', 'SER', 'TYR', 'TRP',
			'PRO', 'PHE', 'ASN', 'GLN', 'ASP', 'GLU', 'ARG', 'LYS', 'HIS', 'CYS', 'MET']
			
				
		het_arr = []
		dic = defaultdict(list)
		for i in arr:
			i1 = i.strip()
			if i1[:4] == 'ATOM' and i1[13:16].strip() == 'CA':
				dic[i1[17:26]].append([ i1[28:38].strip(), i1[38:46].strip(), i1[46:54].strip() ])
		for i in lig_aline:
			i1 = i.strip()		
			if i1[:6] == 'HETATM':
				het_arr.append([ i1[6:11].strip(), i1[28:38].strip(), i1[38:46].strip(), i1[46:54].strip() ])

		prot_arr = []
		for i,j in dic.items():
			prot_arr.append([i,np.asarray(j, dtype='float').mean(axis=0)])
			
		Matrix = []
		for i in het_arr:
			x, y, z = list(map(float,i[1:]))
			
			mat = []
			
			for i1 in range(len(ResArr)):
				Ans = [0]*10
				index = 0
				for j in prot_arr:
					x1, y1, z1 = j[1]
					if j[0][:3] == ResArr[i1]:
						x_ans = pow(( x - x1),2)
						y_ans = pow(( y - y1),2)
						z_ans = pow(( z - z1),2)
						ans = math.sqrt(x_ans + y_ans + z_ans)
						ans = int(ans)
						if index < 9:
							Ans[index] = ans
						index += 1	
				
				mat.append(Ans)
			
			Matrix.append(mat)
		Matrix = np.asarray(Matrix, dtype='int')
		return Matrix
		
	
	


	'''	
	def Calc(self, MatrixVal, TotalMatrix):
		c = 0
		ln = []
		for k in range(1):
			for s in range(len(TotalMatrix)):
				c1, c2 = 0, 0
				for i in range(len(TotalMatrix[s])):
					c1_local = 0
					for j in range(20):
						 
						if TotalMatrix[s][i][j] != 0:
							if abs(MatrixVal[i][j] - TotalMatrix[s][i][j]) < .35:
								#c1 += 1	
								c1_local += 1
							else:
								c2 += abs(MatrixVal[i][j] - TotalMatrix[s][i][j])
								#c2 += 1
					#print c1, c2	
					c1 += c1_local**c1_local
				ln.append([c1, int(c2)])
				
		
		
		ln = sorted(ln, key=lambda x:int(x[0]), reverse=True)
		#print ln[-10:]
		ln = [ [ln[i][0], ln[i][1], i] for i in range(len(ln)) if ln[i][0] > 5 ]
		#print c
		
		#print ln[:10], len(ln)	
		#print TotalMatrix.shape
		if not ln:
			return 0.00001
			
		ln = sorted(ln, key=lambda x:int(x[0]), reverse=True)
		matched_val = float(sum([ i[0] for i in ln ]))/len(ln)
		matched_unval = float(sum([ i[1] for i in ln ]))/len(ln)
		
		#print ln[:10],matched_val, matched_unval
		#print np.percentile([ i[0] for i in ln ], 90) / np.percentile([ i[1] for i in ln ], 90)
		
		#print np.percentile([ i[0] for i in ln ], 90)
		##print np.percentile([ i[1] for i in ln ], 90)
		#print np.percentile([ i[0] for i in ln ], 90) / np.percentile([ i[1] for i in ln ], 90) 
		#print '\n'
		try:
			if np.percentile([ i[1] for i in ln ], 90) == 0:
				return 0.00001
			return np.percentile([ i[0] for i in ln ], 90) / np.percentile([ i[1] for i in ln ], 90) # this is best
		except:
			return 0.00001

	'''


	def SortScore(self, SiteArr, arr, rev):	
		ranks = []
		
		if rev:
			#print(arr, 'from Sort Score Function')
			arr1 = np.array(arr, dtype='float')
			arr1 = np.reshape(arr1, (-1,1))
			#print arr1.shape
			scaler = MinMaxScaler()
			normalized_X = scaler.fit_transform(arr1)
			for i in range(len(normalized_X)):
				ranks.append([ SiteArr[i] , 100-int(normalized_X[i][0]*100) ])
		else:
			arr1 = np.array(arr, dtype='float')
			arr1 = np.reshape(arr1, (-1,1))
			scaler = MinMaxScaler()
			normalized_X = scaler.fit_transform(arr1)	
			for i in range(len(normalized_X)):
				ranks.append([ SiteArr[i] , int(normalized_X[i][0]*100) ])
		#print ranks
		dic1 = {}
		for i in ranks:
			dic1[i[0]] = i[1]
		return dic1	
		
			
		
	def file_read(self, aline, Folder):
		lig_aline = open(dire+'/'+Folder+'/lig.pdb', 'r').readlines()
		TotalMatrix = np.load(dire+'/'+Folder+'/Association2.npy')
		
						 
		res_count_stretch = {'GLY':1 ,'ALA':1 ,'LEU':1 ,'ILE':1 ,'TRP':1 ,'SER':1 ,'THR':1 ,'TYR':1 ,'PHE':1,\
						 'PRO':1 ,'ASP':1,'GLU':1 ,'HIS':1 ,'CYS':1 ,'MET':1 ,'VAL':1 ,'ASN':1 ,\
						 'GLN':1 ,'ARG':1 ,'LYS':1 }
						 
		SiteArr, AssociationRank, GridRank, EnergyRank, GeometryRank = [], [], [], [], []
		count = 0				 
		
		rank = RankSite(Folder)
		arr = rank.run(aline)
		#print(arr)
		dic = {}
		for i in arr:
			dic[int(i[0].split('.')[0])] = i[1]		
		RankArr = []
		
		for line in aline:
			l = line
			#print (l)
			i1 = l[0].split('_')
			arr = []
			for j in i1:
				#print j
				for j1 in self.res_coord[j]:
					arr.append(j1)
			#print i1
			GridVal = sum([ max(np.asarray(self.KMeanClusterRes[j])) * res_count_stretch[j[:3]] for j in i1 if j in self.KMeanClusterRes ])
			GridVal = sum([ sum(self.KMeanClusterRes[j]) for j in i1 if j in self.KMeanClusterRes ]) # show better improvements
			MatrixVal = self.GroupResidues(arr, lig_aline)
			AssVal = Calc(MatrixVal, TotalMatrix)
			GeometryRank.append(self.GeometryFn(arr))
			GridRank.append(GridVal)
			EnergyRank.append(l[1])
			AssociationRank.append(AssVal)
			SiteArr.append(count)
			RankArr.append(dic[count])
			
			count += 1
			
			#time.sleep(11)
			
		GridRank = np.asarray(GridRank, dtype='int')	
		AssociationRank = np.asarray(AssociationRank, dtype='float')	
		EnergyRank = np.asarray(EnergyRank, dtype='float')		
		GeometryRank = np.asarray(GeometryRank, dtype='int')
		RankArr = np.asarray(RankArr, dtype='float')
		
		
		#print (sorted(AssociationRank)[1] - sorted(AssociationRank)[0]) / 2, np.mean(AssociationRank)
		#print (sorted(GridRank)[1] - sorted(GridRank)[0]) / 2, np.mean(GridRank)
		
		
		#print AssociationRank
		
		GridDic = self.SortScore(SiteArr, GridRank, True)	
		EnergyDic = self.SortScore(SiteArr, EnergyRank, False)	
		AssociationDic = self.SortScore(SiteArr, AssociationRank, True)
		GeometryDic = self.SortScore(SiteArr, GeometryRank, True)
		RankDic = self.SortScore(SiteArr, RankArr, False)
		
		
		#print(GeometryDic)
		#print 'done'
		#time.sleep(11)
		Total = []
		SiteInt = defaultdict(list)
		c = 0
		for i in SiteArr:
			SiteInt[c].append(AssociationRank[i])
			SiteInt[c].append(GridRank[i])
			SiteInt[c].append(EnergyRank[i])
			#Total.append([i, GridDic[i]+AssociationDic[i]+EnergyDic[i]+EnergyDic[i]])
			#Total.append([i, GridDic[i]+EnergyDic[i]+EnergyDic[i]])
			Total.append([i, GridDic[i]+AssociationDic[i]+EnergyDic[i]+GeometryDic[i]])
			c += 1
		
		Final = []
		for i in range(len(Total)):
			s1 = Total[i][0]
			#print Total[i]
			#ranks = str(GridDic[i])+' '+str(AssociationDic[i])+' '+str(EnergyDic[i])+' '+str(EnergyDic[i])
			#print GridDic[i], EnergyDic[i], AssociationDic[i]
			ranks = float(GridDic[i]) + float(EnergyDic[i]) + float(AssociationDic[i]) + float(EnergyDic[i])
			#time.sleep(1)
			Final.append(str(s1)+'.pdb'+'\t'+str(Total[i][1])+'\t'+'_'.join(map(str,SiteInt[s1]))+'\t'+str(AssociationDic[i])+' '+str(GridDic[i])+' '+str(EnergyDic[i]))
		Final = sorted(Final, key = lambda x:float(x.split('\t')[1]))
		Final = Final[:30]
		random.shuffle(Final)
		return Final



'''

fitness = Fitness('ZITOutput')
dire = os.getcwd()

Arr = []
for i in os.listdir(dire+'/ZITOutput/site1')[:10]:
	#print i
	aline = open(dire+'/ZITOutput/site1/'+i, 'r').readlines()
	arr = sorted(set([ i[17:26] for i in aline if i[:4] == 'ATOM' ]))
	arr = '_'.join(arr)
	Arr.append([arr,-2])
	
Final = fitness.file_read(Arr,'ZITOutput')	

#print Final
'''

'''
files = os.listdir(dire+'/site1')
arr = []
for i in zip(files, Final):
	#print i
	arr.append([i[0], i[1]])
arr = sorted(arr, key = lambda x:int(x[1].split('\t')[1]), reverse=False)
for i in arr[:20]:
	print i
print '\n'

'''










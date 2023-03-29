import time
import random
import math
import numpy as np
import sys
from collections import Counter, defaultdict
import copy
import re
from numba import jit
from numba import cuda
import numba
import os
from openbabel import openbabel
import pickle
import shutil
import subprocess
from tqdm import tqdm
import random
from multiprocessing import Pool

openbabel.obErrorLog.StopLogging()
obConversion = openbabel.OBConversion()
obConversion.SetInAndOutFormats("pdb", "sdf")



dire = os.getcwd()


if len(sys.argv) == 4:
	LigFile = sys.argv[1]
	OutFold = sys.argv[2]
	nprocs = int(sys.argv[3])
else:
	print ('MCSSNumba.py <Ligand-1> <OutFolder> <Nprocs>') # site1.pdb site2.pdb
	sys.exit()


AtomFile = open(dire+'/CustomDatas/AtomIndex', 'r').readlines()
ResidueGrouping = {}
for i in AtomFile:
	i1 = i.strip().split(' ')
	ResidueGrouping[i1[0]] = int(i1[1])

aline = open(dire+'/CustomDatas/CheckMolOutput', 'r').readlines()
CheckMolPDBHash = {}
for line in aline:
	line = line.strip()
	l = line.split('\t')
	if len(l) == 1:
		CheckMolPDBHash[l[0]] = None
	else:
		CheckMolPDBHash[l[0]] = l[1]




def RunCheckmol():
	os.chdir(dire+'/'+OutFold)
	shutil.copy(dire+'/checkmol', dire+'/'+OutFold)
	rnd = 'a'
	aline = open('a.pdb','r').readlines()
	het = [ line[6:11].strip() for line in aline if line[:6] == "HETATM"]
	dic = {}
	for i in range(len(het)):
		dic[str(i+1)] = het[i]
	
	mol = openbabel.OBMol()
	obConversion.ReadFile(mol, rnd+'.pdb')
	obConversion.WriteFile(mol, rnd+'.sdf')
	
	c1 = subprocess.Popen("./checkmol -p "+rnd+".sdf", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	c2 = c1.communicate()[0]
	c2 = c2.decode('UTF-8')
	
	arr = []
	for i in c2.split('\n'):
		#zprint i
		if len(i) > 1:
			if i[0] == '#':
				check = True
				#print i
				for j in i.split(':')[2].split(','):
					j1 = j.strip('-')
					if len(j1.split('-')) > 3:
						#print j1.split('-')
						arr.append(j1)
		
	# print arr # print end
	dic1 = {}
	#print dic
	for i in arr:
		#print i
		i1 = '-'.join([dic[k] for k in i.split('-')])
		#print i1
		for j in i.split('-'):
			dic1[dic[j]] = i1
			
	arr1 = []
	for i in arr:
		j = [dic[k] for k in i.split('-')]
		arr1.append(j)
	check = True
	while True:
		check1 = False
		if not check:
			break
		arr2 = arr1
		i_c, j_c = None, None
		for i in range(0, len(arr2)):
			if check1:
				break
			for j in range(0, len(arr2)):
				if check1:
					break
				if i != j:
					summed = sum([ 1 for i1 in arr2[i] if i1 in arr2[j] ])
					if summed >= 1:
						i_c, j_c = i, j
						check1 = True
					#time.sleep(1)
		if not check1:
			break			
		arr3 = []
		for i in range(0, len(arr2)):
			if i != i_c and i != j_c:
				#print i
				arr3.append(arr2[i])
		arr3.append(sorted(set(arr2[i_c]+arr2[j_c])))	
		arr1 = arr3
				
		
	check = True
	if not arr:
		check = False
	arr1 = [ '-'.join(i) for i in arr1 ]	
	arr1 = ','.join(arr1)
	os.chdir(dire)
	return arr1
	

def CenterResidues(arr):
	arr = [ i for i in arr if i[74:78].strip() != 'H' ]
	coord = []
	for i in arr:
		coord.append([i[28:38].strip(), i[38:46].strip(), i[46:54].strip()])
	coord = np.asarray(coord, dtype='float')	
	coord -= np.mean(coord, axis=0)
	dic = {1:"   ", 2:"  ", 3:" ", 4:""}
	brr = []
	dic1 = {}
	for i in range(len(arr)):
		
		#if arr[i][13:16]+' '+arr[i][17:26] not in dic1:
		#	dic1[arr[i][13:16]+' '+arr[i][17:26]] = 0
		j = [ "%.3f"%j for j in coord[i]]
		val = ''.join([ dic[len(k.split(".")[0])]+k for k in j ])
		brr.append(arr[i][:30]+val+arr[i][54:])
	return brr

def ConstructMatrix(Arr1):

	Arr1 = [ line.strip() for line in Arr1 if line[:6] == 'HETATM' ]
	Arr1 = [ line for line in Arr1 if line[71:].strip() in ResidueGrouping ]
	Arr1 = CenterResidues(Arr1)
	
	Coord = []
	ResList = []
	for i in Arr1:
		Coord.append([i[28:38].strip(), i[38:46].strip(), i[46:54].strip()])
		ResList.append(i[6:11].strip())
	Coord = np.asarray(Coord, dtype='float')		

	L = len(Coord)
	DistMatrix = np.zeros((L, L))		
	PCMatrix = np.zeros((L))		
	
	for i in range(L):
		for j in range(L):
			#print(Coord[i],Coord[j])
			DistMatrix[i][j] = math.sqrt(pow(( Coord[i][0] - Coord[j][0] ),2) + pow(( Coord[i][1] - Coord[j][1] ),2) + pow(( Coord[i][2] - Coord[j][2] ),2))
	
	
	for i in range(len(Arr1)):
		PCMatrix[i] = ResidueGrouping[Arr1[i][71:].strip()]
	
	FullCoord = defaultdict(list)
	FullPDB = defaultdict(list)
	CoordKabsch = []
	AtomKabsch = []
	for i in Arr1:
		#print i
		FullCoord[i[17:26]].append([i[28:38].strip(), i[38:46].strip(), i[46:54].strip()])
		FullPDB[i[17:26]].append(i.strip())
		AtomKabsch.append(i)	
		#if i[13:16].strip() == 'CA':
		CoordKabsch.append([i[28:38].strip(), i[38:46].strip(), i[46:54].strip()])
		
	
	CoordKabsch = np.asarray(CoordKabsch, dtype='float')
	
	return DistMatrix, Arr1, L, PCMatrix, FullCoord, FullPDB, ResList, CoordKabsch, AtomKabsch	





T = 0.35 # Tau
S = 20





@jit(nopython=True)
def kabsch(P, Q):
	#print(P, Q)
	P_mean = np.zeros(3)
	for i in range(P.shape[1]):
		P_mean[i] = P[:,i].mean()
	P -= P_mean
	
	Q_mean = np.zeros(3)
	for i in range(Q.shape[1]):
		Q_mean[i] = Q[:,i].mean()
	Q -= Q_mean
	
	C = np.dot(np.transpose(P), Q)
	V, S, W = np.linalg.svd(C)
	d = (np.linalg.det(V) * np.linalg.det(W)) < 0.0
	if d:
		S[-1] = -S[-1]
		V[:, -1] = -V[:, -1]
	U = np.dot(V, W)
	return U
	




@jit(nopython=True)
def GetAlignments(GlobalSimilarity1, GlobalSimilarity2, CoordKabsch1, CoordKabsch2):
	count_align = -1000
	AlignResTotal = np.zeros((len(GlobalSimilarity1), 1000,2), dtype='int')-1
	GlobalCount = 0
	for i in range(len(GlobalSimilarity1)):
	#for i in range(1):
		arr1 = GlobalSimilarity1[i][np.where(GlobalSimilarity1[i] > -1)]
		arr2 = GlobalSimilarity2[i][np.where(GlobalSimilarity2[i] > -1)]

		arr1 = [ int(j) for j in arr1]
		arr2 = [ int(j) for j in arr2]
		
		if len(set(arr1)) <= 3:
			continue
		if len(arr1) != len(arr2):
			continue
		
		arr1_coord = np.zeros((len(arr1), 3))
		arr1_coord_kabsch = np.zeros((len(arr1), 3))
		#print(arr1, arr2)
		for j in range(len(arr1)):
			arr1_coord[j] = CoordKabsch1[arr1[j]]
			arr1_coord_kabsch[j] = CoordKabsch1[arr1[j]]
			
		arr2_coord = np.zeros((len(arr2), 3))
		arr2_coord_kabsch = np.zeros((len(arr2), 3))
		for j in range(len(arr2)):
			arr2_coord[j] = CoordKabsch2[arr2[j]]
			arr2_coord_kabsch[j] = CoordKabsch2[arr2[j]]
		
		
		U = kabsch(arr1_coord_kabsch, arr2_coord_kabsch)
		
		B1 = np.zeros((len(CoordKabsch1),3))
		for j in range(len(CoordKabsch1)):
			B1[j] = CoordKabsch1[j]
		
		
		c1_res = np.zeros(3)
		for j in range(3):
			c1_res[j] = arr1_coord[:,j].mean()
		c2_res = np.zeros(3)
		for j in range(3):
			c2_res[j] = arr2_coord[:,j].mean()	
		
		#print(B1, c1_res, arr1_coord)
		#print(c1_res)
		#print(arr1_coord.mean(axis=0))
		B1 -= c1_res
		B1 = np.dot(B1, U)
		B1 += c2_res
		AlignRes = np.zeros((1000,2), dtype='int')-1
		#print (B1[:3], CoordKabsch2[:3])
		count = 0
		for j in range(len(B1)):
			for k in range(len(CoordKabsch2)):
				ans = math.sqrt(math.pow(( B1[j][0] - CoordKabsch2[k][0] ),2) + math.pow(( B1[j][1] - CoordKabsch2[k][1] ),2) + math.pow(( B1[j][2] - CoordKabsch2[k][2] ),2))
				if ans < .45:
					AlignRes[count] = [j,k]
					count += 1
		#print(count)
		#print(AlignRes[:10])
		if count > 3:
			AlignResTotal[GlobalCount] = AlignRes
			GlobalCount += 1	
	return 100, AlignResTotal[:GlobalCount]				
		


@jit(nopython=True)
def CheckExtension(S1, S2, P1, P2, DistMatrix1, DistMatrix2):
	
	for i in range(len(S1)):
		if S1[i] != -1:
			#print (S1, S2, P1, P2)
			#print('\n\n')
			
			if abs( DistMatrix1[S1[i]][int(P1)] - DistMatrix2[S2[i]][int(P2)] ) > T:
				return False
	return True



@jit(nopython=True)
def GetCommanElements(DistBin1a, DistBin1b, DistBin1c, DistBin2a, DistBin2b, DistBin2c, DistMatrix1, L1, PCMatrix1, DistMatrix2, L2, PCMatrix2, CoordKabsch1, CoordKabsch2):
	# Residue index, distance of nearby residues, index of nearby residues
	
	#print(DistBin2c)
	count = 0
	Pairs1a = np.zeros(100)
	Pairs1b = np.zeros(100)
	Pairs1a = Pairs1a-1
	Pairs1b = Pairs1b-1
	
	Pairs2a = np.zeros(100)
	Pairs2b = np.zeros(100)
	Pairs2a = Pairs2a-1
	Pairs2b = Pairs2b-1
	
	for i in range(1,len(DistBin1b)):
		for j in range(1,len(DistBin2b)):	
			#print (DistBin1c[i], DistBin2c[j], j, len(DistBin2b))
			if PCMatrix1[DistBin1c[i]] == PCMatrix2[DistBin2c[j]]:
				if DistBin1c[i] > -1 and DistBin2c[j] > -1: # Tag Changed
					if abs(DistBin1b[i] - DistBin2b[j] ) < T:
						#print (i, j,len(DistBin2b),DistBin1a) 
						#time.sleep(1)
						if count >=100:
							break
						Pairs1a[count] = DistBin1a
						Pairs1b[count] = DistBin1c[i]
						Pairs2a[count] = DistBin2a
						Pairs2b[count] = DistBin2c[j] # Tag Changed
						count += 1
	
	
	Pairs1a = Pairs1a[:count]		# res1, res2
	Pairs1b = Pairs1b[:count]	
	Pairs2a = Pairs2a[:count]	
	Pairs2b = Pairs2b[:count]			
	
	#print(Pairs1a)
	
	Pairs = np.zeros((len(Pairs1a),5))
	for i in range(len(Pairs1a)):
		Pairs[i] = [ Pairs1a[i], Pairs1b[i], Pairs2a[i], Pairs2b[i], abs(Pairs1b[i] - Pairs2b[i])   ] # this takes care of order from each residues
	#Pairs = sorted(Pairs, key = lambda x:int(x[4])) # this is not numpy array
	#Pairs = np.array(Pairs, dtype='float')
	ind = np.argsort(Pairs[:,4])
	#print(np.take_along_axis(Pairs, ind, axis=0)) 
	Pairs = Pairs[ind]
	
	
	PairsHash = np.zeros((len(Pairs),5))
	PairsHashCount = -1
	PairsHashIndexCount = np.zeros(1000, dtype='int')
	PairsHashIndexCount -= 1
	
	for i in range(len(Pairs)):
		
		c = 0
		for j in range(len(PairsHash)):
			
			check = True
			for k in range(len(Pairs[i])):
				if Pairs[i][k] != PairsHash[j][k]:
					check = False
			if check:
				c = 1
		if c == 0:
			PairsHashCount += 1
			PairsHash[PairsHashCount] = Pairs[i]
			PairsHashIndexCount[PairsHashCount] = i		
		
	PairsHashIndexCount = PairsHashIndexCount[np.where(PairsHashIndexCount > -1)]	
	Pairs = Pairs[PairsHashIndexCount]	
	
	
	Pairs1a = np.zeros((len(Pairs),2))
	Pairs2a = np.zeros((len(Pairs),2))
	
	for i in range(len(Pairs)):
		Pairs1a[i] = [ Pairs[i][0], Pairs[i][1] ] 
		Pairs2a[i] = [ Pairs[i][2], Pairs[i][3] ] 
	
	if len(Pairs1a) == 0:
		return -10, np.zeros((1,1000,2), dtype='int')-1
	Answer = np.zeros(len(Pairs1a))
	
	GlobalSimilarity1 = np.zeros((len(Pairs1a),100))
	GlobalSimilarity2 = np.zeros((len(Pairs1a),100))
	Gcount = 0
	
	CutOff = 3
	CutOffCount = 0
	for i in range(len(Pairs1a)):
		#i = 36
		#print(Pairs1a)
		
		Visited1 = np.zeros(30, dtype='int')
		Visited1 = Visited1-1
		Visited1[0] = Pairs1a[i][1]
		
		Visited2 = np.zeros(30, dtype='int')
		Visited2 = Visited2-1
		Visited2[0] = Pairs2a[i][1]
		
		
		Similarity1 = np.zeros(100, dtype='int')
		Similarity2 = np.zeros(100, dtype='int')
		Similarity1 = Similarity1-1
		Similarity2 = Similarity2-1
		Similarity1[0] = Pairs1a[i][0]
		Similarity1[1] = Pairs1a[i][1]
		Similarity2[0] = Pairs2a[i][0]
		Similarity2[1] = Pairs2a[i][1]
		
		#print(Pairs1a[i], DistMatrix1[int(Pairs1a[i][0])][int(Pairs1a[i][1])])
		
		if DistMatrix1[int(Pairs1a[i][0])][int(Pairs1a[i][1])] > 2.0:
			continue
		
		CheckCount = 2
		HashCount = 1
		for j in range(len(Pairs2a)):
			#print (Pairs1a[j][1])
			#j = 38
			#if int(Pairs1a[j][1]) not in Visited1 and int(Pairs2a[j][1]) not in Visited2:
				
			if PCMatrix1[int(Pairs1a[j][1])] == PCMatrix2[int(Pairs2a[j][1])]:
				#print('got', int(Pairs1a[j][1]))
				if CheckExtension(Similarity1, Similarity2, Pairs1a[j][1], Pairs2a[j][1], DistMatrix1, DistMatrix2):
					if HashCount > 25:
						break
					if int(Pairs1a[j][1]) not in Visited1:
						if min([ DistMatrix1[d1][int(Pairs1a[j][1])] for d1 in Similarity1[:CheckCount] ]) < 2:
							#print (Similarity1, Similarity1[:CheckCount], Pairs1a[j][1])
							#time.sleep(1)
							Similarity1[CheckCount] = Pairs1a[j][1]
							Similarity2[CheckCount] = Pairs2a[j][1]
							Visited1[HashCount] = int(Pairs1a[j][1])
							Visited2[HashCount] = int(Pairs2a[j][1])
							#print (Visited1, Visited2)
							CheckCount += 1
							HashCount += 1

		Answer[i] = len(np.unique(Similarity1))-1
		if len(np.where(Similarity1 > -1)[0]) > CutOff: # revisit. change to 3
			GlobalSimilarity1[Gcount] = Similarity1
			GlobalSimilarity2[Gcount] = Similarity2
			#print(Similarity1, np.where(Similarity1 > -1)[0], len(np.where(Similarity1 > -1)) )
			Gcount += 1
			CutOffCount += 1
			if CutOffCount > -1:
				CutOff = 3
			
	#print (Answer)
	
	
	#print(GlobalSimilarity1[:3])
	#print(Gcount)
	if Gcount == 0:
		return -1000, np.zeros((1,1000,2), dtype='int')-1
	GlobalSimilarity1 = GlobalSimilarity1[:Gcount]
	GlobalSimilarity2 = GlobalSimilarity2[:Gcount]
	#GIndex = [ i for i in range(len(GlobalSimilarity1)) if len( np.where(GlobalSimilarity1[i] > 0)[0] ) > 1   ] 
	
	#GlobalSimilarity2 = [ GlobalSimilarity2[i] for i in GIndex ]
	#GlobalSimilarity1 = [ GlobalSimilarity1[i] for i in GIndex ]
	#G1 = np.asarray(G1)
	#GlobalSimilarity1 = np.asarray(GlobalSimilarity1)
	#print('\n')
	
	
	Gs = [ sorted(i) for i in GlobalSimilarity1 ]
	#print(Gs)
	Gs = np.asarray(Gs, dtype='int')
	GsZero = np.zeros(Gs.shape, dtype='int')
	
	IndexCount = np.zeros(1000, dtype='int')
	IndexCount -= 1
	
	count = -1
	
	for i in range(len(Gs)):
		c = 1
		#print(i)
		for j in range(len(GsZero)):
			check = True
			for i_k in range(len(Gs[i])):
				if Gs[i][i_k] != GsZero[j][i_k]:
					check = False # no longer same
			if check:
				c = 0
		
		if c == 1:
			#print(i)
			count += 1
			#print(count)
			GsZero[count] = Gs[i]
			IndexCount[count] = i
	
	IndexCount = IndexCount[np.where(IndexCount > -1)]
	

	GlobalSimilarity1 = GlobalSimilarity1[IndexCount]
	GlobalSimilarity2 = GlobalSimilarity2[IndexCount]
	
	#print(GlobalSimilarity1)
	
	#return 100, np.zeros((1,1000,2))
	return GetAlignments(GlobalSimilarity1, GlobalSimilarity2, CoordKabsch1, CoordKabsch2)
	#return max(Answer)	
	#print('\n\n')
			
			

@jit(nopython=True)
def Run(DistMatrix1, L1, PCMatrix1, DistMatrix2, L2, PCMatrix2, CoordKabsch1, CoordKabsch2):
	#for n in numba.prange(10000):
	
	
	ResPos1 = np.asarray(list(range(L1)))

	DistBin1a = np.zeros(L1, dtype='int')
	DistBin1a = DistBin1a-1
	DistBin1b = np.zeros((L1,S), dtype='float')
	DistBin1b = DistBin1b-1
	DistBin1c = np.zeros((L1,S), dtype='int')
	DistBin1c = DistBin1c-1
	
	
	for i in range(L1):
		DistStacked = np.stack((DistMatrix1[i],ResPos1), axis=1)
		DistStacked = sorted(DistStacked, key = lambda x:float(x[0]))[:S]
		DistStacked = [ [j[0], int(j[1])] for j in DistStacked ]
		#print(DistStacked)
		DistBin1a[i] = i # Residue index
		#print(S, len(DistStacked))
		#print(DistBin1a,'dist 1')
		#for j in range(len(DistStacked)):
		#	DistBin1b[i][j] = DistStacked[j][0]
		#	DistBin1c[i][j] = int(DistStacked[j][1])
		
		
		DistBin1b[i][:len(DistStacked)] = [ j[0] for j in DistStacked ] # distance of nearby residues
		DistBin1c[i][:len(DistStacked)] = [ int(j[1]) for j in DistStacked ] # index of nearby residues
	
	#print(len(DistBin1a))
	
	
	ResPos2 = np.asarray(list(range(L2)))
	#print(L2,'--')
	DistBin2a = [-1]*L2
	DistBin2b = [[-1]*S]*L2
	DistBin2c = [[-1]*S]*L2
	#DistBin2a = np.asarray(DistBin2a, dtype='int')
	DistBin2b = np.asarray(DistBin2b, dtype='float')
	DistBin2c = np.asarray(DistBin2c, dtype='int')
	
	DistBin2a = np.zeros(L2, dtype='int')
	DistBin2a = DistBin2a-1
	DistBin2b = np.zeros((L2,S), dtype='float')
	DistBin2b = DistBin2b-1
	DistBin2c = np.zeros((L2,S), dtype='int')
	DistBin2c = DistBin2c-1
	
	
	#print(DistMatrix2[i])
	for i in range(L2):
		#print(i,' is i',)
		DistStacked = np.stack((DistMatrix2[i],ResPos2), axis=1)
		
		DistStacked = sorted(DistStacked, key = lambda x:float(x[0]))[:S]
		#print (DistStacked)
		#print(DistBin2a, len(DistBin2a), len(DistBin2b[i]), len(DistStacked))
		
		DistStacked = [ [j[0], int(j[1])] for j in DistStacked ]
		#print( [ j[0] for j in DistStacked ] )
		#print(i, L2, list(range(L2)))
		DistBin2a[i] = i # Residue index
		
		'''
		for j in range(len(DistStacked)):
			#print(j)
			DistBin2b[i][j] = DistStacked[j][0]
			DistBin2c[i][j] = int(DistStacked[j][1])
		'''
		
		DistBin2b[i][:len(DistStacked)] = [ j[0] for j in DistStacked ] # distance of nearby residues
		DistBin2c[i][:len(DistStacked)] = [ int(j[1]) for j in DistStacked ] # index of nearby residues
		#print('\n\n')
	FinalStore = -1000
	
	
	
	
	AlignResTotal = np.zeros((500, 1000, 2), dtype='int')-1 
	AlignCount = 0
	for i in range(len(DistBin1a)):
	#for i in range(1):
		for j in range(0,len(DistBin2a)):
			if AlignCount > 490:
				break
			if PCMatrix1[i] == PCMatrix2[j]:
				#print (i,j, DistBin1a[i])
				Error, AlignRes = GetCommanElements(DistBin1a[i], DistBin1b[i], DistBin1c[i], DistBin2a[j], DistBin2b[j], DistBin2c[j],
				DistMatrix1, L1, PCMatrix1, DistMatrix2, L2, PCMatrix2, CoordKabsch1, CoordKabsch2)
				if Error > 1:
					for k in AlignRes:
						if AlignCount > 490:
							break
						if len([ 1 for l in k if l[0] > -1 ]) > 3:
							AlignResTotal[AlignCount] = k
							AlignCount += 1
				#break		
	#print(AlignCount)			
	return AlignResTotal[:AlignCount]		



def PickSequence(Alignments, CheckmolQuery, CheckmolTarget):
	if not CheckmolQuery:
		CheckmolQuery = ""
	if not CheckmolTarget:
		CheckmolTarget = ""
		
	query_total = {}
	for i in CheckmolQuery.split(','):
		for j in i.split('-'):
			query_total[j] = 0
			
	target_total = {}
	for i in CheckmolTarget.split(','):
		for j in i.split('-'):
			target_total[j] = 0		
			
	CheckmolQueryArr = CheckmolQuery.split(',')	
	CheckmolTargetArr = CheckmolTarget.split(',')
	#print(len(Alignments))
	NewAlignments = []
	for alignment in Alignments:
		
		query_check, target_check = False, False
		query_local, target_local = {}, {}
		for i in alignment:
			query_local[i[0]] = 0
			target_local[i[1]] = 0
			if i[0] in query_total:
				query_check = True
			if i[1] in target_total:
				target_check = True
				
		if query_check and target_check:
			
			temp = []
			for k in CheckmolQueryArr:
				count = 0
				for l in k.split('-'):
					if l not in query_local:
						count = 1
				temp.append(count)
			T1 = min(temp)
			
			temp = []
			for k in CheckmolTargetArr:
				count = 0
				for l in k.split('-'):
					if l not in target_local:
						count = 1
				temp.append(count)
			T2 = min(temp)
			
			if T1 == 0 and T2 == 0:
				#print(alignment, len(alignment))
				NewAlignments.append(alignment)
		
		if not query_check and not target_check:
			NewAlignments.append(alignment)

	return len(NewAlignments)
	


def FindConnected(mcss, DistMatrix1, DistMatrix2):
	#mcss.append([10,10])
	#mcss.append([11,10])
	mcss_temp = []
	#print(mcss)
	dic = {}
	L = []
	while True:
		mcss_local = [ i for i in mcss if i[0] not in dic ]
		#print(mcss_temp, mcss_local, dic)
		if not mcss_local:
			break
		i1 = mcss_local[0][0]
		#print(i1)
		dic[i1] = 0
		
		ans = []
		if not mcss_temp:
			mcss_temp.append(i1)
			
		else:
			mcss_temp.append(i1)
			check = False
			for j in mcss:
				if j[0] not in dic:
					#print(i1, j[0], DistMatrix1[i1][j[0]])
					ans = DistMatrix1[i1][j[0]]
					if ans < 2.5:
						check = True
						#print(i1,j[0],  DistMatrix1[i1][j[0]])
						mcss_temp.append(j[0])
						dic[j[0]] = 0
						break
			if not check:
				#print(mcss_temp)
				L.append(mcss_temp)
				mcss_temp = []
				#dic = {}
	if not L:
		return [[1,1]]			
	L = sorted(L, key = lambda x:len(x), reverse=True)
	#print(L[0])
	dic = { i:0 for i in L[0] }
	#print(dic)
	mcss = [ i for i in mcss if i[0] in dic ]
	
	return mcss


path_out_fold = dire+'/'+OutFold+'/numbaTemp'
if not os.path.exists(path_out_fold):
	os.mkdir(path_out_fold)
else:
	shutil.rmtree(path_out_fold)
	os.mkdir(path_out_fold)	
	
	
	

aline = open(LigFile, 'r').readlines()

out = open(dire+'/'+OutFold+'/a.pdb','w')
het_aline = [ line for line in aline if line[:6] == 'HETATM' ] 
for line in het_aline:
	out.write(line)
out.close()	

CheckmolQuery = RunCheckmol()
#print(CheckmolQuery)

DistMatrix1t, Arr1, L1t, PCMatrix1t, FullCoord1, FullPDB1, ResList1, CoordKabsch1t, AtomKabsch1 = ConstructMatrix(aline)
CoordKabsch1t = np.asarray(CoordKabsch1t, dtype='float')

#
#Run(DistMatrix1t, L1t, PCMatrix1t, DistMatrix1t, L1t, PCMatrix1t, CoordKabsch1t, CoordKabsch1t)
#print('Running the Fast Scan module of Fragment Match...')

'''
out = open(dire+'/'+OutFold+'/BestSite','w')
for Ligands in tqdm(os.listdir(dire+'/CustomDatas/nrsite')):
	start = time.time()
	bline = open(dire+'/CustomDatas/nrsite/'+Ligands,'r').readlines()
	DistMatrix2t, Arr2, L2t, PCMatrix2t, FullCoord2, FullPDB2, ResList2, CoordKabsch2, AtomKabsch2 = ConstructMatrix(bline)
	CoordKabsch2 = np.asarray(CoordKabsch2, dtype='float')
	
	ResCorres = Run(DistMatrix1t, L1t, PCMatrix1t, DistMatrix2t, L2t, PCMatrix2t, CoordKabsch1t, CoordKabsch2) 
	ResTotal = []
	for i in ResCorres:
		#print(i)
		mcss = [ j for j in i if j[0] > -1 ]
		if len(mcss) > 3:
			#print(mcss)
			mcss = FindConnected(mcss, DistMatrix1t, DistMatrix2t)
			if len(mcss) > 3:
				ResCorres = [ [ResList1[mcss[i][0]], ResList2[mcss[i][1]]] for i in range(len(mcss)) ]
				ResTotal.append(ResCorres)
	
	Ln = PickSequence(ResTotal, CheckmolQuery, CheckMolPDBHash[Ligands])
	end = time.time()
	out.write(Ligands+' '+str(Ln)+' '+str(end-start)+'\n')
out.close()
'''

NrSite = os.listdir(dire+'/CustomDatas/nrsite')
print('Distributing NUMBA codes on {:0.0f} Cores. Please wait.....'.format(nprocs))
print('Compiling the Fragment Scan Module...')

#NrSite = NrSite[:100]

L = len(NrSite)

def parallel_run(Ligands):
	
	bline = open(dire+'/CustomDatas/nrsite/'+Ligands,'r').readlines()
	DistMatrix2t, Arr2, L2t, PCMatrix2t, FullCoord2, FullPDB2, ResList2, CoordKabsch2, AtomKabsch2 = ConstructMatrix(bline)
	CoordKabsch2 = np.asarray(CoordKabsch2, dtype='float')
	ResCorres = Run(DistMatrix1t, L1t, PCMatrix1t, DistMatrix2t, L2t, PCMatrix2t, CoordKabsch1t, CoordKabsch2) 
	ResTotal = []
	for i in ResCorres:
		#print(i)
		mcss = [ j for j in i if j[0] > -1 ]
		if len(mcss) > 3:
			#print(mcss)
			mcss = FindConnected(mcss, DistMatrix1t, DistMatrix2t)
			if len(mcss) > 3:
				ResCorres = [ [ResList1[mcss[i][0]], ResList2[mcss[i][1]]] for i in range(len(mcss)) ]
				ResTotal.append(ResCorres)
	
	Ln = PickSequence(ResTotal, CheckmolQuery, CheckMolPDBHash[Ligands])
	
	out = open(path_out_fold+'/'+Ligands, 'w')
	out.write(Ligands+' '+str(Ln)+'\n')
	out.close()
	
	if random.randint(0,200) < 10:
		L1 = len(os.listdir(path_out_fold))
		print('Total No. of Sites is {:0.0f}. Site that have been processed so far is {:0.0f}. Percent completed : {:3.3f}%'.format(L, L1, (L1/float(L))*100 ), end = '\r')

	
p = Pool(nprocs)
p.map(parallel_run, NrSite)
print('\n')	
out = open(dire+'/'+OutFold+'/BestSite.txt','w')
for i in os.listdir(path_out_fold):
	aline = open(path_out_fold+'/'+i, 'r').readlines()
	out.write(aline[0])
out.close()

#shutil.rmtree(path_out_fold)
	



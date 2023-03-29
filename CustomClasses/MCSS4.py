import numpy as np
import math
import time
from collections import defaultdict, Counter
import copy
import os
import subprocess
import sys
from openbabel import openbabel

dire = os.getcwd()


openbabel.obErrorLog.StopLogging()
obConversion = openbabel.OBConversion()
obConversion.SetInAndOutFormats("pdb", "sdf")


class MCSS:


	def checkmol(self, aline, rnd, Folder):	
		rnd = str(rnd)
		dire_temp = dire+'/'+Folder+'/temp'
		os.chdir(dire_temp)
		out = open(dire+'/'+Folder+"/temp/"+rnd+".pdb", 'w')
		for line in aline:
			if line[:6] == 'HETATM':
				out.write(line)
		out.close()
		het = [ line[6:11].strip() for line in aline if line[:6] == "HETATM"]
		dic = {}
		for i in range(len(het)):
			dic[str(i+1)] = het[i]
			
		#c1 = subprocess.Popen("obabel -ipdb "+rnd+".pdb -osdf > "+rnd+".sdf", shell=True, stdout=subprocess.PIPE)
		#c2 = c1.communicate()[0]
		#print(aline[0].strip())
		mol = openbabel.OBMol()
		obConversion.ReadFile(mol, rnd+'.pdb')
		obConversion.WriteFile(mol, rnd+'.sdf')
		
		c1 = subprocess.Popen("./checkmol -p "+rnd+".sdf", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		c2 = c1.communicate()[0]
		c2 = c2.decode('UTF-8')

		os.chdir(dire)
		check = False
		
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
				#print j, dic[j], dic1[dic[j]] # print end
			#print '\n'	
			
			
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

			
		return check, arr1
	
	
	
	def center_residues(self, arr, coord):
		coord -= np.mean(coord, axis=0)
		dic = {1:"   ", 2:"  ", 3:" ", 4:""}
		brr = []
		for i in range(len(arr)):
			#print arr[i]
			j = [ "%.3f"%j for j in coord[i]]
			val = ''.join([ dic[len(k.split(".")[0])]+k for k in j ])
			#print val
			brr.append(arr[i][:30]+val+arr[i][54:])
			#print '\n' # j[:30]+i+j[54:]
			#time.sleep(1)
		return brr
		
	
	
	def PairWise(self, res, coord, het_atm_label):
		coord = np.asarray(coord, dtype='float')
		dic = {}
		arr = []
		for i in range(len(res)):
			x_ca, y_ca, z_ca = coord[i]
			for j in range(len(res)):
				if i != j:
					dic[res[i]+' '+res[j]] = het_atm_label[res[i]]+' '+het_atm_label[res[j]]
					#time.sleep(1)
					x1_ca, y1_ca, z1_ca = coord[j]
					ans = math.sqrt(pow(( x_ca - x1_ca ),2) + pow(( y_ca - y1_ca ),2) + pow(( z_ca - z1_ca ),2))
					arr.append( [res[i], res[j], ans] )
					
		for i in range(len(res)):
			#print i, res[i], het_atm_label[res[i]]
			dic[res[i]+' '+res[i]] = het_atm_label[res[i]]+' '+het_atm_label[res[i]]
		#time.sleep(11)			
					
		return arr, dic					
			
	
	def file_process(self, arr):
		brr, coord = [], []
		whole_dic, h_dic = {}, {"H":0}
		het_atm_label = {}
		
		for i in arr:
			i = i.strip()
			if i[:6] == "HETATM":
				if i[74:78].strip() not in h_dic:
					var = i[6:11].strip()
					if var not in whole_dic:
						#print i
						het_atm_label[var] = i[71:].strip() 
						#time.sleep(1)
						whole_dic[var] = 0
						brr.append(i)
						coord.append([i[28:38].strip(), i[38:46].strip(), i[46:54].strip()])
	
		coord = np.asarray(coord, dtype='float')
		brr = self.center_residues(brr, coord)
		
		dic1 = defaultdict(list)
		res, coord1 = [], []
		for line in brr:
			if line[:6] == "HETATM":
				val = line[6:11].strip()
				res.append(val)
				coord1.append([line[28:38].strip(), line[38:46].strip(), line[46:54].strip()])
		#print coord1[:3];time.sleep(1)
			
		arr, het_atm_label = self.PairWise(res, coord1, het_atm_label)
		dic, arr1, dic1_pairs = {}, [], defaultdict(list)
		dic, dic1_pairs = defaultdict(list), defaultdict(list)
		arr = sorted(arr, key = lambda x:float(x[2]))
		
		for i in arr:
			#if i[0] == '8213':
			#	print i
			#	time.sleep(.1)
			arr1.append(i[0]+' '+i[1])
			dic[i[0]+' '+i[1]] = i[2]
			dic1_pairs[i[0]].append(i[1])
		#print '\n'
		#for i in arr1:
		#	if i == '8215 8214':
		#		print i
		return arr1, dic, dic1_pairs, brr, het_atm_label
		
		
		
	def SortedArr(self):
		return True
		#print SequenceArrays1
		#print SequenceArrays2
		a1 = []
		if len(SequenceArrays1) <= 3:
			return True
		else:
			for i in zip(SequenceArrays1, SequenceArrays2):
				a1.append(i[0].split(' ')[1]+'-'+i[1].split(' ')[1])
			a1 = sorted(set(a1))
			a1 = ' '.join(a1)
			if a1 not in SortedArrDic:
				SortedArrDic[a1] = 0
				return True
			else:
				return False
				
		
		
	
		
		
	def compare_matrices(self, r1, r2):
		#print  het_atm_label1[r1], het_atm_label2[r2] 
		if het_atm_label1[r1] == het_atm_label2[r2]:
			#print  het_atm_label1[r1], het_atm_label2[r2] # if (mat2[0]-1) < mat1[0] < (mat2[0]+1):
			#print res_dic2[r2], res_dic1[r1]
			if res_dic2[r2]-.25 < res_dic1[r1] < res_dic2[r2]+.25:
				#print r1, r2, res_dic2[r2], res_dic1[r1],'+++	'
				return True
			else:
				return False
		else:
			#print 'het error'
			return False
		return False	
		
	
	
	def CheckDistance(self, S1, S2, v1, v2):	
		#print S1, S2, v1, v2	,'aa'
		if self.compare_matrices(v1, v2):
			#print 'yes'
			v1a = v1.split(' ')[1]
			v2a = v2.split(' ')[1]
			for i in zip(S1, S2):
				p1, p2 = i
				p1 = p1.split(' ')[1]
				p2 = p2.split(' ')[1]
				if p1 == v1a or p2 == v2a:
					return False
				if not self.compare_matrices(p1+' '+v1a, p2+' '+v2a):
					return False
			return True		
		else:
			#print 'nope'
			return False	
			
		
		
		
	def Recursion(self, res_pair1, sequence):	
		if sequence == 'start':
			for i in res_arr2:
				if self.compare_matrices(res_pair1, i):
					#print (i, res_dic1[res_pair1], res_dic2[i] )
					return res_pair1, i, 'First'
			return None, None, 'First'
			
		else:
			res1a, res1b = res_pair1.split(' ')
			res2a, res2b = sequence.split(' ')
			#print dic_single1
			for i in res_pairs_dic1[res1b]:
				#print res1b,i;time.sleep(.1)
				if i != res1a and i not in dic_single1:
					#print 'went'
					for j in res_pairs_dic2[res2b]:
						#print i,j 
						#time.sleep(1)
						if j != res2a and j not in dic_single2:
							#print SequenceArrays1, SequenceArrays2, res1b+' '+i, res2b+' '+j
							Check = self.CheckDistance(SequenceArrays1, SequenceArrays2, res1b+' '+i, res2b+' '+j)
							#time.sleep(.1)
							if Check:
								#print 'wen1'
								if i+'\t'+j not in dic_pair_captch:
									seq1 = ' '.join(SequenceArrays1)+' '+res1b+' '+i
									seq2 = ' '.join(SequenceArrays2)+' '+res2b+' '+j
									if seq1+'\t'+seq2 not in dic_whole:	 
										return res1b+' '+i, res2b+' '+j, 'Next'
								#else:
								#	print (dic pairs)			
			return None, None, 'Next'
			
		
	def PairNext(self, S1, S2):
		dic1, dic2 = {}, {}
		if len(S1) <= 1:
			return None, None, None, None
		for i in zip(S1[:-1], S2[:-1]):	
			dic1[i[0].split(' ')[0]] = 0
			dic1[i[0].split(' ')[1]] = 0
			dic2[i[1].split(' ')[0]] = 0
			dic2[i[1].split(' ')[1]] = 0
		
		S1, S2 = S1[:-1], S2[:-1]
		p1 = S1[-1].split(' ')[1]
		p2 = S2[-1].split(' ')[1]
		
		for i in res_pairs_dic1[p1]:
			if i not in dic1:
				for j in res_pairs_dic2[p2]:
					if j not in dic2:
						Check = self.CheckDistance(S1, S2, p1+' '+i, p2+' '+j)
						if Check:
							seq1 = ' '.join(S1)+' '+p1+' '+i
							seq2 = ' '.join(S2)+' '+p2+' '+j	
							#print(seq1, seq2,'--')
							if seq1+'\t'+seq2 not in dic_whole:
								return(S1, S2, p1+' '+i, p2+' '+j)
		return self.PairNext(S1, S2)		
		
		
		
			
		
		
	def run(self):	
		global dic_single1, dic_single2, SequenceArrays1, SequenceArrays2, dic_pair_captch, dic_whole, SortedArrDic

		dic_loop1, dic_loop2 = {}, {}
		Final1, Final2 = [], []
		BreakLoop = False
		dic_whole, SortedArrDic = {}, {}
		#print res_arr1[:3]
		#print res_arr2[:3]
		for i in res_arr1:
			if res_dic1[i] > 4:
				continue
			if BreakLoop:
				#print 'triger'
				break
			#start = time.time()
			#print i,'--'	
			#if i != '8208 8217':
			#	continue
			#print i;time.sleep(.1)	
			ans = 'start'
			SequenceArrays1, SequenceArrays2 = [], []
			dic_loop1, dic_loop2 = {}, {}
			dic_single1, dic_single2 = {}, {}	
			
			InitiateFirstBreak = False
			ObtainedCount = 0
		
			while True:
				#print 'while1'
				dic_pair_captch = {}
				if ans != 'start':
					if not SequenceArrays1:
						break
				if InitiateFirstBreak:
					if not SequenceArrays1:
						break
					#if ObtainedCount <= 1:
					#	print 'BREAK Obtain'
					#	break	
					if len(SequenceArrays1) <= 2:
						seq1 = ' '.join(SequenceArrays1)
						seq2 = ' '.join(SequenceArrays2)
						dic_whole[seq1+'\t'+seq2] = 0
						#print 'BREAK CALLED'
						break		
					
					ObtainedCount = 0
					#print dic_whole
					#print SequenceArrays1, SequenceArrays2,'------outer while inputs'
					SequenceArrays1, SequenceArrays2, i1, ans = self.PairNext(copy.deepcopy(SequenceArrays1), copy.deepcopy(SequenceArrays2))
					if not i1:
						break	
					SequenceArrays1.append(i1)
					SequenceArrays2.append(ans)
					if i1.split(' ')[0] in dic_single1:	
						del dic_single1[i1.split(' ')[0]]
					if ans.split(' ')[0] in dic_single2:
						del dic_single2[ans.split(' ')[0]]
					#print SequenceArrays1,  i1,'-',ans,'-----outer while'
					
					for j in zip(SequenceArrays1, SequenceArrays2):
						dic_pair_captch[j[0].split(' ')[0]+'\t'+j[1].split(' ')[0]] = 0

					i = i1	
									
				while True:
					#print 'While 2'
					#print dic_single1,'bwefore'
					if not ans:
						#print 'ans1 break'
						break
					#print i, ans,'inputs',	SequenceArrays1
					i, ans, CheckPoint = self.Recursion(i, ans)
					#print i, ans,'--'
					#time.sleep(1)
					ObtainedCount += 1
					if not ans:
						#print 'ans2 break'
						break
					if i+'\t'+ans in dic_whole:
						#print 'dic break'
						break	
						
					InitiateFirstBreak = True	
					dic_single1[i.split(' ')[0]] = 0	
					dic_single2[ans.split(' ')[0]] = 0
					#dic_pair_captch[i.split(' ')[1]+'\t'+ans.split(' ')[1]] = 0
					dic_pair_captch[i.split(' ')[0]+'\t'+ans.split(' ')[0]] = 0
					#print (dic_pair_captch)
					SequenceArrays1.append(i)	
					SequenceArrays2.append(ans)	
					#print dic_single1
					#print SequenceArrays1
				'''
				print len(SequenceArrays1)#, SequenceArrays2
				time.sleep(.1)
				'''
				
				seq1 = ' '.join(SequenceArrays1)
				seq2 = ' '.join(SequenceArrays2)	
				dic_whole[seq1+'\t'+seq2] = 0
				Final1.append(SequenceArrays1)
				Final2.append(SequenceArrays2)
				break
				#if not self.SortedArr():
				#	break
		#print len(Final1)			
		#print max([ len(i) for i in Final1 ])
		#time.sleep(11)
		return Final1, Final2			
					
					
	def process_hits1(self, Final1, Final2):	
		
		Final = []
		for i in range(len(Final1)):
			Final.append([Final1[i], Final2[i], len(Final1[i])])
		Final = sorted(Final, key = lambda x:int(x[2]), reverse=True)
		Final1, Final2 = [], []
		for i in Final:
			Final1.append(i[0])
			Final2.append(i[1])
		
		arr = []
		for i in range(len(Final1)):
			#print Final1[i], len(Final1[i])
			if len(Final1[i]) >= 3:
				#print Final1[i],'--'
				arr1 = []
				for j in zip(Final1[i], Final2[i]):
					#print j[0].split(' ')[1], j[1].split(' ')[1]
					arr1.append(j[0].split(' ')[1]+'-'+j[1].split(' ')[1])
				#print 	Final1[i], Final2[i]
				#if len(arr1) > 3:
				#	arr1 = arr1[:-1]
				#arr1 = arr1[:-1]	
				#print arr1;time.sleep(1);print '\n'
				arr2 = sorted(set(arr1))
				arr2 = ' '.join(arr2)
				arr.append(arr2)
				
				arr1 = arr1[:-1]
				arr2 = sorted(set(arr1))
				arr2 = ' '.join(arr2)
				arr.append(arr2)
				
		arr = sorted(set(arr))		
		arr = sorted(arr, key = lambda x:len(x.split(' ')), reverse=True)
		#print arr
		
		
		arr1 = []
		for i_local in arr:
			iM = '-'.join([ j.split('-')[1] for j in i_local.split(' ') ])+' '+'-'.join([ j.split('-')[0] for j in i_local.split(' ') ]) 
			iM = self.GetConnected(iM)
			if iM:
				for i in iM:
					im1, im2 = i.split(' ')
					im1, im2 = im1.split('-'), im2.split('-')
					if len(im1) >= 3:
						arr1.append(' '.join([ k[1]+'-'+k[0] for k in zip(im1, im2)]))
		arr1 =  sorted(arr1, key = lambda x:len(x.split(' ')), reverse=True)
		return arr1
			
			
			
			
	def process_hits(self, Final1, Final2):
		arr = []
		for i in zip(Final1, Final2):
			#print (i, len(i[0]))
			arr.append([i[0], i[1], len(i[0])])

		arr = sorted(arr, key = lambda x:int(x[2]), reverse=True)
		NewArr = []
		#print len(arr)
		#print arr[14:16]
		for i in arr:
			#print i
			if int(i[2]) > 3:
				#print i
				#time.sleep(11)
				val1 = [ j.split(' ')[1] for j in i[0][:-1] ]
				val2 = [ j.split(' ')[1] for j in i[1][:-1] ]
				NewArr.append([val1, val2, len(val1)])

		#print len(NewArr)
		NewArray = []
		if not NewArr:
			return None
		NewArray.append(NewArr[0]) # base condition
		for i in NewArr:
			check = True
			for j in NewArray:
				dic = { j[0][k]+' '+j[1][k]:0 for k in range(len(j[0])) }
				#print j
				#print [ j[0][k]+' '+j[1][k] for k in range(len(j[0])) ]
				#time.sleep(11)
				dic_count = sum([ 1 for k in range(len(i[0])) if i[0][k]+' '+i[1][k] in dic  ])
				if dic_count == len(i[0]):
					check = False
					break
			if check:
				NewArray.append(i)
		#print len(NewArray)		
		return NewArray
		

	
	# Alignment 
	def rmsd(self, V, W):
		D = len(V[0])
		N = len(V)
		result = 0.0
		for v, w in zip(V, W):
		    result += sum([(v[i] - w[i])**2.0 for i in range(D)])
		return np.sqrt(result/N)


	def kabsch_rmsd(self, P, Q, translate=False):
		P = self.kabsch_rotate(P, Q)
		return self.rmsd(P, Q)


	def kabsch_rotate(self, P, Q):
		U = self.kabsch(P, Q)
		P = np.dot(P, U)
		return P	


	def kabsch(self, P, Q):
		P -= np.mean(P, axis=0)
		Q -= np.mean(Q, axis=0)
		C = np.dot(np.transpose(P), Q)
		V, S, W = np.linalg.svd(C)
		d = (np.linalg.det(V) * np.linalg.det(W)) < 0.0
		if d:
		    #print "reflection detected"
		    S[-1] = -S[-1]
		    V[:, -1] = -V[:, -1]
		# Create Rotation matrix U
		U = np.dot(V, W)
		return U


	def centroid(self, X):
		C = X.mean(axis=0)
		return C
   
   
	def SiteGen(self, check):
		arr = copy.deepcopy(B_all) 
		ans_dic = defaultdict(list)
		ans_arr = []
		if check:
			val = .5
		else:
			val = 1
		val = .5	
		for i in range(0, len(arr)):
			x, y, z = arr[i]
			name1 = het_atm_label1[pdb1_res_info[i]+' '+pdb1_res_info[i]]
			info1 = pdb1_res_info[i]
			for j in range(0, len(site2_coord)):
				name2 = het_atm_label2[pdb2_res_info[j]+' '+pdb2_res_info[j]]
				if name1 == name2:
					x1, y1, z1 = site2_coord[j]
					x_ans = pow((x-x1),2)
					y_ans = pow((y-y1),2)
					z_ans = pow((z-z1),2)
					ans = math.sqrt(x_ans + y_ans + z_ans)
					#print ans;time.sleep(.1)
					if ans < val:
						info2 = pdb2_res_info[j]
						#ans_dic[info1].append(info2)
						ans_arr.append([info1, info2, ans])
					
		ans_arr = sorted(ans_arr, key = lambda x:float(x[2]))		
		dic = {}	
		arr = []	
		for i in ans_arr:
			if i[1] not in dic:
				#print i
				arr.append([i[0], i[1]])
				dic[i[1]] = 0		
		return arr		
		
		
   
	def SiteGen1(self, arr):
	
		#print len(arr), len(site1a), len(site2a)
		arr1 = []
		dic = {1:"   ", 2:"  ", 3:" ", 4:""}
		for i in arr:
			var = ""
			for j in i:
				j1 = "%.3f"%j
				var += dic[len(j1.split(".")[0])]+j1
			arr1.append(var)
	   
		out = open("frag_lig.pdb", 'w')
		#print len(arr1), len(site1a),' ---'
		for i1 in range(len(arr1)):
			i = arr1[i1]
			j = site1a[i1]
			out.write(j[:30]+i+j[54:]+"\n")		
		out.close()	
		
		out = open('fixed_lig.pdb', 'w')
		for i in site2a:
			out.write(i+"\n")
		out.close()	
	   
   
   
	def FilterHits(self, arr, arr1):# New_res, aline_list
		#arr = [ i[0] for i in arr ]
		dic = {}
		for i in arr:
			dic[i[0]] = i[1]
		break_arr = []
		#print arr1,' --arr1s'
		for i in range(0, len(arr1)):
			#print arr1[i]
			checkbreak = False
			for j in arr1[i]:
				if j not in dic:
					#break_arr.append(i)
					checkbreak = True
					break
			if checkbreak:
				break_arr.append(arr1[i]) 		
		#print break_arr,'==='			
					
		dic1 = {}			
		for i in break_arr:
			for j in i:
				dic1[j] = 0
		#print dic1	
		arr2 = []		
		for i in arr:
			if i[0] not in dic1:
				#print i,'--'
				arr2.append(i)
		arr3 = []	
		for i in arr2:
			x, y, z = map(float, pdb1_het_dic[i[0]])
			check = False
			#print i
			for j in arr2:
				if i[0] != j[0]:
					x1, y1, z1 = map(float, pdb1_het_dic[j[0]])
					ans = math.sqrt(pow(( x - x1 ),2) + pow(( y - y1 ),2) + pow(( z - z1 ),2))
					if ans < 2.5:
						check = True
						break
			if check:
				arr3.append(i)		
		#time.sleep(11)
		#print arr3		
		#for i in arr3:
		#	if '8194' in i:
		#		print 'asasas'
		#time.sleep(11)
		return arr3
	

	def GetConnected(self, val):	
		F2, F1 = val.split(' ')
		F1 = F1.split('-')
		F2 = F2.split('-')		
		
		Map = {}
		for i in range(len(F1)):
			Map[F1[i]] = F2[i]
		
		Graph = defaultdict(set)
		for i in F1:
			#print(list(map(float,pdb1_het_dic[i])))
			x, y, z = list(map(float,pdb1_het_dic[i]))
			for j in F1:
				if i != j:
					x1, y1, z1 = list(map(float,pdb1_het_dic[j]))
					ans = math.sqrt(pow(( x - x1 ),2) + pow(( y - y1 ),2) + pow(( z - z1 ),2))
					if ans < 2.2:
						Graph[i].add(j)
						Graph[j].add(i)
		
		#print(F1)
		Result = []
		dic = {}
		for i in F1:
			if i not in dic:
				queue = []
				queue.append(i)
				Answer = set()
				while queue:
					node = queue.pop()
					Answer.add(node)
					for nodes in Graph[node]:
						if nodes not in dic:
							queue.append(nodes)
							dic[nodes] = 0
				#print(i,list(Answer))
				if len(Answer) >= 4:
					Result.append([ list(Answer), len(Answer)  ])
		Result = sorted(Result, key = lambda x:x[1], reverse=True)

		if not Result:
			return None
		
		if len(Result[0][0]) <= 3:
			return None
	
		ResultArr = []
		for R in Result:
			ResultArr.append('-'.join( [ Map[i] for i in R[0] ]  )+' '+'-'.join(R[0]))
		return ResultArr
			


	def file_read(self, aline, bline, rnd, Folder):
		global res_dic1, res_dic2, res_arr1, res_arr2, res_pairs_dic1, res_pair2_dic2
		global dic_loop1, dic_loop2, dic_whole, SortedArrDic
		global dic_single1, dic_single2, site1a, site2a, pdb1_het_dic, pdb2_het_dic
		global pdb1_ln, pdb2_ln, het_atm_label1, het_atm_label2

		global res_arr1, res_arr2, res_dic1, res_dic2, res_pairs_dic1, res_pairs_dic2, B_all
		global pdb1_res_info, pdb2_res_info, site1_coord, site2_coord
		#print aline
		#print bline
		aline_checkmol, aline_list = self.checkmol(aline, rnd, Folder)
		bline_checkmol, bline_list = self.checkmol(bline, rnd, Folder)
		
		#print aline_list
		#print bline_list
		
		res_arr1, res_dic1, res_pairs_dic1, pdb1_lines, het_atm_label1 = self.file_process(aline)
		res_arr2, res_dic2, res_pairs_dic2, pdb2_lines, het_atm_label2 = self.file_process(bline)

		Final1, Final2 = self.run()	
		

		
		#for i in zip(Final1, Final2):
		#	if len(i[0]) > 3:
		#		print i
		#time.sleep(11)
		pdb1_trans_coord, pdb1_res_info, pdb1_generated_coord, pdb1_het_dic = [], [], [], defaultdict(list)
		pdb2_trans_coord, pdb2_res_info, pdb2_generated_coord, pdb2_het_dic = [], [], [], defaultdict(list)
		for line in pdb1_lines:
			if line[:6] == "HETATM":
				val = line[6:11].strip()
				pdb1_res_info.append(val)
				pdb1_generated_coord.append(line)
				pdb1_trans_coord.append([line[28:38].strip(), line[38:46].strip(), line[46:54].strip()])
				pdb1_het_dic[val].append(line[28:38].strip())
				pdb1_het_dic[val].append(line[38:46].strip())
				pdb1_het_dic[val].append(line[46:54].strip())
					
		for line in pdb2_lines:
			if line[:6] == "HETATM":
				val = line[6:11].strip()
				pdb2_res_info.append(val)
				pdb2_generated_coord.append(line)
				pdb2_trans_coord.append([line[28:38].strip(), line[38:46].strip(), line[46:54].strip()])	
				pdb2_het_dic[val].append(line[28:38].strip())
				pdb2_het_dic[val].append(line[38:46].strip())
				pdb2_het_dic[val].append(line[46:54].strip())	
		
		pdb1_trans_coord = np.asarray(pdb1_trans_coord, dtype='float')	
		pdb2_trans_coord = np.asarray(pdb2_trans_coord, dtype='float')
		
		
		NewArray = self.process_hits1(Final1, Final2)
		
		ResLists = []
		maxi, index_ln = 0, 0
		Visit = {}
		if not NewArray:
			return 'SANT'
		for i in NewArray:
			#print i;time.sleep(1)
			#print i
			VisitCheck = False
			visitCount = 0
			for j in i.split(' '):
				if j not in Visit:
					VisitCheck = True
					visitCount += 1
			#if not VisitCheck:
			#	continue
			if visitCount <= 2:
				#print 'c', visitCount
				continue
				
			#print i	
			site1_arr, site2_arr = [], []
			site1_coord, site2_coord = copy.deepcopy(pdb1_trans_coord), copy.deepcopy(pdb2_trans_coord)
			
			for j in i.split(' '):
				site1_arr.append(pdb1_het_dic[  j.split('-')[0]  ])
				site2_arr.append(pdb2_het_dic[  j.split('-')[1]  ])
				
			site1_arr, site2_arr = np.asarray(site1_arr, dtype=float), np.asarray(site2_arr, dtype=float)
			
			U = self.kabsch(copy.deepcopy(site1_arr), copy.deepcopy(site2_arr))	
			B_all = copy.deepcopy(site1_coord)
			B_all -= site1_arr.mean(axis=0)
			B_all = np.dot(B_all, U)
			B_all += site2_arr.mean(axis=0)
			
			
			site1a = pdb1_generated_coord
			site2a = pdb2_generated_coord
			
			new_res = self.SiteGen(True)
			#print new_res,'0000'
			#self.SiteGen1(B_all)
			new_res = self.FilterHits(new_res, aline_list)
			#if len(new_res) > 3:
			#	print new_res,'--',i
			#	self.SiteGen1(B_all)
			#	#time.sleep(11)
			visitCount = 0
			for j in [ '-'.join(j) for j in new_res ]:
				if j not in Visit:
					visitCount += 1
			if visitCount <= 2:
				#print 'v', new_res
				continue
			
			#time.sleep(11)
			#print new_res,'++'
			if len(new_res) > 3:
				#print new_res, i
				ResLists.append(new_res)
				for j in new_res:
					Visit[j[0]+'-'+j[1]] = 0
				#time.sleep(1)
				#print '\n'
		#print 'done'
		#time.sleep(11)
		ResLists = sorted(ResLists, key = lambda x:len(x), reverse=True)
		if not ResLists:
			return None
		dic = {}
		arr1 = []
		
		for i in ResLists:
			#print i,' ----+'
			arr = []
			for j in i:
				if j[0] not in dic:
					#dic[j[0]] = 0 check this later
					arr.append(j)
			if arr:
				#print i,'---', arr
				if len(arr) > 3:
					#print '---'
					arr1.append(arr)
		ans = []			
		for i in arr1:
			#print i	
			a1, a2 = [], []
			for j in i:
				a1.append(j[0])
				a2.append(j[1])
			#ans.append('-'.join(a1)+' '+'-'.join(a2))
			ans.append('-'.join(a2)+' '+'-'.join(a1))
		Final = []
		
		Dic = {}
		#print aline_list, bline_list
		for i in aline_list:
			for j in i:
				Dic[j] = 0
		Dic1 = {}
		for i in bline_list:
			for j in i:
				Dic1[j] = 0
		
		ans1 = []
		for i in ans:
			check = self.GetConnected(i)
			if check:
				for j in check:
					ans1.append(j)
		
		for i in ans1:
			check1 = False
			for j in i.split(' ')[0].split('-'):
				if j in Dic:
					check1 = True

			#check1 = True
			#print '\n'
			if check1:
				dic2 = { j:0 for j in i.split(' ')[0].split('-') }
				dic1 = { j:0 for j in i.split(' ')[1].split('-') }
				#print dic2		
				#print dic1		
				dic1a, dic2a = {}, {}
				for j in res_dic1.items():
					#print j,'--'
					if j[1] < 2 and j[0].split(' ')[1] not in dic1 and j[0].split(' ')[0] in dic1:
						dic1a[j[0].split(' ')[0]] = het_atm_label1[j[0]]
				for j in res_dic2.items():
					#print j;time.sleep(1)
					if j[1] < 2 and j[0].split(' ')[1] not in dic2 and j[0].split(' ')[0] in dic2:
						dic2a[j[0].split(' ')[0]] = het_atm_label2[j[0]]
				#print dic1a, dic2a
				check = True
				#print dic2a
				for j in zip(i.split(' ')[0].split('-'), i.split(' ')[1].split('-')):
					#print j
					if j[0] in dic1a and j[1] not in dic2a:
						check = False
	
					
					if j[0] in dic1a and j[1] in dic2a:
						if dic1a[j[0]] != dic2a[j[1]]:
							check = False
				#print check			
				if check:
					Final.append(i)
				#print '\n'	
			else:
				Final.append(i)	
			
			
		Final1 = copy.deepcopy(Final)
		Final = []
		#print Final1
		#print Dic
		#print Dic1
		for i in Final1:
			i1 = i.split(' ')
			#print i,'00'
			if sum([ 1 for j in i1[0].split('-') if j in Dic1 ]) == sum([ 1 for j in i1[1].split('-') if j in Dic ]):
				#print sum([ 1 for j in i1[0].split('-') if j in Dic1 ]), sum([ 1 for j in i1[1].split('-') if j in Dic ])
				Final.append(i)
		#print Final,'--'
		return Final



'''
mcss = MCSS()

aline = open('ATP.pdb', 'r').readlines() # ref
bline = open('ATP.pdb', 'r').readlines() # fixed

# 3648-3649-3645-3646-3642 2339-2342-2347-2341-2333 5
#aline = open('HEM.pdb', 'r').readlines() # ref
#bline = open('HEM.pdb', 'r').readlines() # fixed
m1 = mcss.file_read(aline, bline, 1000)
print(m1,'ans')
'''








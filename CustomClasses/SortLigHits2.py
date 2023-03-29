import time
import copy
from sklearn import preprocessing
import numpy as np
import copy
import os

dire = os.getcwd()

class SortLigHits:

	def get_maximum(self, aline):
		SitesDic = {}
		for line in aline:
			line = line.strip()
			l = line.split(' ')
			if len(l) > 2:
				nos = len(l[1].split('-'))
				if l[0] not in SitesDic:
					SitesDic[l[0]] = 0
				if l[0] in SitesDic:
					if SitesDic[l[0]] < nos:
						SitesDic[l[0]] = nos
		return SitesDic
			

	def file_read(self, lig_aline, Folder):
		aline = open(dire+'/'+Folder+'/ligand_mcss_hits.txt', 'r').readlines()
		SitesDic = self.get_maximum(aline)
		arr = []
		#print len(aline)
		aline = sorted(set(aline))
		#print len(aline)
		for line in aline:
			line = line.strip()
			l = line.split(' ')
			if len(l) == 3:
				if len(l[1].split('-')) > 3:
					ln = len(l[1].split('-'))+SitesDic[l[0]]
					arr.append(line+"\t"+str(ln))
		#arr = sorted(arr, key = lambda x:int(x.split('\t')[1]))
		arr = sorted(arr, key = lambda x:int(x.split('\t')[1]), reverse=True)
		#print arr[:3]
		arr1 = copy.deepcopy(arr)

		bline = lig_aline
		bline = [ line[6:11].strip() for line in bline if line[:6] == 'HETATM']
		#print aline
		dic = {}
		for i in bline:
			dic[i] = 0
		count = 0

		hit_dic = {}
		for line in arr:
			l = line.split(' ')
			for i in l[1].split('-'):
				if i not in hit_dic:
					hit_dic[i] = 0
				if i in hit_dic:
					val = hit_dic[i]
					val += 1
					hit_dic[i] = val
		
		sum1 = 0
		for i,j in hit_dic.items():
			#print i,j
			sum1 += j
		arr = []
		for i,j in hit_dic.items():
			arr.append([i,j,float(j)/sum1])
		arr = sorted(arr, key = lambda x:float(x[2]), reverse=True)
		
		arr = np.asarray(arr)
		rng = arr[:,2]
		rng = np.asarray(rng, dtype='float').reshape(-1,1)
		#print rng
		min_max_scaler = preprocessing.MinMaxScaler()
		x_scaled = min_max_scaler.fit_transform(rng)
		rng = [round(i[0],3) for i in x_scaled]
		#print rng
		hit_parse = []
		het_dic_cutoff = {}
		
		minim = []
		for i in range(len(arr)):
			#print arr[i], rng[i]
			#print i
			percent = 0
			check = rng[i]
			#print check
			if check >= .8:
				percent = .15
			elif check >= .6:
				percent = .2
			elif check >= .4:
				percent = 0.3
			elif check >= .2:
				percent = 0.4
			elif check >= 0.0:
				#print check
				percent = 0.5
			#print arr[i], rng[i], percent, int(arr[i][1])*percent	
			minim.append(float(int(arr[i][1])*percent))
		minim = sorted(minim)	
		mini = min(minim)
		
		
		#print mini,'+++',sum(minim[:3])/3
		print(mini,' sort lig')
		if mini < 1000: # keep every thing
			rng1 = range(1,100)
			rng1 = rng1[::-1]
			val = 1
			for j in rng1:
				if float(mini)/j > 3000:
					val = j
					break
			if val == 1:
				val = 1	
			#print(arr)
			for i in range(len(arr)):
				percent = 0
				check = rng[i]
				if check >= .8:
					percent = .3
				elif check >= .6:
					percent = .3
				elif check >= .4:
					percent = .4
				elif check >= .2:
					percent = .5
				elif check >= 0.0:
					#print check
					percent = .7
				percent = 1.0
				#print(arr[i][1])
				het_dic_cutoff[arr[i][0]] = int((int(arr[i][1])*percent)/val )
				percent = 0
				
		else:
			rng1 = range(1,100)
			rng1 = rng1[::-1]
			val = 1
			for j in rng1:
				if float(mini)/j > 3000:
					val = j
					break
			if val == 1:
				val = 1		
			
				
			for i in range(len(arr)):
				percent = 0
				check = rng[i]
				if check >= .8:
					percent = .15
				elif check >= .6:
					percent = .15
				elif check >= .4:
					percent = 0.3
				elif check >= .2:
					percent = 0.4
				elif check >= 0.0:
					#print check
					percent = 0.5	
				#print arr[i], rng[i], percent, int(arr[i][1])*percent, (int(arr[i][1])*percent)/val 
				het_dic_cutoff[arr[i][0]] = int((int(arr[i][1])*percent)/val )
			
			
			
				
			#print 'done'
	
		print(het_dic_cutoff)	
		out = open(dire+'/'+Folder+'/top_pick_lengths.txt', 'w')
		for line in arr1:
			count += 1
			line = line.strip()
			l = line.split(' ')
			for i in l[1].split('-'):
				if i in dic:
					val = dic[i]
					val += 1
					dic[i] = val
			#print l
			res1 = l[1].split('-')
			res2 = l[2].split('\t')[0].split('-')
			res1a, res2a = [], []
			for i in range(len(res1)):
				if dic[res1[i]] < het_dic_cutoff[res1[i]]:
					res1a.append(res1[i])
					res2a.append(res2[i])
			#print dic		
			if len(res1a) > 3:		
				#print line
				out.write(l[0]+" "+"-".join(res1a)+" "+"-".join(res2a)+" "+line.split('\t')[1]+'\n')
			#else:
			#	print line	
		out.close()
		
		

'''
sort_lig_hits = SortLigHits()
aline = open('ATP.pdb', 'r').readlines()
sort_lig_hits.file_read(aline, 'ATPSingleOut')
'''





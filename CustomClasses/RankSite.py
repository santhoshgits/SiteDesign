import time
import os
from collections import Counter, defaultdict
import sys
import shutil
import math
from tqdm import tqdm


dire = os.getcwd()


class RankSite():

	def __init__(self,Folder):
		
		aline = open(dire+'/'+Folder+'/final2.pdb', 'r').readlines()
		self.res_coord = defaultdict(list)
		for line in aline:
			line = line.strip()
			if line[:4] == 'ATOM' and line[13:16].strip() == 'CA':
				self.res_coord[line[17:26]].append(float(line[28:38].strip()))
				self.res_coord[line[17:26]].append(float(line[38:46].strip()))
				self.res_coord[line[17:26]].append(float(line[46:54].strip()))
				
				
	def compare(self, s1, s2):
		
		count = 0
		for i in s1[0]:
			x1, y1, z1 = i[0]
			r1 = i[1][:3]
			for j in s2[0]:
				x2, y2, z2 = j[0]
				r2 = j[1][:3]	 	
				if r1 == r2:
					ans = math.sqrt(pow(( x1 - x2 ),2) + pow(( y1 - y2 ),2) + pow(( z1 - z2 ),2))
					if ans < 1:
						count += 1
		return round(count/float(min(len(s2[0]), len(s1[0]))),2)
		
				
	def run(self, aline):	
		
		site_dic = defaultdict(list)
		count = 0 
		for line in aline:
			res = line[0].split('_')
			res_coord = [ self.res_coord[i] for i in res ]
			site_dic[count].append(list(zip(res_coord, res)))
			count += 1
			
		score_dic = {}	
		for s1 in site_dic.keys():
			score = 0
			for s2 in site_dic.keys():
				score += self.compare(site_dic[s1], site_dic[s2])
			score_dic[s1] = score
		
		arr = []
		for i,j in score_dic.items():
			arr.append([str(i)+'.pdb', j])
			
		arr = sorted(arr, key = lambda x:float(x[1]))	
		
		return arr
		
			
			
'''

rank = RankSite('SAM1')

aline = open(dire+'/SAM1/DockResult2','r').readlines()[:50]

Arr = []
for line in aline:
	line = line.strip()
	l = line.split('\t')
	#print(l)
	bline = open(dire+'/SAM1/site/'+l[0], 'r').readlines()
	arr = '_'.join(sorted(set([ i[17:26] for i in bline if i[:4] == 'ATOM' ])))
	Arr.append([arr, 0])

rank.run(Arr)

'''



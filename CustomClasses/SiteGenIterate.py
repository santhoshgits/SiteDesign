import time
import os
#from SiteGen15a import SiteGen
from CustomClasses.SiteGen16 import SiteGen
import sys
from multiprocessing import Pool
from CustomClasses.SiteGenFit import SiteGenFit
from collections import Counter, defaultdict
from CustomClasses.Processed import Processed


dire = os.getcwd()

class SiteGenIter():
	
	'''
	def processed(self, Total):
		#print Total[0], Total[2]
		x, ligand_aline, minim, maxim = Total
		return None
		#return site_gen.file_read(1, ligand_aline, minim, maxim)
	'''
	
	
	def file_read(self, res_point_line, ligand_aline, minim, maxim, total_no, nproc, Folder):
		process = Processed()
		site_gen_fit = SiteGenFit()
		total_no1 = total_no*3

		#site_gen.file_read(1, ligand_aline, 20, 29)
		arr = []
		ans = []
		#ligand_aline = '\n'.join(ligand_aline)
		p = Pool(nproc)
		for i in range(total_no1):
			arr.append([1, ligand_aline, minim, maxim, Folder])		
		
		ans = p.map(process.file_read, arr)
		#print (ans)
		Final = []
		for i in ans:
			#print(i, 'site gen iter')
			#if maxim >= len(i) > minim:
			v1, v2 = site_gen_fit.file_read(i, Folder)
			Final.append([i, v1, v2])
		Final = sorted(Final, key = lambda x:int(x[1]), reverse=True)	
		#Final = sorted(Final, key = lambda x:int(len(x[0])), reverse=True)
		Final = Final[:int(total_no1/2)]
		Final = sorted(Final, key = lambda x:int(x[2]), reverse=True)
		#for i in Final[:total_no]:
		#	print i,'\n'
		#return Final[:total_no]
		return Final

'''
site_gen_iter = SiteGenIter()
res_point_line = open('residue_position', 'r').readlines()	
ligand_aline = open('lig.pdb', 'r').readlines()	
minim, maxim = 12, 18
ans = site_gen_iter.file_read(res_point_line, ligand_aline, minim, maxim, 20, 4)
print (ans)
'''
'''
dic = defaultdict(list)
dire = os.getcwd()
aline = open('final2.pdb', 'r').readlines()
for line in aline:
	line = line.strip()
	dic[line[17:26]].append(line)
	
for i in range(len(ans)):
	out = open(dire+'/temp1/'+str(i)+'a.pdb', 'w')
	#print ans[i][-1]
	for j in ans[i][0]:
		#print j
		for k in dic[j]:
			out.write(k+'\n')
	#print '\n'		
	out.close()		

'''







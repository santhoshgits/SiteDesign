import numpy as np
import random
import os
import time
import shutil



dire = os.getcwd()


	
class RenameSite():

	
	def __init__(self):
		self.alphabet = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O",\
		   "P", "Q", "R", "S", "T", "U", "W", "X", "Y", "Z" ]

		self.rng_nos = range(1,9999)	

	def b(self, aline):
		arr = []
		for line in aline:
			if line[:4] == "ATOM":
				arr.append(line[17:26])
		arr = sorted(set(arr))
		return arr

	def c(self, aline, dic):
		arr = []
		for line in aline:
			#line = line.strip()
			if line[:4] == "ATOM":
				#print line
				arr.append(line[:17]+dic[line[17:26]]+line[26:])
				#print "\n"
			#if line[:6] == "HETATM":
			#	arr.append(line)
		return arr	
		
	def file_read(self, Folder):
		dic_nos = {}
		dic_space = {1:"   ", 2:"  ", 3:" ", 4:""}
		total_ln = len(os.listdir(dire+'/'+Folder+"/matched_hits"))
		count1 = 0
		#for i in os.listdir(dire+"/matched_hits"):
		bline = open(dire+'/'+Folder+'/align_site_mapper.txt', 'r').readlines()
		outf = open(dire+'/'+Folder+'/StrideMap.txt','w')
		for line in bline:
			i = line.strip()
			i = i.split(' ')
			#print i
			aline = open(dire+'/'+Folder+"/matched_hits/"+i[-1], 'r').readlines()
			arr = self.b(aline)
			brr = []
			count1 += 1
			#print count1," ------ ", total_ln
			dic_residue = {}
			#print arr," arr"
			for j in arr:
				for k in self.alphabet:
					count = 0
					for l in self.rng_nos:
						var = j[:4]+k+" "+str(l)
						if var not in dic_nos:		
							dic_nos[var] = 0
							brr.append(k+dic_space[len(str(l))]+str(l))
							count = 1
							break
					if count == 1:
						break
			if arr:		
				#print (arr,brr)
				for j in zip(arr,brr):
					dic_residue[j[0]] = j[0][:4]+j[1]
				#print(i[0],dic_residue)	
				
				crr = []
				for j in dic_residue.items():
					crr.append('+'.join(j))
				crr = '_'.join(crr)
				outf.write(i[0]+'\t'+crr+'\n')	
				
				final = self.c(aline, dic_residue)
				out = open(dire+'/'+Folder+"/matched_hits_renamed/"+i[-1], 'w')
				out.write("".join(final))
				out.close()
		outf.close()		
				
'''
rename_site = RenameSite()
rename_site.file_read('ZITOutput')
'''
	
	
		





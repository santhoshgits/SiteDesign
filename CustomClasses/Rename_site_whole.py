import numpy as np
import random
import os
import time
import shutil



dire = os.getcwd()


	
class RenameSiteWhole():

	
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
		#total_ln = len(os.listdir(dire+"/matched_hits"))
		count1 = 0
		arr = []
		try:
			with open(dire+'/'+Folder+"/residue_overlap_removed.pdb") as fh:
				while True:
					line = next(fh)
					if line[13:16] == "CA ":
						arr.append(line[17:26])
		
		except StopIteration:
			pass
		#print arr[:10]
		dic_renamed = {}
		
		for j in arr:
			for k in self.alphabet:
				count = 0
				for l in self.rng_nos:
					var = k+" "+str(l)
					if var not in dic_nos:		
						dic_nos[var] = 0
						#print j[:3]+" "+k+dic_space[len(str(l))]+str(l)
						dic_renamed[j] = j[:3]+" "+k+dic_space[len(str(l))]+str(l)
						#print j
						count = 1
						#print "\n"
						#time.sleep(1)
						break
				if count == 1:
					break	
		
		dic_rename = {}
		aline = open(dire+'/'+Folder+'/residue_frequencies.txt', 'r').readlines()
		for line in aline:
			line = line.strip()
			l = line.split('\t')
			dic_rename[l[0]] = l[1]
		
		
		out = open(dire+'/'+Folder+'/reside_rename_map.txt', 'w')
		#print dic_rename
		for i,j in dic_renamed.items():
			#if i in dic_rename:
			out.write(i+'\t'+j+'\t'+dic_rename[i]+'\n')
			#else:
			#	print i,j
		out.close()
		
		out = open(dire+'/'+Folder+"/final_residue_pool.pdb", 'w')
		try:
			with open(dire+'/'+Folder+"/residue_overlap_removed.pdb") as fh:
				while True:
					line = next(fh)
					#print line
					out.write(line[:17]+dic_renamed[line[17:26]]+line[26:])
					#print "\n"
		except StopIteration:
			pass	
		out.close()
		
'''		
rename_site_whole = RenameSiteWhole()
rename_site_whole.file_read()
'''
	
		





import os
import time
import numpy as np
from collections import defaultdict
import math

dire = os.getcwd()

class ResidueStretch():
	
	def ligand_process(self, lig_line):
		het_dic = defaultdict(list)
		for line in lig_line:
			line = line.strip()
			#print line
			#print line[6:11]
			# append([  float(line[27:38].strip()), float(line[38:46].strip()), float(line[46:54].strip()) ])
			het_dic[line[6:11].strip()].append(float(line[27:38].strip()))
			het_dic[line[6:11].strip()].append(float(line[38:46].strip()))
			het_dic[line[6:11].strip()].append(float(line[46:54].strip()))
		return het_dic	

	
	def association(self, aline, het_dic_coord):
		dic = defaultdict(list)		
		for line in aline:
			line = line.strip()
			dic[line[17:26]].append([  float(line[27:38].strip()), float(line[38:46].strip()), float(line[46:54].strip()) ])
		arr = []
		for i in het_dic_coord.items():
			#print i[0]
			x, y, z = i[1]
			for j in dic.items():
				ln = len(j[1])
				count = 0
				for k in j[1]:
					x1, y1, z1 = k
					x_ans = pow(( x - x1 ),2)
					y_ans = pow(( y - y1 ),2)
					z_ans = pow(( z - z1 ),2)
					ans = math.sqrt(x_ans + y_ans + z_ans)
					if ans < 5:
						count += 1
				if float(count)/ln > .39:
					#print 'done'
					#print j[0],i[0],'done', count, ln,'--',float(count)/ln
					arr.append(i[0]+'-'+j[0])
				#else:	
				#	if count > 1:
				#		print j[0],i[0],'done', count, ln, float(count)/ln
				#time.sleep(11)
		return arr		
	
	def file_read(self, lig_line, Folder):
		het_dic_coord = self.ligand_process(lig_line)
		#out = open(dire+'/'+Folder+'/AssociationOutput', 'w')
		for i in os.listdir(dire+'/'+Folder+'/matched_hits_renamed1'):
			#print (i)
			aline = open(dire+'/'+Folder+'/matched_hits_renamed1/'+i, 'r').readlines()
			result1 = "_".join(sorted(set([j[17:26] for j in aline if j[:4] == 'ATOM'])))
			result = self.association(aline, het_dic_coord)
			#for j in result:
			#	out.write(j+'-'+result1+'\n')
			#print 'done'
			#time.sleep(11)
		#out.close()
		
	
'''	
residue_stretch = ResidueStretch()
ligand_aline = open('ATP.pdb', 'r').readlines()
residue_stretch.file_read(ligand_aline)
'''





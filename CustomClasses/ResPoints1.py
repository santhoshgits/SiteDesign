import time
import numpy as np
import os
from collections import Counter, defaultdict
import math

dire = os.getcwd()

class ResPoints:
	
	
	def __init__(self, ligand_aline, Folder):
		self.res_dict = {'GLY':'G' ,'ALA':'A' ,'LEU':'L' ,'ILE':'I' ,'TRP':'W' ,'SER':'S' ,'THR':'T' ,'TYR':'Y' ,'PHE':'F' ,\
						 'PRO':'P' ,'ASP':'D','GLU':'E' ,'HIS':'H' ,'CYS':'C' ,'MET':'M' ,'VAL':'V' ,'ASN':'N' ,\
						 'GLN':'Q' ,'ARG':'R' ,'LYS':'K' }
		self.lig_coord_line = defaultdict(list)	
		for line in ligand_aline:
			line = line.strip()
			if line[:6] == 'HETATM':
				#[line[28:38].strip(), line[38:46].strip(), line[46:54].strip()]
				self.lig_coord_line[line[6:11].strip()].append(float(line[28:38].strip()))
				self.lig_coord_line[line[6:11].strip()].append(float(line[38:46].strip()))
				self.lig_coord_line[line[6:11].strip()].append(float(line[46:54].strip()))	
				
		aline = open(dire+'/'+Folder+'/reside_rename_map.txt', 'r').readlines()		
		self.res_pairs = {}
		for line in aline:
			line = line.strip()
			l = line.split('\t')
			self.res_pairs[l[0]] = 0
		#print 'done'
		#time.sleep(11)
		
		aline = open(dire+'/'+Folder+'/final2.pdb', 'r').readlines()
		self.res_coord = defaultdict(list)
		for line in aline:
			line = line.strip()
			if line[17:20] in self.res_dict:
				self.res_coord[line[17:26]].append([line[28:38].strip(), line[38:46].strip(), line[46:54].strip()])	

			
	
	
	def DistMeasure(self, dic_coord):
		arr = []
		dic = defaultdict(list)
		dic1 = defaultdict(list)
		for i,j in self.lig_coord_line.items():
			x,y,z = j
			for i1, j1 in dic_coord.items():
				x1, y1, z1 = j1
				k = "_".join(map(str, j1))
				x_ans = pow(( x - x1),2)
				y_ans = pow(( y - y1),2)
				z_ans = pow(( z - z1),2)
				ans = math.sqrt(x_ans + y_ans + z_ans)
				if ans < 6:
					#print ans, i1, iILE A 505
					dic[i].append(i1)
					dic1[i].append(k)
					#arr.append(str(ans)+'\t'+)
		#print dic1			
		return dic, dic1
	
	
	def angle_measure(self, arr):
		new_arr = []
		#angle_dic = defaultdict(list)
		angle_dic = {}
		for i in arr:
			angle_list = []
			for j in i[1].split('\t'):
				#print i[0]
				coord = np.asarray(self.res_coord[j], dtype='float')
				#print(coord)
				check = True
				try:
					ca = coord[1]
				except:
					check = False
					pass	
				if not check:
					continue	
				coord_mean = coord.mean(axis=0)
				# angles --> ca-lig-cn
				a = ca
				b = self.lig_coord_line[i[0]]
				c = coord_mean
				
				ba = a - b
				bc = c - b
				cosine_angle = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
				angle = np.arccos(cosine_angle)
				angle = np.degrees(angle)			
				angle_list.append(round(angle,1))
			angle_list = map(str, angle_list)	
			angle_dic[i[0]] = ' '.join(angle_list)
			new_arr.append('\t'.join(i)+'\t'+' '.join(angle_list))
		return angle_dic	
				
	
	def file_read(self, Folder):
		aline = open(dire+'/'+Folder+'/align_site_mapper.txt', 'r').readlines()		
		out = open(dire+'/'+Folder+'/residue_position.txt', 'w')
		for line in aline:
			line = line.strip()
			l = line.split(' ')
			#print l
			#time.sleep(1)
			if os.path.isfile(dire+'/'+Folder+'/matched_hits_renamed/'+l[-1]):
				pdb_line = open(dire+'/'+Folder+'/matched_hits_renamed/'+l[-1], 'r').readlines()
				dic = {}
				dic_coord = defaultdict(list)
				#print l[-1]
				for i in pdb_line:
					i = i.strip()
					#print i
					if i[:4] == "ATOM" and i[12:15].strip() == "CA" and i[17:26] in self.res_pairs:
						if i[17:26] not in dic:
							#print i
							dic[i[17:26]] = 0
							coord = [i[28:38].strip(), i[38:46].strip(), i[46:54].strip()]
							coord = list(map(float, coord))
							coord = list(map(int, coord))
							#print i[17:26], self.res_dict[i[17:26][:3]], coord
							dic_coord[i[17:26]].append(coord[0])
							dic_coord[i[17:26]].append(coord[1])
							dic_coord[i[17:26]].append(coord[2])
					#else:
					#	print i
				lig_dic, lig_dic1 = self.DistMeasure(dic_coord)
				#print lig_dic
				angle_list = []
				for i in lig_dic:
					#print i
					angle_list.append([i, '\t'.join(lig_dic[i])])
				#print angle_list
				angle_dic = self.angle_measure(angle_list)
				
				for i in lig_dic:
					#print i, lig_dic[i]
					out.write(i+"\t"+"-".join(lig_dic[i])+"\t"+" ".join([ self.res_dict[j[:3]] for j in lig_dic[i] if j[:3] in self.res_dict ])+"\t"+" ".join(lig_dic1[i])+'\t'+angle_dic[i]+'\n')
					#print i+"\t"+"-".join(lig_dic[i])+"\t"+" ".join([ self.res_dict[j[:3]] for j in lig_dic[i] if j[:3] in self.res_dict ])+"\t"+" ".join(lig_dic1[i])+'\t'+angle_dic[i]
					#out.write(i+"\t"+"-".join(lig_dic[i])+"\t"+" ".join([ self.res_dict[j[:3]] for j in lig_dic[i]])+"\t"+" ".join(lig_dic1[i])+'\n')
					#print '\n'
				#print 'done\n'
				#time.sleep(1)
			
		out.close()
		
		
'''
aline = open('MF8.pdb', 'r').readlines()
res_points = ResPoints(aline,'mf8')
res_points.file_read('mf8')
'''







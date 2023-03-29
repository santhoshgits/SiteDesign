import math
import os
import numpy as np
from collections import defaultdict
import time
import copy

dire = os.getcwd()

class AlignSites:

	def ligand_frag(self, aline):
		ligand_atom_dic = defaultdict(list)
		ligand_atom, ligand_coord = [], []
		for line in aline:
			line = line.strip()
			#print line
			ligand_atom.append(line[6:11].strip())
			ligand_coord.append([line[28:38].strip(), line[38:46].strip(), line[46:54].strip()])
		#print ligand_coord	
		ligand_coord = np.asarray(ligand_coord, dtype='float')	
		for i in range(len(ligand_atom)):
			ligand_atom_dic[ligand_atom[i]].append(ligand_coord[i])
		return ligand_atom_dic
	
	
	def kabsch(self, P, Q):
		#print P, Q
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
	
	def site_gen(self, arr, target_atom, Folder):
		arr1 = []
		dic = {1:"   ", 2:"  ", 3:" ", 4:""}
		for i in arr:
			var = ""
			for j in i:
				j1 = "%.3f"%j
				var += dic[len(j1.split(".")[0])]+j1
			arr1.append(var)
		arr2 = []	
		for i in range(len(arr1)):
			arr2.append(target_atom[i][:30]+arr1[i]+target_atom[i][54:])
		if os.path.isfile(dire+'/'+Folder+"/a.pdb"):
			os.system("rm "+dire+'/'+Folder+"/a.pdb")
		return arr2
	
	
	
	
	def CheckClashes(self, sites, bline, ligand_atom_dic, ligand_aline):
		prot_coord_dic = defaultdict(list)
		het_atom_dic = defaultdict(list)
		het_atom_line = defaultdict(list)
		prot_atom_line = defaultdict(list)
		for line in sites:
			line = line.strip()
			if line[:6] == 'HETATM':
				#print line,'------'
				het_atom_dic[line[6:11].strip()].append([line[28:38].strip(), line[38:46].strip(), line[46:54].strip()])
				het_atom_line[line[6:11].strip()].append(line)
			if line[:4] == "ATOM":
				prot_coord_dic[line[17:26]].append([line[28:38].strip(), line[38:46].strip(), line[46:54].strip()])
				prot_atom_line[line[17:26]].append(line)
				#target_coord.append([line[28:38].strip(), line[38:46].strip(), line[46:54].strip()])
				#target_atom.append(line)
		remove_res = {}		
		for i in prot_coord_dic.items():
			#print i[0]
			#print i[1]
			i1 = np.asarray(i[1], dtype='float')
			minim = []
			for j1 in i1:
				x, y, z = j1
				for j in ligand_aline:
					j = j.strip()
					x1, y1, z1 = j[28:38].strip(), j[38:46].strip(), j[46:54].strip()
					x1, y1, z1 = float(x1), float(y1), float(z1)
					ans = math.sqrt(pow(( x - x1 ),2) + pow(( y - y1 ),2) + pow(( z - z1 ),2))
					minim.append(ans)
			#if min(minim) < 2.45:
			if min(minim) < 2.3:
				remove_res[i[0]] = 0
			
		coord = [ het_atom_dic[i][0] for i in bline.split(' ')[2].split('-') ]
		coord = np.asarray(coord, dtype='float')
		#out = open('het1.pdb', 'w')
		#out.write('\n'.join([ het_atom_line[i][0] for i in bline.split(' ')[2].split('-') ]))
		#out.close()
		#time.sleep(11)
		arr = []
		site = []
		for i in prot_coord_dic.items():
			#if i[0] != 'GLY A 397':
			#	continue
			i1 = np.asarray(i[1], dtype='float')
			minim = []
			count = 0
			for j1 in i1:
				x, y, z = j1
				for j in coord:
					#print j
					x1, y1, z1 = j
					ans = math.sqrt(pow(( x - x1 ),2) + pow(( y - y1 ),2) + pow(( z - z1 ),2))
					#if ans < 4.5:
					#	print 
					minim.append(ans)
					if ans < 4.1:
						count += 1
					if ans < 3.5:
						count += 2
					if ans < 3:
						count += 2
			#minim = sorted(minim)		
			#print minim[:5]
			if count >= 3:
				if i[0] not in remove_res:
					arr.append(i[0])
					
		arr = sorted(set(arr))
		for i in arr:
			#print i
			for j in prot_atom_line[i]:
				site.append(j)
		check = False
		if site:
			check = True
			
		#out = open('prot.pdb', 'w')
		#for i in site:
		#	out.write(i)
		#out.close()
		for i in bline.split(' ')[2].split('-'):
			site.append(het_atom_line[i][0])
			
		return site, check
		
		
		
	
	
	def prot_coord(self, prot_aline, bline, ligand_atom_dic, ligand_aline, Folder):
		lig_dir = dire+'/CustomDatas/nrsite'
		l = bline.split(' ')
		aline = open(lig_dir+"/"+l[0], 'r').readlines()
		het_atom_dic = defaultdict(list)
		target_coord, target_atom = [], []
		for line in prot_aline:
			line = line.strip()
			if line[:6] == 'HETATM':
				#print line
				het_atom_dic[line[6:11].strip()].append([line[28:38].strip(), line[38:46].strip(), line[46:54].strip()])
			if len(line) > 2:
				target_coord.append([line[28:38].strip(), line[38:46].strip(), line[46:54].strip()])
				target_atom.append(line)
		prot_het = l[2].split('-')
		query_het = l[1].split('-')
		#print prot_het
		#print query_het
		prot_het_coord, query_het_coord = [], []
		#print het_atom_dic.keys()
		#print prot_het
		##print('\n')
		for i in range(0, len(prot_het)):
			prot_het_coord.append(het_atom_dic[prot_het[i]][0])
			query_het_coord.append(ligand_atom_dic[query_het[i]][0])
		coord1 = np.asarray(prot_het_coord, dtype='float')
		coord2 = np.asarray(query_het_coord, dtype='float')
		target_coord = np.asarray(target_coord, dtype='float')
		
		
		coord1_cnt = coord1.mean(axis=0)
		coord2_cnt = coord2.mean(axis=0)
		coord1 -= coord1_cnt
		coord2 -= coord2_cnt
		
		target_coord -= coord1_cnt
		U = self.kabsch(coord1, coord2)
		U = np.dot(target_coord, U)
		U += coord2_cnt
		
		sites = self.site_gen(U, target_atom, Folder)
		
		
		sites, check = self.CheckClashes(sites, bline, ligand_atom_dic, ligand_aline)
		
		
		out = open(dire+'/'+Folder+'/a.pdb', 'w')
		for i in sites:
			out.write(i+'\n')
		out.close()	
		
		return sites, check
		
		
		
	
	
	
	def file_read(self, ligand_aline, Folder):
		aline = open(dire+'/'+Folder+'/top_pick_lengths.txt', 'r').readlines()
		ligand_atom_dic = self.ligand_frag(ligand_aline)
		lig_dir = dire+'/CustomDatas/nrsite'
		out1 = open(dire+'/'+Folder+'/align_site_mapper.txt', 'w')
		out_count = 0
		for line in aline:
			line = line.strip()
			l = line.split(' ')
			#print l
			out_count += 1
			prot_aline = open(lig_dir+"/"+l[0], 'r').readlines()
			prot_aline = [ i for i in prot_aline if i[74:].strip() != 'H' ]
			#print prot_aline
			sites, site_check = self.prot_coord(prot_aline, line, ligand_atom_dic, ligand_aline, Folder)	
			
			if site_check:
				out1.write(line+' '+str(out_count)+'_'+l[0]+"\n")
				out = open(dire+'/'+Folder+"/matched_hits/"+str(out_count)+'_'+l[0], 'w')
				#print str(out_count)+'_'+l[0]
				for i in sites:
					out.write(i+"\n")
				out.close()	
			
			
			
		out1.close()	
		
		
		
		
		
		
		

'''
align_sites = AlignSites()
ligand_aline = open('ATP.pdb', 'r').readlines()
ligand_aline = [line for line in ligand_aline if line[:6] == 'HETATM']
align_sites.file_read(ligand_aline, 'ATPout')
'''





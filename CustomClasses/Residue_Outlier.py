import numpy as np
import time
import math
import os

dire = os.getcwd()

class ResidueOutlier:
	
	
	def compare_coord(self, lig_coord, aline, Folder):
		coord = []
		out = open(dire+'/'+Folder+"/a.pdb", 'w')
		prot_id = aline[0][17:26]
		for line in aline:
			#if line[:6] == "HETATM":
			coord.append([line[28:38].strip(), line[38:46].strip(), line[46:54].strip()])
			out.write(line)
		out.close()	
		coord = np.asarray(coord, dtype="float")
		ans = []
		#print len(coord), len(lig_coord)
		for i in coord:
			x, y, z = i
			for j in lig_coord:
				x1, y1, z1 = j
				x_ans = pow(( x - x1),2)
				y_ans = pow(( y - y1),2)
				z_ans = pow(( z - z1),2)
				ans.append(math.sqrt(x_ans + y_ans + z_ans))
		if ans:
			if min(ans) < 2.5:
				return True
			if min(ans) > 5.5:
				return True
			return False
		else:
			return True
		
	
	def file_read(self, prot_collated, lig, Folder):
		arr, dic, lig_coord = [], {}, []
		count = 0
		for line in lig:
			if line[:6] == "HETATM":
				lig_coord.append([line[28:38].strip(), line[38:46].strip(), line[46:54].strip()])
		lig_coord = np.asarray(lig_coord, dtype="float")		
		#print lig_coord		
		#time.sleep(11)
		delete_residue = []
		out = open(dire+'/'+Folder+"/residue_overlap_removed.pdb", "w")
		try:
			with open(dire+'/'+Folder+'/'+prot_collated) as fh:
				while True:
					line = next(fh)
					if line[:4] == "ATOM":
						#dic[] = 0
						#print line
						arr.append(line)
						if line[17:26] not in dic:
							count += 1
							dic[line[17:26]] = 0
						if count > 1:
							#print arr[:-1],"\n\n"
							if arr[:-1]:
								check = self.compare_coord(lig_coord, arr[:-1], Folder)
								if check:
									delete_residue.append(arr[0][17:26])
								else:
									out.write("".join(arr[:-1]))
								val = arr[-1]
								arr = []
								arr.append(val)
								dic = {}
								count = 0
							#time.sleep(1)
		except StopIteration:
			pass	
		out.close()
		#print delete_residue
		
'''		
residue_outlier = ResidueOutlier()

aline = open("ATP.pdb", 'r').readlines()
bline = "collated.pdb"

residue_outlier.file_read(bline, aline, 'ATPout')
'''


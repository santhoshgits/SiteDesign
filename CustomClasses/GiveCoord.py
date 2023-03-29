from collections import defaultdict
import os
# ['LEU A1547', 'SER C7027', 'VAL A3354', 'ARG B2923', 'ALA B5428', 'PRO A4832', 'SER C9785', 'ARG A4956']

dire = os.getcwd()

class GiveCoord:
	
	def __init__(self, Folder):
		aline = open(dire+'/'+Folder+'/final2.pdb', 'r').readlines()
		self.res_line_dic = defaultdict(list)
		for line in aline:
			self.res_line_dic[line[17:26]].append(line)
		
	def file_read(self, arr):
		arr1 = []
		for i in arr:
			for j in self.res_line_dic[i]:
				arr1.append(j)
		arr1 = "".join(arr1)
		return arr1
		
'''	
give_coord = GiveCoord()
arr = ['LEU A1547', 'SER C7027', 'VAL A3354', 'ARG B2923', 'ALA B5428', 'PRO A4832', 'SER C9785', 'ARG A4956']
give_coord.file_read(arr)
'''



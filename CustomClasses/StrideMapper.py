import time
import os

dire = os.getcwd()



class StrideMapper():

	
	def stride_dict(self, aline):
		dic = {}
		for line in aline:
			line = line.strip()
			l = line.split('\t')
			dic[l[0]] = l[1:]
		return dic

	def file_read(self, Folder):
		aline = open(dire+'/'+Folder+'/StrideMap.txt', 'r').readlines()
		out = open(dire+'/'+Folder+'/StrideMap1.txt', 'w')
		for line in aline:
			line = line.strip()
			l = line.split('\t')
			if os.path.isfile(dire+'/CustomDatas/StrideOut/'+l[0][:4]+'.pdb'):
				stride_line = open(dire+'/CustomDatas/StrideOut/'+l[0][:4]+'.pdb', 'r').readlines()			
				dic = self.stride_dict(stride_line)	
				for i in l[1].split('_'):
					i1 = i.split('+')
					#print(i1)
					if i1[0] in dic:
						out.write(i1[1]+'\t'+'\t'.join(dic[i1[0]])+'\n')
						#print(i1[1]+'\t'+'\t'.join(dic[i1[0]]))
		out.close()			
				
		
'''
stride_mapper = StrideMapper()
stride_mapper.file_read('ZITOutput')
'''



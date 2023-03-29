import time
import os

dire = os.getcwd()

class GeometryCoordinate():

	def stride(self, Folder):
		aline = open(dire+'/'+Folder+'/StrideMap1.txt', 'r').readlines()
		dic = {}
		for line in aline:
			line = line.strip()
			l = line.split('\t')
			dic[l[0]] = l[1:]
		return dic
		
	def file_read(self, Folder):
		stride_dic = self.stride(Folder)
		
		aline = open(dire+'/'+Folder+'/final2.pdb', 'r').readlines()
		res_dict = {}
		for line in aline:
			line = line.strip()
			if line[13:16].strip() == 'CA':
				val = [line[28:38].strip(), line[38:46].strip(), line[46:54].strip()]
				val = list(map(float, val))
				val = list(map(int, val))
				val = list(map(str, val))
				val = '_'.join(val)
				res_dict[line[17:26]] = val
				
		
		aline = open(dire+'/'+Folder+'/StrideMap.txt', 'r').readlines()
		
		out = open(dire+'/'+Folder+'/StrideMap2.txt', 'w')
		for line in aline:
			line = line.strip()
			l = line.split('\t')
			coord, ss, phi, psi, res = [], [], [], [], []
			for i in l[1].split('_'):
				i1 = i.split('+')
				if i1[1] in stride_dic:
					coord.append(res_dict[i1[1]])
					res.append(i1[1][:3])
					ss.append(stride_dic[i1[1]][0])
					phi.append(stride_dic[i1[1]][1])
					psi.append(stride_dic[i1[1]][2])
			out.write(l[0]+'\t'+'/'.join(coord)+'\t'+'/'.join(res)+'\t'+'/'.join(ss)+'\t'+'/'.join(phi)+'\t'+'/'.join(psi)+'\n')
		out.close()	
		
'''	
geomCoord = GeometryCoordinate()
geomCoord.file_read('ZITOutput')
'''



import os
import time
import shutil

dire = os.getcwd()


class RenameSecond():
	
		
	
	def file_read(self, Folder):
		aline = open(dire+'/'+Folder+'/final_residue_pool.pdb', 'r').readlines()
		arr = []
		for line in aline:
			if line[:4] == 'ATOM' and line[13:15] == 'CA':
				arr.append(line[17:26])
		final_res = sorted(set(arr))
		#print final_res[:10]
		aline = open(dire+'/'+Folder+'/reside_rename_map.txt','r').readlines()
		res_change = {}
		for line in aline:
			line = line.strip()
			l = line.split('\t')
			#print l
			res_change[l[0]] = l[1]
			
		path_out = dire+'/'+Folder+'/matched_hits_renamed1'
		if not os.path.exists(path_out):
			os.mkdir(path_out)
		else:
			shutil.rmtree(path_out)
			os.mkdir(path_out)	

		
		for i in os.listdir(dire+'/'+Folder+'/matched_hits_renamed'):
			#print (i)
			aline = open(dire+'/'+Folder+'/matched_hits_renamed/'+i, 'r').readlines()
			out = open(dire+'/'+Folder+'/matched_hits_renamed1/'+i, 'w')
			for line in aline:
				if line[:4] == "ATOM":
					if line[17:26] in res_change:
						#print res_change[line[17:26]]
						#print line.strip()
						out.write(line[:17]+res_change[line[17:26]]+line[26:])
			out.close()			
			#print 'done'
			#time.sleep(11)
		
		
	
#rename_second = RenameSecond()	

#rename_second.file_read()





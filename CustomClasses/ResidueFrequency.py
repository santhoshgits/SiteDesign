import numpy as np
import random
import os
import time

dire = os.getcwd()

class ResidueFrequency:
	
	
	def file_read(self, Folder):
		aline = open(dire+'/'+Folder+'/align_site_mapper.txt', 'r').readlines()
		count = 0
		out = open(dire+'/'+Folder+'/residue_frequencies.txt', 'w')
		#print 'aaaa', aline
		for line in aline:
			line = line.strip()
			l = line.split(' ')
			#print l
			count += 1
			val = str(count)+'_'+l[0]
			val = l[-1]
			if os.path.isfile(dire+'/'+Folder+"/matched_hits_renamed/"+val):
				bline = open(dire+'/'+Folder+"/matched_hits_renamed/"+val, 'r').readlines()
				#print bline
				for k in sorted(set([ j[17:26] for j in bline if j[:4] == 'ATOM'])):
					out.write(k+'\t'+l[-2]+"\n")
					#print l[-1]
					#print '\n'
					#print 'done'
					#time.sleep(11)
				#else:
				#	print val
				#	time.sleep(1)
		out.close()		
		
		
'''	
residue_frequency = ResidueFrequency()
residue_frequency.file_read('ATPout')
'''



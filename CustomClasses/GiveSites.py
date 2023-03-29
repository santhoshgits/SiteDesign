import time
import os
from collections import Counter, defaultdict
import sys
import shutil
#from DockAnalyse import DockAnalyse

dire = os.getcwd()

if len(sys.argv) == 2:
	Folder = sys.argv[1]
else:
	print('GiveSites.py <Folder>')	
	sys.exit()

aline = open(dire+'/'+Folder+'/final2.pdb', 'r').readlines()

dic = defaultdict(list)
for line in aline:
	if line[:4] == "ATOM":
		dic[line[17:26]].append(line)

dire = os.getcwd()
path_out = dire+'/'+Folder+'/site'
if not os.path.exists(path_out):
	os.mkdir(path_out)
else:
	shutil.rmtree(path_out)
	os.mkdir(path_out)

'''
path_out1 = dire+'/'+Folder+'/site1'
if not os.path.exists(path_out1):
	os.mkdir(path_out1)
else:
	shutil.rmtree(path_out1)
	os.mkdir(path_out1)	
dock_analyse = DockAnalyse()	
'''

aline = open(dire+'/'+Folder+'/PocketDesignOutput.txt', 'r').readlines()	
#aline = aline[:-1]
count = 0
count1 = 1
for line in aline:
	count += 1
	#print count
	line = line.strip()
	l = line.split('\t')
	out = open(dire+'/'+Folder+'/site/'+l[0]+'-'+l[1]+'-'+str(count)+'.pdb', 'w')
	#print str(count1)+'-'+l[1]+'-'+str(count)+'.pdb'
	#time.sleep(.5)
	for i in l[2].split('_'):
		out.write(''.join(dic[i]))
	out.close()	
	if count == 15:
		count = 0
		count1 += 1
	#time.sleep(11)


'''

dock_analyse.file_read(Folder)

aline = open(dire+'/'+Folder+'/clustered_site', 'r').readlines()
count = 0
for line in aline:
	line = line.strip()
	l = line.split(' ')
	#print l
	shutil.copy(dire+'/'+Folder+'/site/'+l[0], dire+'/'+Folder+'/site1')
'''







import time
import os
import sys
import shutil
from tqdm import tqdm


dire = os.getcwd()
direExact = os.getcwd()

if len(sys.argv) == 5:
	Fold = sys.argv[1]
	nprocs = sys.argv[2]
	outFold = sys.argv[3]
	lig = sys.argv[4]
else:
	print('NRScan.py <LigFoldId> <nprocs> <OutFold> <Lig.pdb>')
	sys.exit()
	


dire = dire+'/'+Fold	

dire1 = direExact+'/CustomDatas/nrsiteActual'
	
path_out_fold = dire+'/SitesNr'
if not os.path.exists(path_out_fold):
	os.mkdir(path_out_fold)
else:
	shutil.rmtree(path_out_fold)
	os.mkdir(path_out_fold)	
	
	
print('Copying files....')
for i in tqdm(os.listdir(dire1)):
	if '_'+lig+'_' in i:
		shutil.copy(dire1+'/'+i, path_out_fold)
	

os.system('python '+direExact+'/CustomClasses/CreateVectors.py SitesNr SitesNrVec '+outFold)
time.sleep(1)


for i in os.listdir(dire+'/SiteVec'):
	shutil.copy(dire+'/SiteVec/'+i, dire+'/SitesNrVec')
	
	
out = open(dire+'/PairsSelf.txt','w')
for i in os.listdir(dire+'/SitesNrVec'):
	if '-' in i:
		for j in os.listdir(dire+'/SitesNrVec'):	
			if '-' not in j:
				out.write(i+'\t'+j+'\n')
out.close()


os.system('python '+direExact+'/CustomClasses/PocketApproxAlignParallel.py SitesNrVec PairsSelf.txt SiteSimSelf.txt '+str(nprocs)+' '+outFold)


shutil.rmtree(path_out_fold)
shutil.rmtree(dire+'/SitesNrVec')




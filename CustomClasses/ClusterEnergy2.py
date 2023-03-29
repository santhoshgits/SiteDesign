import numpy as np
import random
from collections import Counter, defaultdict
import time
import math
import os
from Dock import Dock
import shutil
#from pathos.multiprocessing import ProcessingPool as Pool
#from pathos.multiprocessing import Pool
from multiprocessing import Pool, freeze_support
#from pathos.pools import ProcessPool as Pool
#from multiprocessing.pool import ThreadPool as Pool
import sys
from tqdm import tqdm


dock = Dock()
dire = os.getcwd()

# # export PATH=/media/santhosh/7add6b54-492d-474d-b478-871833f96d09/home/santhosh/Downloads/mgltools_x86_64Linux2_1.5.6/bin:$PATH

#print dire
#time.sleep(1)

if len(sys.argv) == 4:
	lig_pdb = sys.argv[1]
	nproc = int(sys.argv[2])
	Folder = sys.argv[3]
else:
	print ("ClusterEnergy2.py <pdb> <nprocs> <OutFolder>")
	sys.exit()
#print lig_pdb

lig_aline = open(dire+'/'+Folder+'/'+lig_pdb, 'r').readlines()


nproc = int(nproc/2)
if nproc <= 1:
	nproc = 1
	
#print (nproc,' no of procs')


def dock_pool(Total):
	val, nos = Total
	#print nos
	try:
		if os.path.isfile(dire+'/'+Folder+"/temp/"+str(nos)+".pdbqt"):
			os.system('rm '+dire+'/'+Folder+"/temp/"+str(nos)+".pdbqt")
		if os.path.isfile(dire+'/'+Folder+"/temp/"+str(nos)+".pdb"):
			os.system('rm '+dire+'/'+Folder+"/temp/"+str(nos)+".pdb")	
		return dock.file_read(val, nos, Folder)
		#return 10
	except:
		return '0.0'



if __name__ == '__main__':
	def file_read(nprocs, lig_aline):
		#print('HI')
		'''
		out = open(dire+'/'+Folder+'/lig.pdb', 'w')
		for i in lig_aline:
			out.write(i)
		out.close()
		'''
		
		shutil.copy(dire+'/CustomClasses/prepare_ligand4.py', dire+'/'+Folder)
		os.chdir(dire+'/'+Folder)
		os.system('python2.5 prepare_ligand4.py -l lig.pdb -o lig.pdbqt')
		os.chdir(dire)
		
		if os.path.isfile(dire+'/'+Folder+"/temp/lig.pdbqt"):
			os.system('rm '+dire+'/'+Folder+"/temp/lig.pdbqt")
		if os.path.isfile(dire+'/'+Folder+"/temp/AD_score.py"):
			os.system('rm '+dire+'/'+Folder+"/temp/AD_score.py")
		if os.path.isfile(dire+'/'+Folder+"/temp/prepare_receptor4.py"):
			os.system('rm '+dire+'/'+Folder+"/temp/prepare_receptor4.py")		
		
		shutil.move(dire+'/'+Folder+"/lig.pdbqt", dire+'/'+Folder+"/temp/lig.pdbqt")
		shutil.copy(dire+"/CustomClasses/AD_score.py", dire+'/'+Folder+"/temp")
		shutil.copy(dire+"/CustomClasses/prepare_receptor4.py", dire+'/'+Folder+"/temp")
		
		res_len_dict = {'GLY':4 ,'ALA':5 ,'LEU':8 ,'ILE':8 ,'TRP':14 ,'SER':6 ,'THR':7 ,'TYR':12 ,'PHE':11 ,'PRO':7 ,'ASP':8 ,\
		   'GLU':9 ,'HIS':10 ,'CYS':6 ,'MET':8 ,'VAL':7 ,'ASN':8 ,'GLN':9 ,'ARG':11 ,'LYS':9 }
		dic = defaultdict(list)
		aline = open(dire+'/'+Folder+"/final_residue_pool.pdb", 'r').readlines()
		dic_ca = defaultdict(list)
		dic_cn = defaultdict(list)
		for line in aline:
			dic[line[17:26]].append(line)
			if line[13:16] == "CA ":
				dic_ca[line[17:26]].append([  float(line[27:38].strip()), float(line[38:46].strip()), float(line[46:54].strip()) ])
			dic_cn[line[17:26]].append([  float(line[27:38].strip()), float(line[38:46].strip()), float(line[46:54].strip()) ])	
		
		aline = open(dire+'/'+Folder+"/residue_cluster.txt", 'r').readlines()
		residue_cluster = []
		for line in aline:
			line = line.strip()
			residue_cluster.append(line.split("-"))

		#residue_cluster = [ i for i in residue_cluster if i[0] == 'TRP A3299' ]
		residue_cluster = sorted(residue_cluster, key = lambda x:len(x), reverse=True)	

		residue_map_dic = {}
		aline = open(dire+'/'+Folder+'/reside_rename_map.txt', 'r').readlines()
		for line in aline:
			line = line.strip()
			l = line.split('\t')
			residue_map_dic[l[1]] = l[2]

		check = True
		cluster_dic = {}
		residue_picked = []
		count = 0	
		

		for i in tqdm(residue_cluster):
			#break
			if i[0] not in cluster_dic:
				count += 1
				#if i[0] == "ASP A2967":
				#print(i[0],'---')
				#time.sleep(1)
				cluster_dic[i[0]] = 0
				i1 = random.choice(i)
				#print i1
				x, y, z = dic_ca[i1][0]
				x2, y2, z2 = np.asarray(dic_cn[i1], dtype='float').mean(axis=0)
				residues = []
				residues.append(i)
				for j in residue_cluster:
					if j[0] not in cluster_dic:
						#print j[0]				
						j1 = random.choice(j)
						x1, y1, z1 = dic_ca[j1][0]
						x3, y3, z3 = np.asarray(dic_cn[j1], dtype='float').mean(axis=0)
						x_ans = pow(( x - x1 ),2)
						y_ans = pow(( y - y1 ),2)
						z_ans = pow(( z - z1 ),2)
						ans_cn = math.sqrt(x_ans + y_ans + z_ans)
						if ans_cn < 1:
							x_ans = pow(( x2 - x3 ),2)
							y_ans = pow(( y2 - y3 ),2)
							z_ans = pow(( z2 - z3 ),2)
							ans_cn = math.sqrt(x_ans + y_ans + z_ans)
							if ans_cn < 2:
								cluster_dic[j[0]] = 0
								residues.append(j)
				ans = []
				val_arr = []
				Total = []
				kc = 0
				#print residues
				for k in residues:
					kc += 1
					k1 = random.choice(k)
					val = "".join(dic[k1])
					val_arr.append(val)
					Total.append([val, kc])
				#print len(val_arr)
				#print range(len(val_arr))
				p = Pool(nprocs)	
				#print Total[0]
				start = time.time()
				#freeze_support()
				p_ans = p.map(dock_pool, Total)
				time.sleep(.1)
				end = time.time()
				#print '\n\n' ps -o nlwp pid
				#print end-start
				#print p_ans
				#print '\n\n'
				#p.join()
				p.close()
				#time.sleep(1)
				for k2 in p_ans:
					#print k2
					ans.append(k2)
				residue_result = []	
				#print '\n\n'
				for l in range(0, len(ans)):
					#print ans[l], residues[l], residues[l][0]
					if float(ans[l]) < -0.05:
						map_freq = 0
						for i in residues[l]:
							map_freq += int(residue_map_dic[i])
						#print '\n'	
						if residues[l][0][:3] in res_len_dict:
							residue_result.append(residues[l][0]+"\t"+str(ans[l])+"\t"+\
												  str((float(ans[l])/res_len_dict[residues[l][0][:3]])*(res_len_dict[residues[l][0][:3]]*(len(residues[l])+map_freq)))+
												  "\t"+str(float(ans[l])/res_len_dict[residues[l][0][:3]])+"\t"+str(len(residues[l])))

				residue_result = sorted(residue_result, key = lambda x:float(x.split("\t")[2]))
				#print ans[:3], len(ans),'--'
				#print (residue_result)
				#print '\n'
				#time.sleep(11)
				if residue_result:
					cutoff = np.asarray([i.split('\t')[2] for i in residue_result], dtype='float').mean()/2
					#cutoff1 = np.asarray([i.split('\t')[2] for i in residue_result], dtype='float').mean()/1.5
					#cutoff2 = np.asarray([i.split('\t')[2] for i in residue_result], dtype='float').mean()/2
					#print cutoff
					for m in residue_result:
						m1 = m.split("\t")
						#print m1[2]
						#if float(m1[2]) < float(cutoff):
							#print m1,' asa'
						residue_picked.append(m1[0])
				#print 'domne'
				#time.sleep(11)
			#break	
						
		#print residue_picked,'picked'
		dic = {}
		for i in residue_picked:
			dic[i] = 0
		aline = open(dire+'/'+Folder+"/residue_cluster.txt", 'r').readlines()	
		out = open(dire+'/'+Folder+'/residue_cluster1.txt', 'w')	
		for line in aline:
			line = line.strip()
			l = line.split('-')
			if l[0] in dic:
				out.write(line+"\n")
		out.close()			


file_read(nproc, lig_aline)









	
	
	
	
	
	
	
	

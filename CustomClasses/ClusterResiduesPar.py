import os
import numpy as np
import math
from collections import Counter, defaultdict
import random
import time
from multiprocessing import Pool
import sys
import copy

if len(sys.argv) == 4:
	final = sys.argv[1]
	nprocs = int(sys.argv[2])
	Folder = sys.argv[3]
else:
	print ('ClusterResiduesPar.py <final-residue-pool.pdb> nprocs')
	sys.exit()

global res_whole_dic, res_coord_dic

dire = os.getcwd()

def neighbours(Total):
	res_taken, res_list = Total
	res_coord = np.asarray(res_coord_dic[res_taken], dtype='float')
	x_ca, y_ca, z_ca = res_coord[1]
	x_cn, y_cn, z_cn = res_coord.mean(axis=0)
	#print res_taken," --"
	res_final = []
	res_final.append(res_taken)
	for i in res_list:
		if i != res_taken:
			coord = np.asarray(res_coord_dic[i], dtype='float')
			x_ca1, y_ca1, z_ca1 = coord[1]
			x_ans = pow(( x_ca1 - x_ca ),2)
			y_ans = pow(( y_ca1 - y_ca ),2)
			z_ans = pow(( z_ca1 - z_ca ),2)
			ans_ca = math.sqrt(x_ans + y_ans + z_ans)
			# if ans_ca < .5:
			if ans_ca < .75:
				x_cn1, y_cn1, z_cn1 = coord.mean(axis=0)
				x_ans = pow(( x_cn1 - x_cn ),2)
				y_ans = pow(( y_cn1 - y_cn ),2)
				z_ans = pow(( z_cn1 - z_cn ),2)
				ans_cn = math.sqrt(x_ans + y_ans + z_ans)
				#if ans_cn < 1:
				if ans_cn < 1.5:
					res_final.append(i)
	return res_final


def check_clash(res1, res2):
	#print res1, res2
	res_coord = np.asarray(res_coord_dic[res1], dtype='float')
	x_ca, y_ca, z_ca = res_coord[1]
	x_cn, y_cn, z_cn = res_coord.mean(axis=0)
	
	res_coord1 = np.asarray(res_coord_dic[res2], dtype='float')
	x_ca1, y_ca1, z_ca1 = res_coord1[1]
	x_cn1, y_cn1, z_cn1 = res_coord1.mean(axis=0)
	
	x_ans = pow(( x_ca1 - x_ca ),2)
	y_ans = pow(( y_ca1 - y_ca ),2)
	z_ans = pow(( z_ca1 - z_ca ),2)
	ans_ca = math.sqrt(x_ans + y_ans + z_ans)
	
	if ans_ca > .7:
		return False

	x_ans = pow(( x_cn1 - x_cn ),2)
	y_ans = pow(( y_cn1 - y_cn ),2)
	z_ans = pow(( z_cn1 - z_cn ),2)
	ans_cn = math.sqrt(x_ans + y_ans + z_ans)
	if ans_cn > 1.5:
		return False
	return True


def ResMerge(ResInp):
	res_single = []
	ResInp.append(ResInp[0])
	for i in ResInp:
		res_single.append(i[0])

	dic = {}
	merged_unit = []
	for i in range(len(res_single)):
		if i not in dic:
			for j in range(len(res_single)):
				if j not in dic:
					if i != j:
						#print i,j
						case = check_clash(res_single[i], res_single[j])
						if case:
							dic[j] = 0
							merged_unit.append([i,j])
							
	dic = {}				
	ResInp1 = []
	for i in merged_unit:
		arr = []
		for j in ResInp[i[0]]:
			arr.append(j)
		for j in ResInp[i[1]]:
			arr.append(j)
		arr = sorted(set(arr))	
		ResInp1.append(arr)
		dic[i[0]] = 0
		dic[i[1]] = 0
		
	for i in range(len(ResInp)):
		if i not in dic:
			ResInp1.append(ResInp[i])
	return ResInp1
		

def file_read(aline):
	print (aline)
	aline = open(dire+'/'+Folder+'/'+aline, 'r').readlines()
	global res_whole_dic, res_coord_dic
	res_whole_dic, res_coord_dic = defaultdict(list), defaultdict(list)
	res_list = []
	for line in aline:
		line = line.strip()
		#print line
		res_whole_dic[line[17:26]].append(line)
		res_coord_dic[line[17:26]].append([  float(line[27:38].strip()), float(line[38:46].strip()), float(line[46:54].strip()) ])
		res_list.append(line[17:26])
	res_list = sorted(set(res_list))
	
	check = True
	res_total = []
	p = Pool(nprocs)
	cc = 0
	while check:
		cc += 1
		#if cc > 20:
		#	check = False
		print (len(res_list), end='\r')
		random.shuffle(res_list)
		res_taken = [ random.choice(res_list) for _ in range(nprocs)]
		#res_taken = ['PHE A 297']
		
		res_lists = [ res_list for _ in range(nprocs)]
		res_fuse = []
		for i in range(0, len(res_taken)):
			res_fuse.append([res_taken[i], res_lists[i]])
		res_out = p.map(neighbours, res_fuse)
		#print cc
		
		#print res_out
		#time.sleep(11)
		res_final_total = ResMerge(res_out)
		#print res_final_total
		#time.sleep(11)
		for res_final in res_final_total:
			dic = {}
			res_total.append(res_final)
			#print len(res_list)
			for i in res_final:
				dic[i] = 0
			res_list_new = []	
			for i in res_list:
				if i not in dic:
					res_list_new.append(i)
			res_list = copy.deepcopy(res_list_new)	
			if len(res_list_new) < 10:
				check = False
			
	
	out = open(dire+'/'+Folder+"/residue_cluster.txt", 'w')
	for i in res_total:
		out.write("-".join(i)+"\n")
	out.close()
	


file_read(final)






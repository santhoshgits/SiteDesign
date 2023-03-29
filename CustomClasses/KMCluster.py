import time
import math
import numpy as np
from sklearn.cluster import KMeans
from sklearn.preprocessing import MinMaxScaler
import sys
from collections import Counter, defaultdict
import os


if len(sys.argv) == 2:
	Folder = sys.argv[1]
else:
	print('python3.8 KMCluster.py <Folder>')	
	sys.exit()

dire = os.getcwd()

aline = open(dire+'/'+Folder+'/final2.pdb', 'r').readlines()

'''
grids = []
for line in aline:
    line = line.strip()
    grid = list(map(float, [line[28:38].strip(), line[38:46].strip(), line[46:54].strip() ] ))
    grid = list(map(int, grid))
    grids.append(grid)

grids = np.asarray(grids, dtype='int')

prev = KMeans(n_clusters=1).fit(grids).inertia_
for i in range(2,20):
    kmeans = KMeans(n_clusters=i).fit(grids)
    print kmeans.inertia_, kmeans.inertia_/float(prev)
    prev = kmeans.inertia_
'''

final_dic = defaultdict(list)
aline = open(dire+'/'+Folder+'/final2.pdb', 'r').readlines()
for j in aline:
    j1 = j.strip()
    final_dic[j1[17:26]].append([j1[28:38].strip(), j1[38:46].strip(), j1[46:54].strip()])





def ResidueNearest(het_coord):
    x, y, z = list(map(float,het_coord))
    aline = open(dire+'/'+Folder+'/final2.pdb', 'r').readlines()
    arr = []
    for j, j2 in final_dic.items():
        for j1 in j2:
            x1, y1, z1 = list(map(float, j1))
            x_ans = pow(( x - x1),2)
            y_ans = pow(( y - y1),2)
            z_ans = pow(( z - z1),2)
            ans = math.sqrt(x_ans + y_ans + z_ans)
            if ans < 3:
                arr.append(j)
    arr = sorted(set(arr))
    arr1 = []
    for i in arr:
        a,b,c = np.asarray(final_dic[i], dtype='float').mean(axis=0)
        #time.sleep(1)
        x, y, z = list(map(int,map(float,final_dic[i][1])))
        #arr1.append(list(map(int,map(float,final_dic[i][1]))))
        arr1.append([x, y, z, a, b, c])
    return arr1



def KmeanNearest(KmeanCoord):
    ans1 = []
    for i in KmeanCoord:
        #print i
        x, y, z, a, b, c = i
        arr = []
        for j, j2 in final_dic.items():
            for j1 in j2:
                x1, y1, z1 = list(map(float, j1))
                x_ans = pow(( x - x1),2)
                y_ans = pow(( y - y1),2)
                z_ans = pow(( z - z1),2)
                ans = math.sqrt(x_ans + y_ans + z_ans)
                if ans < 2:
                    arr.append(j)
        arr = sorted(set(arr))
        arr2 = []
        for j in arr:
            #print final_dic[j]
            #print np.asarray(final_dic[j], dtype='float')
            a1, b1, c1 = np.asarray(final_dic[j], dtype='float').mean(axis=0)
            x_ans = pow(( a - a1),2)
            y_ans = pow(( b - b1),2)
            z_ans = pow(( c - c1),2)
            ans = math.sqrt(x_ans + y_ans + z_ans)
            if ans < 3:
                arr2.append(j)
        ans1.append(arr2)
    lnsa = [ len(i) for i in ans1 ]
    lns = np.asarray([ len(i) for i in ans1 ]).reshape(-1,1)
    lns = np.asarray(lns, dtype='float')
    lns_transformed = [ i*100+10 for i in MinMaxScaler().fit_transform(lns) ]
    lns_transformed = np.asarray(lns_transformed, dtype='float')
    lns_transformed = np.asarray(lns_transformed, dtype='int')
    #lns_transformed = list(map(str, lns_transformed))
    arr2 = []
    for i in range(len(lns_transformed)):
        arr2.append('-'.join(ans1[i])+'\t'+str(lns_transformed[i][0])+'\t'+str(lnsa[i]) )
    #print 'done'
    #print [ len(i) for i in ans1 ]
    #time.sleep(11)
    return arr2




def KmeanNearest_Patch(KmeanCoord):
    ans1 = []
    for i in KmeanCoord:
        #print i
        x, y, z, a, b, c = i
        arr = []
        for j, j2 in final_dic.items():
            for j1 in j2:
                x1, y1, z1 = list(map(float, j1))
                x_ans = pow(( x - x1),2)
                y_ans = pow(( y - y1),2)
                z_ans = pow(( z - z1),2)
                ans = math.sqrt(x_ans + y_ans + z_ans)
                if ans < .75:
                    #print ans
                    arr.append(j)
        arr = sorted(set(arr))
        arr2 = []
        for j in arr:
            #print final_dic[j]
            #print np.asarray(final_dic[j], dtype='float')
            a1, b1, c1 = np.asarray(final_dic[j], dtype='float').mean(axis=0)
            x_ans = pow(( a - a1),2)
            y_ans = pow(( b - b1),2)
            z_ans = pow(( c - c1),2)
            ans = math.sqrt(x_ans + y_ans + z_ans)
            if ans < 3:
                arr2.append(j)
        ans1.append(arr2)
    lnsa = [ len(i) for i in ans1 ]
    lns = np.asarray([ len(i) for i in ans1 ]).reshape(-1,1)
    lns = np.asarray(lns, dtype='float')
    lns_transformed = [ i*100+10 for i in MinMaxScaler().fit_transform(lns) ]
    lns_transformed = np.asarray(lns_transformed, dtype='float')
    lns_transformed = np.asarray(lns_transformed, dtype='int')
    #lns_transformed = list(map(str, lns_transformed))
    arr2 = []
    for i in range(len(lns_transformed)):
        arr2.append('-'.join(ans1[i])+'\t'+str(lns_transformed[i][0])+'\t'+str(lnsa[i]) )
    #print 'done'
    #print [ len(i) for i in ans1 ]
    #time.sleep(11)
    return arr2




def patch(arr):
	res_arr = arr.split('\t')[0].split('-')
	# final_dic
	coord = []
	for i in res_arr:
		x, y, z = final_dic[i][1]
		#print  np.asarray(final_dic[i], dtype='float').mean(axis=0)
		a, b, c = np.asarray(final_dic[i], dtype='float').mean(axis=0)
		#print a,b,c
		#coord.append(final_dic[i][1])
		coord.append([x, y, z, a, b, c])
		#print i, final_dic[i][1]
	coord = np.array(coord, dtype='float')
	#print len(coord)
	prev = KMeans(n_clusters=1).fit(coord).inertia_
	for i in range(2,200):
		try:
			kmeans = KMeans(n_clusters=i).fit(coord)
			#print kmeans.inertia_, kmeans.inertia_/float(prev), prev-kmeans.inertia_
			if prev-kmeans.inertia_ < 10:
			    break
			if kmeans.inertia_/float(prev) > .85:
			    break
			prev = kmeans.inertia_
		except:
			pass
			break	
	#print i-1
	try:
		res_pool = KmeanNearest_Patch(kmeans.cluster_centers_)
		#for i in res_pool:
		#	print len(i.split('\t')[0].split('-'))
		#print len(res_pool)
		return res_pool
	except:
		return None	




aline = open(dire+'/'+Folder+'/lig.pdb', 'r').readlines()
het_line = []
for line in aline:
	line = line.strip()
	if line[:6] == 'HETATM':
		if len(line) > 2:
		    het_line.append( [line[7:11], [float(line[28:38].strip()), float(line[38:46].strip()), float(line[46:54].strip())] ] )

out = open(dire+'/'+Folder+'/clusters.txt', 'w')
for het in het_line:
    #print het[0]
	grids = ResidueNearest(het[1])
	if len(grids) > 0:
		#print het, grids
		prev = KMeans(n_clusters=1).fit(grids).inertia_
		for i in range(2,200):
			try:
				kmeans = KMeans(n_clusters=i).fit(grids)
				#print kmeans.inertia_, kmeans.inertia_/float(prev), prev-kmeans.inertia_
				if prev-kmeans.inertia_ < 10:
				    break
				if kmeans.inertia_/float(prev) > .85:
				    break
				prev = kmeans.inertia_
			except:
				pass
				break	
		kmeans = KMeans(n_clusters=i-1).fit(grids)
		#print len(kmeans.cluster_centers_)
		res_pool = KmeanNearest(kmeans.cluster_centers_)
		for i in res_pool:
		    #print i
			res_pool_patch = patch(i)
			if res_pool_patch:
				for j in res_pool_patch:
					out.write(j+'\n')
		#print len(res_pool)
		#print '\n'
		#time.sleep(11)
out.close()









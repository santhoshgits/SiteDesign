import os
import subprocess
import time
import copy


class Dock():
	
	def __init__(self):
		#self.cwd = os.getcwd()
		#self.temp_dir = self.cwd+"/temp"
		#print self.cwd
		#print self.temp_dir
		pass
	
	def file_read(self, aline, rno, Folder):
		cwd = os.getcwd()
		#print cwd,'--cwd'
		temp_dir = copy.deepcopy(cwd)+'/'+Folder+"/temp"
		p2 = str(1000.0)
		#print self.cwd
		#print self.temp_dir
		os.chdir(temp_dir)
		
		if len(aline.split("\n")) < 1:
			#print "hii"
			return p2
		rno = str(rno)
		out = open(rno+".pdb", 'w')
		#print rno,"--------------------------------"
		out.write(aline)
		out.close()
		#time.sleep(.01)
		#print aline[:30],rno
		#print compute_truth
		#os.system("python prepare_receptor4.py -r "+rno+".pdb -o "+rno+".pdbqt")
		#if compute_truth:
		#print os.getcwd()
		
		'''
		if not os.path.isfile(rno+".pdb"):
			print aline.split('\n')[:3],rno
			print '\n'	
			print os.getcwd()
			#print self.temp_dir
			time.sleep(11)
		'''	
		
		try:
			#os.system("python prepare_receptor4.py -r "+rno+".pdb -o "+rno+".pdbqt")
			c1 = subprocess.Popen("python2.5 prepare_receptor4.py -r "+rno+".pdb -o "+rno+".pdbqt", shell=True, stdout=subprocess.PIPE)
			c2 = c1.communicate()[0]

			p1 = subprocess.Popen("python2.5 AD_score.py "+rno+".pdbqt lig.pdbqt", shell=True, stdout=subprocess.PIPE)
			p2 = p1.communicate()[0]
		except:
			#sys.exit()
			pass
		p2 = p2.strip()	
		p2 = float(p2)
		p2 = str(p2)
		#p2 = "1000"
		#print os.getcwd(),' -- 1'
		os.chdir(cwd)
		#print os.getcwd(),' --- 2'
		#print(p2)
		return p2
	
'''	
dock = Dock()
aline = open("obj01.pdb", 'r').readlines()
aline = "".join(aline)
print aline
print dock.file_read(aline, 25)
'''





	

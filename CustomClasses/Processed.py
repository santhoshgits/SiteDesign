import time
import os
#from SiteGen15a import SiteGen
from CustomClasses.SiteGen16 import SiteGen

dire = os.getcwd()

class Processed():
	
	def file_read(self, Total):
		x, ligand_aline, minim, maxim, Folder = Total
		res_point_line = open(dire+'/'+Folder+'/residue_position.txt', 'r').readlines()
		#print (ligand_aline,'-----')
		site_gen = SiteGen(res_point_line, ligand_aline, Folder)
		#site_gen_fit = SiteGenFit()
		return site_gen.file_read(1, ligand_aline, minim, maxim)
	

	

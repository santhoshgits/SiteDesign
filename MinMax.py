import time
import os
import random
import numpy as np
from pickle import dump, load
import sys
import subprocess
from sklearn.preprocessing import LabelBinarizer
from sklearn.preprocessing import MinMaxScaler

from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense

from collections import Counter

from PyBioMed.PyMolecule import connectivity, topology
from PyBioMed import Pymolecule
from PyBioMed.PyMolecule.cats2d import CATS2D

from rdkit import Chem
import random

import tensorflow_probability as tfp
import tensorflow as tf
from tqdm import tqdm

from openbabel import openbabel

tfd = tfp.distributions
tfk = tf.keras

dire = os.getcwd()

if len(sys.argv) == 2:
	pdb = sys.argv[1]
else:
	print('\n\n')
	print("INSTRUCTION TO EXECUTE THE CODE")
	print('For Example: To run MinMax.py for ATP ligand, the syntax would be like the below.')
	print('python3.8 MinMax.py ATP.pdb')
	print('\n')
	print('Please do not supply additional arguments. Thank You')
	print('\n')
	sys.exit()

print('Running Open Babel for the Ligand \n')

openbabel.obErrorLog.StopLogging()
obConversion = openbabel.OBConversion()
obConversion.SetInAndOutFormats("pdb", "smi")

mol = openbabel.OBMol()
rnd = '1'
obConversion.ReadFile(mol, pdb)
obConversion.WriteFile(mol, rnd+'.smi')
		
print(rnd+'.smi')		
	
		
#os.system('obabel -ipdb '+pdb+' -osmi > a.smi')
aline = open(rnd+'.smi', 'r').readlines()
print(aline)
os.system('rm '+rnd+'.smi')
mol = aline[0].split('\t')[0]


mol = Chem.MolFromSmiles(mol, sanitize=False)
mol.UpdatePropertyCache(strict=False)
mol_py = Pymolecule.PyMolecule()
mol_py.ReadMolFromSmile(aline[0].split('\t')[0])
EstateDesc = list(mol_py.GetEstate().values())
ChargeDesc = list(mol_py.GetCharge().values())
CATDesc = list(CATS2D(mol,PathLength = 10,scale = 3).values())
ConDesc = list(connectivity.GetConnectivity(mol).values())
TopDesc = list(topology.GetTopology(mol).values())


#print(len(EstateDesc), len(ChargeDesc), len(CATDesc), len(ConDesc), len(TopDesc))
#print(len(EstateDesc)+len(ChargeDesc)+len(CATDesc)+len(ConDesc)+len(TopDesc))

EstateDesc = list(map(str, EstateDesc))
ChargeDesc = list(map(str, ChargeDesc))
CATDesc = list(map(str, CATDesc))
ConDesc = list(map(str, ConDesc))
TopDesc = list(map(str, TopDesc))

print(len(EstateDesc), len(ChargeDesc), len(CATDesc), len(ConDesc), len(TopDesc)) #237 25 150 44 25



LigandDescriptors = '\t'.join(EstateDesc)+'\t'+'\t'.join(ChargeDesc)+'\t'+'\t'.join(CATDesc)+'\t'+'\t'.join(ConDesc)+'\t'+'\t'.join(TopDesc)
LigandDescriptors = LigandDescriptors.split('\t')

LigandDescriptors = np.asarray(LigandDescriptors, dtype='float')
LigandDescriptors = np.reshape(LigandDescriptors,(1,-1))
LigandDescriptors = np.asarray(LigandDescriptors, dtype='float')
#print(LigandDescriptors.shape)


# Below is a probalistic model 

if not os.path.exists(dire+'/CustomDatas/Descriptors'):
	os.chdir(dire+'/CustomDatas')
	c1 = subprocess.Popen('conda activate SiteDesign ; python LigDesc.py', shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	c2 = c1.communicate()[0]
	os.chdir(dire)
	os.system('python MinMax.py '+pdb)	



aline = open(dire+'/CustomDatas/Descriptors', 'r').readlines()
x, y, lig = [], [], []
for line in aline:
	line = line.strip()
	l = line.split('\t')
	if int(l[-1]) > 4:
		x.append(l[1:-1])
		y.append(l[-1])
		lig.append(l[0])
	
x1, y1 = [], []		
c = 0
for i in x:
	check = True 
	for j in range(len(i)):
		if i[j] == 'nan' or i[j] == 'inf' or i[j] == '-inf':
			i[j] = 0	
			check = False
	if check:
		for _ in range(1):
			x1.append(i)
			y1.append(y[c])
	c += 1	

x = x1				
y = y1	
	
#print(x)
if len(x) == 0:
	os.chdir(dire+'/CustomDatas')
	print('Vector shape mismatch is detected. We need to rerun LigDesc.py')
	print('This is one time operation. please wait 6 hours')
	c1 = subprocess.Popen('conda activate SiteDesign ; python LigDesc.py', shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	c2 = c1.communicate()[0]
	os.chdir(dire)
	os.system('python MinMax.py '+pdb)	
	
	
x = np.array(x, dtype='float')	
y = np.array(y, dtype='float')	


scaler_x = MinMaxScaler()
scaler_x.fit(x)	
x = scaler_x.transform(x)
print(LigandDescriptors.shape[1], x.shape[1])
if LigandDescriptors.shape[1] != x.shape[1]:
	os.chdir(dire+'/CustomDatas')
	c1 = subprocess.Popen('conda activate SiteDesign ; python LigDesc.py', shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	c2 = c1.communicate()[0]
	os.chdir(dire)
	os.system('python MinMax.py '+pdb)


LEN = len(x[0])
print('\n\n')
print('Total number of features is : ',LigandDescriptors.shape[1])
print(x.shape)
print('\n\n')


# basic keras regression
input_shape = (LEN,) 




def posterior_mean_field(kernel_size, bias_size=0, dtype=None):
  n = kernel_size + bias_size
  c = np.log(np.expm1(1.))
  return tf.keras.Sequential([
      tfp.layers.VariableLayer(2 * n, dtype=dtype),
      tfp.layers.DistributionLambda(lambda t: tfd.Independent(  # pylint: disable=g-long-lambda
          tfd.Normal(loc=t[..., :n],
                     scale=1e-5 + tf.nn.softplus(c + t[..., n:])),
          reinterpreted_batch_ndims=1)),
  ])
  
  

def prior_trainable(kernel_size, bias_size=0, dtype=None):
  n = kernel_size + bias_size
  return tf.keras.Sequential([
      tfp.layers.VariableLayer(n, dtype=dtype),
      tfp.layers.DistributionLambda(
          lambda t: tfd.Independent(tfd.Normal(loc=t, scale=1),  # pylint: disable=g-long-lambda
                                    reinterpreted_batch_ndims=1)),
  ])

  
# Aleatoric undertainity CI range

negloglik = lambda y, p_y: -p_y.log_prob(y)


model = tfk.Sequential([
  tf.keras.layers.Dense(1 + 1),
  tfp.layers.DistributionLambda(
      lambda t: tfd.Normal(loc=t[..., :1],
                           scale=1e-3 + tf.math.softplus(0.05 * t[..., 1:]))),
]) # this performs better than the most


# Do inference.
model.compile(optimizer=tf.optimizers.Adam(learning_rate=0.01), loss=negloglik)
print('\n\n')
print('Tensoflow Model is Running. Please Wait')
print('\n\n')
model.fit(x, y, epochs=15, verbose=False)

LigandDescriptors = scaler_x.transform(LigandDescriptors)

LigandDescriptorsMulti = []
print('Making Predictions : --- ')
for i in range(30000):
	#print(LigandDescriptors.shape)
	LigandDescriptorsMulti.append(LigandDescriptors.flatten())
LigandDescriptorsMulti = np.asarray(LigandDescriptorsMulti, dtype='float')



MinMaxRange = [ i[0] for i in model.predict(LigandDescriptorsMulti) ]

#print('a')	

MinMaxRange = list(map(int, MinMaxRange))
#maxi = Counter(MinMaxRange).most_common(1)[0][1]/10.0
#MinMaxRange = [ i[0] for i in Counter(MinMaxRange).most_common(20) if i[1] > maxi ]
#print(MinMaxRange)

PredictedRange = []
while True:
	PredictedRange.append(random.choice(MinMaxRange))
	if len(PredictedRange) > 1000:
		break
	

'''
print(np.quantile(PredictedRange, 0.1))
print(np.quantile(PredictedRange, 0.25),'  .25')
print(np.quantile(PredictedRange, 0.5),' median')
print(np.quantile(PredictedRange, 0.75),'  .75')
print(np.quantile(PredictedRange, 0.9))
'''
print('\n\n')
print('The optimal Minimum and the Maximum no. of residues needed to design a site for your ligand is : {:0.0f} and {:0.0f} '.format(np.quantile(PredictedRange, 0.55), np.quantile(PredictedRange, 0.95)))
print('Please Note Down this Range. It will be needed later. Thank You!!!')
print('\n\n')








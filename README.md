# CRD (Cognate Receptor Discovery): a De novo Design algorithm for prediction of Cognate Protein Receptors for small molecule ligands.
## NOTE TO USER: After you clone this repo, please download the associated binding site files from the zenodo (https://zenodo.org/record/7780052), extract it and copy 'CustomDatas' folder to SiteDesign.

## Instruction to install CRD
The algorithm has been implemented in Python-3.9 (Anaconda distribution). It is therefore recommended to install the same version of python in your system to run our code.

Once Python-3.9 is installed, we then have to install many additional libraries such as Numba, Tensorflow, AutoDock and NumPy. To enable easy installation, A YAML file called environment.yml is provided which will invoke a seprate environments to run CRD safely.

****
**Use the below command to create environments**
```
 
1. conda env create -f environment.yml
The above command will create an environment called 'SiteDesign' and install all pre-requisite libararies in your machine.

2. conda activate SiteDesign
To activate and load all modules required to run CRD. 

3. conda deactivate.
To deactivate SiteDesign.
```


---


# Steps for Running CRD

The only input required by CRD is the PDB file of the ligand. As an example, I have already provided a sample ligand file (ATP.pdb) for the user to test.


## STEP-1: Implementing LigFrag, FragMatch, FragAlign, ResRes and ResFreq Functions.
### python SiteDesignInitiator.py Argument-1 Argument-2 Argument-3

Usage: python SiteDesignInitiator.py ATP.pdb 8 ATPOutput  

1. **Argument-1 -** Provide ligand coordinate for which sites needs to be designed.

2. **Argument-2 -** Specify the number of CPU cores to utilize for the run. 

3. **Argument-2 -** Specify the folder that will store the output.

## STEP-2: Implementing SiteSize Functions.
### python MinMax.py Argument-1

Usage: python MinMax.py ATP.pdb  

1. **Argument-1 -** Provide ligand coordinate for which sites needs to be designed.

<em> Please note down the two numbers reported by MinMax.py. We need to supply them as input for next step. </em>

## STEP-3: Running the SiteDesign module
This implements functions such as SiteGen, Triplet and use them to generate site based on the GA optimization.
### python SiteDesignGenerator.py Argument-1 Argument-2 Argument-3 Argument-4 Argument-5

Usage: python SiteDesignGenerator.py ATP.pdb 8 ATPOutput 18 25
  
1. **Argument-1 -** Provide ligand coordinate for which sites needs to be designed.

2. **Argument-2 -** Specify the number of CPU cores to utilize for the run. 

3. **Argument-3 -** Specify the folder that will store the output.

4. **Argument-4 -** Min number from MinMax.py

5. **Argument-5 -** Max number from MinMax.py

---
# After you complete all three steps, please refer to the folder 'site', this will be created inside the Argument-2 folder. The folder 'site' will contain 750 designed sites for your ligand.

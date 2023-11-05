import numpy
from rdkit.Chem.EState import EStateIndices
from rdkit.Chem.EState import AtomTypes
from rdkit import Chem
from rdkit.Chem import ChemicalFeatures
from rdkit.Chem.Pharm2D.SigFactory import SigFactory
from rdkit.Chem.Pharm2D import Generate

from rdkit import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit.Chem import MACCSkeys
from rdkit.Chem.AtomPairs import Torsions
from rdkit.Chem import AllChem
from rdkit.Chem.AtomPairs import Pairs
from rdkit import Chem
from rdkit.Chem.Pharm2D import Gobbi_Pharm2D,Generate
import pandas as pd
from fingerprint_list import *

data=['CCC','CO','CCCC']
dataa=pd.DataFrame(data,columns=['smile'])
# finger prints
FingerprintMol_fing(dataa['smile'])       #64
maccskey_fing(dataa['smile'])             #167
Morgan_fing(dataa['smile'])               #1024
Gobbi_pharma_2D(dataa['smile'])           #39972
Minimal_features(dataa['smile'])          #885
base_features(dataa['smile'])             #2988
Atomtype_fing(dataa['smile'])             #79
EStateIndices_fing(dataa['smile'])
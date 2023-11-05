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


def Minimal_features(smile):
    fdefName = 'MinimalFeatures.fdef'
    featFactory = ChemicalFeatures.BuildFeatureFactory(fdefName)
    sigFactory = SigFactory(featFactory,minPointCount=2,maxPointCount=3)  # min and max point (varie form 2,3 and 3,4)
    sigFactory.SetBins([(0,2),(2,5),(5,8)])#,(8,12),(12,17)])             #set bins for distance between two featres
    sigFactory.Init()
    sigFactory.GetSigSize()
    fingerprint_df=pd.DataFrame()
    for i in smile:
        mol = Chem.MolFromSmiles(i)
        fp = Generate.Gen2DFingerprint(mol,sigFactory)
        fp.GetNumOnBits()
        fing=fp.ToBitString()
        fing_df=pd.DataFrame(convert_bit_to_list(fing))
        fing_df=fing_df.T
        fingerprint_df=pd.concat([fingerprint_df,fing_df],axis=0)
    return fingerprint_df.reset_index(drop=True)

def base_features(smile):
    fdefName = 'BaseFeatures.fdef'
    featFactory = ChemicalFeatures.BuildFeatureFactory(fdefName)
    sigFactory = SigFactory(featFactory,minPointCount=2,maxPointCount=3)     # min and max point (varie form 2,3 and 3,4)
    sigFactory.SetBins([(0,2),(2,5),(5,8)])#,(8,12),(12,17)])                 #set bins for distance between two featres
    sigFactory.Init()
    sigFactory.GetSigSize()
    fingerprint_df=pd.DataFrame()
    for i in smile:
        mol = Chem.MolFromSmiles(i)
        fp = Generate.Gen2DFingerprint(mol,sigFactory)
        fp.GetNumOnBits()
        fing=fp.ToBitString()
        fing_df=pd.DataFrame(convert_bit_to_list(fing))
        fing_df=fing_df.T
        fingerprint_df=pd.concat([fingerprint_df,fing_df],axis=0)
    return fingerprint_df.reset_index(drop=True)

def Atomtype_FingerprintMol(mol):
    if AtomTypes.esPatterns is None:
      AtomTypes.BuildPatts()
    nPatts = len(AtomTypes.esPatterns)
    counts = numpy.zeros(nPatts, dtype=numpy.int64)
    sums = numpy.zeros(nPatts, dtype=numpy.float64)

    esIndices = EStateIndices(mol)
    for i, (_, pattern) in enumerate(AtomTypes.esPatterns):
        matches = mol.GetSubstructMatches(pattern, uniquify=1)
        counts[i] = len(matches)
        for match in matches:
            sums[i] += esIndices[match[0]]
    return counts.tolist()

def EStateIndices_FingerprintMol(mol):
    if AtomTypes.esPatterns is None:
      AtomTypes.BuildPatts()
    nPatts = len(AtomTypes.esPatterns)
    counts = numpy.zeros(nPatts, dtype=numpy.int64)
    sums = numpy.zeros(nPatts, dtype=numpy.float64)

    esIndices = EStateIndices(mol)
    for i, (_, pattern) in enumerate(AtomTypes.esPatterns):
        matches = mol.GetSubstructMatches(pattern, uniquify=1)
        counts[i] = len(matches)
        for match in matches:
            sums[i] += esIndices[match[0]]
    return sums.tolist()

def convert_bit_to_list(list):
    digits = str(list)
    split_digits = []
    for digit in digits:
        split_digits.append(int(digit))
    # print(split_digits)
    return split_digits

def FingerprintMol_fing(smile):
    fingerprint_df=pd.DataFrame()
    for i in smile:
        ms=Chem.MolFromSmiles(i)
        fing=FingerprintMols.FingerprintMol(ms).ToBitString()
        fing_df=pd.DataFrame(convert_bit_to_list(fing))
        fing_df=fing_df.T
        fingerprint_df=pd.concat([fingerprint_df,fing_df],axis=0)
    return fingerprint_df.reset_index(drop=True)


def  maccskey_fing(smile):
    fingerprint_df=pd.DataFrame()
    for i in smile:
        ms=Chem.MolFromSmiles(i)
        fing = MACCSkeys.GenMACCSKeys(ms) .ToBitString()
        fing_df=pd.DataFrame(convert_bit_to_list(fing))
        fing_df=fing_df.T
        fingerprint_df=pd.concat([fingerprint_df,fing_df],axis=0)
    return fingerprint_df.reset_index(drop=True)

# pharma list

def Morgan_fing (smile):
    fingerprint_df=pd.DataFrame()
    for i in smile:
        ms=Chem.MolFromSmiles(i)
        fing = AllChem.GetMorganFingerprintAsBitVect(ms,2,nBits=1024).ToBitString()   #(1 is to find all internal details(high quality screening) , 2 general details(general) , 3 is used for  global (ligands) )
        fing_df=pd.DataFrame(convert_bit_to_list(fing))
        fing_df=fing_df.T
        fingerprint_df=pd.concat([fingerprint_df,fing_df],axis=0)
    return fingerprint_df.reset_index(drop=True)
    # tts = [Torsions.GetTopologicalTorsionFingerprintAsIntVect(x) for x in ms]

def Gobbi_pharma_2D(smile):
    fingerprint_df=pd.DataFrame()
    for i in smile:
        ms=Chem.MolFromSmiles(i)
        fing = Generate.Gen2DFingerprint(ms,Gobbi_Pharm2D.factory).ToBitString()
        # fing = AllChem.GetMorganFingerprint(ms,2).ToBitString()   #(1 is to find all internal details(high quality screening) , 2 general details(general) , 3 is used for  global (ligands) )
        fing_df=pd.DataFrame(convert_bit_to_list(fing))
        fing_df=fing_df.T
        fingerprint_df=pd.concat([fingerprint_df,fing_df],axis=0)
    return fingerprint_df.reset_index(drop=True)


def Atomtype_fing(smile):
    fingerprint_df=pd.DataFrame()
    for i in smile:
        ms=Chem.MolFromSmiles(i)
        fing =Atomtype_FingerprintMol(ms)
        # fing = AllChem.GetMorganFingerprint(ms,2).ToBitString()   #(1 is to find all internal details(high quality screening) , 2 general details(general) , 3 is used for  global (ligands) )
        # print('fing',fing)
        fing_df=pd.DataFrame(fing)
        fing_df=fing_df.T
        fingerprint_df=pd.concat([fingerprint_df,fing_df],axis=0)
    return fingerprint_df.reset_index(drop=True)


def EStateIndices_fing(smile):
    fingerprint_df=pd.DataFrame()
    for i in smile:
        ms=Chem.MolFromSmiles(i)
        fing =EStateIndices_FingerprintMol(ms)
        # fing = AllChem.GetMorganFingerprint(ms,2).ToBitString()   #(1 is to find all internal details(high quality screening) , 2 general details(general) , 3 is used for  global (ligands) )
        fing_df=pd.DataFrame(fing)
        fing_df=fing_df.T
        fingerprint_df=pd.concat([fingerprint_df,fing_df],axis=0)
    return fingerprint_df.reset_index(drop=True)
import mol_conversion
from rdkit.Chem import PandasTools
import pandas as pd


def mols_2_smiles_test():
    assert sdf_smiles == mol_conversion.mols_2_smiles(sdf_mols)
    
def smiles_2_mols_test():
    assert excel_smiles == mol_conversion.mols_2_smiles(mol_conversion.smiles_2_mols(excel_smiles))
    
def mols_2_inchi_keys_test():
    assert excel_inchi_keys == mol_conversion.mols_2_inchi_keys(mol_conversion.smiles_2_mols(excel_smiles))
    
def smiles_2_inchi_keys_test():
    assert excel_inchi_keys == mol_conversion.smiles_2_inchi_keys(excel_smiles)


# Load excel test data
test_excel = pd.read_excel('excel_test_case.xlsx')
excel_smiles = list(test_excel.Smiles)
excel_inchi_keys = list(test_excel.InChi)

# Load SDF test data
test_sdf = PandasTools.LoadSDF('sdf_test_case.sdf', molColName='Molecule')
sdf_mols = list(test_sdf.Molecule)
sdf_smiles = list(test_sdf.Smiles)

# Run test functions
mols_2_smiles_test()
smiles_2_mols_test()
mols_2_inchi_keys_test()
smiles_2_inchi_keys_test()
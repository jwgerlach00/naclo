import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import PandasTools

from __test_skeleton import BasicTest
import mol_conversion
from df_ops import nan_ops
from df_ops import column_ops
from df_ops import header_ops
from df_ops import duplicates
from df_ops import sdf_tools


class NanTests(BasicTest):
    def __init__(self) -> None:
        self.test_df = pd.DataFrame({'a': ['', '5', '6'], 'b': ['a', '8', 'None']})
        self.blank_solution = pd.DataFrame({'a': [np.nan, '5', '6'], 'b': ['a', '8', np.nan]})
        self.non_num_solution = pd.DataFrame({'a': [np.nan, '5', '6'], 'b': [np.nan, '8', np.nan]})
            
    def convert_to_nan_test(self):
        """Tests convert list of nas to nan."""
        
        # Test replace nan
        na_df = nan_ops.convert_to_nan(self.test_df)
        assert na_df.equals(self.blank_solution)
        
        # Test case sensitivity
        na_df = nan_ops.convert_to_nan(self.test_df, na=('A'))
        assert na_df.equals(pd.DataFrame({'a': ['', '5', '6'], 'b': [np.nan, '8', 'None']}))
        
    def non_num_to_nan(self):
        """Tests ability to convert all non-numeric values to np.nan."""
        
        # Test on only first column
        na_df = nan_ops.non_num_to_nan(self.test_df, self.test_df.columns[0])
        assert na_df.equals(pd.DataFrame({'a': [np.nan, '5', '6'], 'b': ['a', '8', 'None']}))
        
        # Test on only second column (unchanged)
        na_df = nan_ops.non_num_to_nan(self.test_df, self.test_df.columns[1])
        assert na_df.equals(self.test_df)
        
        # Test on both columns
        na_df = nan_ops.non_num_to_nan(self.test_df, self.test_df.columns)
        assert na_df.equals(self.non_num_solution)

    def nan_col_indices_test (self):
        """Tests nan index finder."""
        
        # Test each col separately
        assert nan_ops.nan_col_indices(self.blank_solution, self.blank_solution.columns[0]) == [0]
        assert nan_ops.nan_col_indices(self.blank_solution, self.blank_solution.columns[1]) == [2]
        
    def remove_nan_cols_test(self):
        nan_col_df = self.test_df.copy()
        
        # Append new columns
        nan_col_df['c'] = 3*['None']
        nan_col_df['d'] = 3*['nan']
        nan_col_df['e'] = 3*['']
        nan_col_df['f'] = 3*[np.nan]
        
        # Check new columns were added successfully
        assert len(nan_col_df.columns) == 6
        
        # Test blank column removal
        assert nan_ops.remove_nan_cols(nan_col_df).equals(self.test_df)
        
    def remove_nan_rows_test(self):
        self.non_num_solution = pd.DataFrame({'a': [np.nan, '5', '6'], 'b': [np.nan, '8', np.nan]})
        
        col_a = nan_ops.remove_nan_rows(self.blank_solution, ['a'])
        assert col_a.equals(pd.DataFrame({'a': ['5', '6'], 'b': ['8', np.nan]}))
        
        col_b = nan_ops.remove_nan_rows(self.blank_solution, ['b'])
        assert col_b.equals(pd.DataFrame({'a': [np.nan, '5'], 'b': ['a', '8']}))
        
        both_cols = nan_ops.remove_nan_rows(self.blank_solution, self.blank_solution.columns)
        assert both_cols.equals(pd.DataFrame({'a': ['5'], 'b': ['8']}))


class ColumnOpsTests(BasicTest):
    def __init__(self) -> None:
        self.test_df = pd.DataFrame({'a': [4, '5', '6'], 'b': ['6', '5', '4']})
        self.mol_test_df = pd.DataFrame({
            'SMILES': ['O=C(O)C1(Sc2ccnc3ccc(Br)cc23)CCC1', 
                       'O=C(O)C1(Sc2ccnc3ccc(C4CC4)cc23)CCC1', 
                       'CC(C)(Sc1ccncc1-c1ccc(C#N)c2ccccc12)C(=O)O',
                       None],
            'Molecule': [Chem.MolFromSmiles('O=C(O)C1(Sc2ccnc3ccc(Br)cc23)CCC1'),
                         Chem.MolFromSmiles('O=C(O)C1(Sc2ccnc3ccc(C4CC4)cc23)CCC1'),
                         Chem.MolFromSmiles('CC(C)(Sc1ccncc1-c1ccc(C#N)c2ccccc12)C(=O)O'),
                         None],
            'InChi': [Chem.MolToInchiKey(Chem.MolFromSmiles('O=C(O)C1(Sc2ccnc3ccc(Br)cc23)CCC1')),
                      Chem.MolToInchiKey(Chem.MolFromSmiles('O=C(O)C1(Sc2ccnc3ccc(C4CC4)cc23)CCC1')),
                      Chem.MolToInchiKey(Chem.MolFromSmiles('CC(C)(Sc1ccncc1-c1ccc(C#N)c2ccccc12)C(=O)O')),
                      None]
        })

    def drop_val_test(self):
        # Test number
        dropped1 = column_ops.drop_val(self.test_df, 'a', 4)
        assert dropped1.equals(pd.DataFrame({'a': ['5', '6'], 'b': ['5', '4']}))
        
        # Test string
        dropped2 = column_ops.drop_val(self.test_df, 'a', '4')
        assert dropped2.equals(self.test_df)
        
        # Test b column
        dropped3 = column_ops.drop_val(self.test_df, 'b', '4')
        assert dropped3.equals(pd.DataFrame({'a': [4, '5'], 'b': ['6', '5']}))
        
    def pull_val_test(self):
        # Test number
        pulled1 = column_ops.pull_val(self.test_df, 'a', 4)
        assert pulled1.equals(pd.DataFrame({'a': [4], 'b': ['6']}).astype(object))
        
        # Test string
        pulled2 = column_ops.pull_val(self.test_df, 'a', '6')
        assert pulled2.equals(pd.DataFrame({'a': ['6'], 'b': ['4']}, index=[2]))
        
        # Test b column
        pulled3 = column_ops.pull_val(self.test_df, 'b', '5')
        assert pulled3.equals(pd.DataFrame({'a': ['5'], 'b': ['5']}, index=[1]))
        
    def pull_not_val_test(self):
        # Test number
        pulled1 = column_ops.pull_not_val(self.test_df, 'a', 4)
        assert pulled1.equals(pd.DataFrame({'a': ['5', '6'], 'b': ['5', '4']}, index=[1, 2]).astype(object))
        
        # Test string
        pulled2 = column_ops.pull_not_val(self.test_df, 'a', '6')
        assert pulled2.equals(pd.DataFrame({'a': [4, '5'], 'b': ['6', '5']}))
        
        # Test b column
        pulled3 = column_ops.pull_not_val(self.test_df, 'b', '5')
        assert pulled3.equals(pd.DataFrame({'a': [4, '6'], 'b': ['6', '4']}, index=[0, 2]))
        
        # Test pull invalid number --> returns what was input
        pulled4 = column_ops.pull_not_val(self.test_df, 'b', '3')
        assert pulled4.equals(self.test_df)
        
    def smiles_to_mols_test(self):
        # Test with dropna
        out1 = column_ops.smiles_to_mols(self.mol_test_df, 'SMILES', 'ROMol')
        assert out1.ROMol.map(Chem.MolToSmiles, na_action='ignore').equals(
            self.mol_test_df.Molecule.map(Chem.MolToSmiles, na_action='ignore').dropna())
        
        # Test without dropna
        out2 = column_ops.smiles_to_mols(self.mol_test_df, 'SMILES', 'ROMol', dropna=False)
        assert out2.ROMol.map(Chem.MolToSmiles, na_action='ignore').equals(
            self.mol_test_df.Molecule.map(Chem.MolToSmiles, na_action='ignore'))
        
    def mols_to_inchis_test(self):
        # Test with dropna
        out1 = column_ops.mols_to_inchis(self.mol_test_df, 'Molecule', 'inchi')
        assert out1.inchi.equals(self.mol_test_df.InChi.dropna())
        
        # Test without dropna
        out2 = column_ops.mols_to_inchis(self.mol_test_df, 'Molecule', 'inchi', dropna=False)
        assert out2.inchi.equals(self.mol_test_df.InChi)
        
    def mols_to_smiles_test(self):
        # Test with dropna
        out1 = column_ops.mols_to_smiles(self.mol_test_df, 'Molecule', 'smiles')
        assert out1.smiles.equals(self.mol_test_df.SMILES.dropna())
        
        # Test without dropna
        out2 = column_ops.mols_to_smiles(self.mol_test_df, 'Molecule', 'smiles', dropna=False)
        assert out2.smiles.equals(self.mol_test_df.SMILES)
        
        
class HeaderOpsTests(BasicTest):
    def __init__(self) -> None:
        self.mol_test_df = pd.DataFrame({
            'SMILES': ['O=C(O)C1(Sc2ccnc3ccc(Br)cc23)CCC1', 
                       'O=C(O)C1(Sc2ccnc3ccc(C4CC4)cc23)CCC1', 
                       'CC(C)(Sc1ccncc1-c1ccc(C#N)c2ccccc12)C(=O)O',
                       'c'],
            'Molecule': [Chem.MolFromSmiles('O=C(O)C1(Sc2ccnc3ccc(Br)cc23)CCC1'),
                         Chem.MolFromSmiles('O=C(O)C1(Sc2ccnc3ccc(C4CC4)cc23)CCC1'),
                         Chem.MolFromSmiles('CC(C)(Sc1ccncc1-c1ccc(C#N)c2ccccc12)C(=O)O'),
                         None],
            'InChi': [Chem.MolToInchi(Chem.MolFromSmiles('O=C(O)C1(Sc2ccnc3ccc(Br)cc23)CCC1')),
                      Chem.MolToInchi(Chem.MolFromSmiles('O=C(O)C1(Sc2ccnc3ccc(C4CC4)cc23)CCC1')),
                      Chem.MolToInchi(Chem.MolFromSmiles('CC(C)(Sc1ccncc1-c1ccc(C#N)c2ccccc12)C(=O)O')),
                      None]
        })
        
    def id_nearest_col_test(self):
        # Test similar variations
        assert header_ops.id_nearest_col(self.mol_test_df, 'smiles') == 'SMILES'
        assert header_ops.id_nearest_col(self.mol_test_df, 'smIles') == 'SMILES'
        assert header_ops.id_nearest_col(self.mol_test_df, 'smes') == 'SMILES'
        
        # Test null return
        assert header_ops.id_nearest_col(self.mol_test_df, 's') == None
        
    def remove_header_chars_test(self):
        test_df1 = self.mol_test_df.copy()
        test_df2 = self.mol_test_df.copy()
        test_df3 = self.mol_test_df.copy()
        test_df1.columns = ['MILE', 'Molecule', 'InChi']
        test_df2.columns = ['SMLES', 'Molecule', 'nCh']
        test_df3.columns = ['SMILES', 'Molecule', 'InCh']
        
        # Test w/ case sensitive on
        assert header_ops.remove_header_chars(self.mol_test_df, 's').equals(test_df1)
        assert header_ops.remove_header_chars(self.mol_test_df, 'S').equals(test_df1)
        
        # Test w/ case sensitive off
        assert header_ops.remove_header_chars(self.mol_test_df, 's', case_insensitive=False).equals(self.mol_test_df)
        assert header_ops.remove_header_chars(self.mol_test_df, 'S', case_insensitive=False).equals(test_df1)
        
        # Test multiple instances
        assert header_ops.remove_header_chars(self.mol_test_df, 'i').equals(test_df2)
        assert header_ops.remove_header_chars(self.mol_test_df, 'i', case_insensitive=False).equals(test_df3)
        

class DuplicatesTests(BasicTest):
    def __init__(self) -> None:
        self.test_df = pd.DataFrame({'aa': ['5', '5', '6', 5, 6, '6'], 'b': [1, 2, 3, 4, 5, 6], 
                                     'c': ['a', 'a', 'b', 'b', 'b', 'b']})
        self.filt_not_dup_sol = pd.DataFrame({'aa': ['5', '6', 5, 6], 'b': [1, 3, 4, 5], 'c': ['a', 'b', 'b', 'b']}, 
                                             index=[0, 2, 3, 4])
        self.filt_dup_sol = pd.DataFrame({'aa': ['5', '6'], 'b': [2, 6], 'c': ['a', 'b']}, index=[1, 5])
        self.filt_multi_not_dup_sol = pd.DataFrame({'aa': ['5', '6'], 'b': [1, 3], 'c': ['a', 'b']}, 
                                                   index=[0, 2])
        self.filt_multi_dup_sol = pd.DataFrame({'aa': ['5', 5, 6, '6'], 'b': [2, 4, 5, 6], 
                                                'c': ['a', 'b', 'b', 'b']}, index=[1, 3, 4, 5])
        self.average_sol = pd.DataFrame({'aa': ['5', '6', 5, 6], 'b': [1.5, 4.5, 4, 5], 'c': ['a', 'b', 'b', 'b']}, 
                                        index=[0, 2, 3, 4])
    
    def filter_test(self):
        # Test filtering out duplicates
        assert duplicates.filter(self.test_df, ['aa']).equals(self.filt_not_dup_sol)
        
        # Test filtering out everything EXCEPT duplicates
        assert duplicates.filter(self.test_df, ['aa'], return_duplicates=True).equals(self.filt_dup_sol)
        
        # Test filtering across multiple subsets
        assert duplicates.filter(self.test_df, ['aa', 'c'], how='any').equals(self.filt_multi_not_dup_sol)
        assert duplicates.filter(self.test_df, ['aa', 'c'], how='any', return_duplicates=True).equals(
            self.filt_multi_dup_sol)

    def indices_test(self):
        assert duplicates.indices(self.test_df, 'c') == [1, 3, 4, 5]
        
    def remove_test(self):
        assert duplicates.remove(self.test_df, ['aa']).equals(self.filt_not_dup_sol)
        
    def average_test(self):
        assert duplicates.average(self.test_df, ['aa'], average_by='b').equals(self.average_sol) 
        
        
class SDFToolsTests(BasicTest):
    def __init__(self) -> None:
        self.excel_df = pd.read_excel('excel_test_case.xlsx')
        self.excel_df = column_ops.smiles_to_mols(self.excel_df, 'Smiles', 'Molecule')
        self.sdf_df = PandasTools.LoadSDF('sdf_test_case.sdf', molColName='Molecule')
        self.sdf_df.ID = list(range(4))
    
    def write_sdf_test(self):
        # Write new and load new
        sdf_tools.write_sdf(self.excel_df, 'sdf_test_out.sdf', mol_col_name='Molecule')
        load_sdf = PandasTools.LoadSDF('./sdf_test_out.sdf', molColName='Molecule')
        
        # Ensure Mols are the same
        assert mol_conversion.mols_2_smiles(load_sdf.Molecule) == \
            mol_conversion.mols_2_smiles(self.sdf_df.Molecule)
        
        # Ensure original SMILES are the same
        assert load_sdf.Smiles.equals(self.sdf_df.Smiles)
        
        
nan_tests = NanTests()
column_ops_tests = ColumnOpsTests()
header_ops_tests = HeaderOpsTests()
duplicates_tests = DuplicatesTests()
sdf_tools_tests = SDFToolsTests()

nan_tests.run_all_tests()
column_ops_tests.run_all_tests()
header_ops_tests.run_all_tests()
duplicates_tests.run_all_tests()
sdf_tools_tests.run_all_tests()

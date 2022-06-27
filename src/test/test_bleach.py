from multiprocessing.sharedctypes import Value
import unittest
import json
import pandas as pd
import numpy as np
from naclo import Bleach
import warnings
from rdkit import Chem
import copy


class TestBleach(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        with open('test/default_params.json') as f:
            cls.default_params = json.load(f)
        
        with open('test/default_options.json') as f:
            cls.default_options = json.load(f)
            
        cls.smiles_df = pd.DataFrame({
            'SMILES': ['Cc1cc(/C=C/C#N)cc(C)c1Nc1nc(Nc2ccc(C#N)cc2)ncc1N',
                       'Cc1cc(/C=C/C#N)cc(C)c1Nc1ncc(N)c(Nc2c(C)cc(/C=C/C#N)cc2C)n1',
                       'CCC',
                       'CCC.Cl',
                       'C.NO.S',
                       'Br',
                       '']
        })
        
        warnings.filterwarnings('error')  # Catch in try except
        
        return super().setUpClass()
    
    def test_init_error_checker(self):
        params = copy.deepcopy(self.default_params)
        options = copy.deepcopy(self.default_options)
        
        try:
            Bleach(self.smiles_df, params, options)
        except ValueError as e:
            self.assertEqual(
                e.args[0],
                'NO_STRUCTURE_COLUMN'
            )
        
        params['structure_col'] = 'test'
        try:
            Bleach(self.smiles_df, params, options)
        except ValueError as e:
            self.assertEqual(
                e.args[0],
                'STRUCTURE_COLUMN_NOT_FOUND'
            )
        
        params['structure_col'] = 'SMILES'
        try:
            Bleach(self.smiles_df, params, options)
        except ValueError as e:
            self.assertEqual(
                e.args[0],
                'INVALID_STRUCTURE_TYPE'
            )
        
        # Should run fine
        params['structure_type'] = 'smiles'
        Bleach(self.smiles_df, params, options)
        
    def test_recognized_options_checker(self):
        params = copy.deepcopy(self.default_params)
        options = copy.deepcopy(self.default_options)
        params['structure_col'] = 'SMILES'
        params['structure_type'] = 'smiles'
        options['molecule_settings']['neutralize_charges']['run'] = 42
        
        try:
            Bleach(self.smiles_df, params, options)
        except ValueError as e:
            self.assertEqual(
                list(e.args[0].keys())[0],
                'BAD_OPTION_MOLECULE_SETTINGS_NEUTRALIZE_CHARGES_RUN'
            )

    def test_drop_na(self):
        df = pd.DataFrame({
            'SMILES': [
                pd.NA,
                '',
                'nan',
                'none',
                'CCC',
                np.nan
            ],
            'ROMol': [
                pd.NA,
                '',
                'nan',
                'C',
                'none',
                np.nan
            ],
            'drop_empty': 6*[np.nan]
        })
        
        params = copy.deepcopy(self.default_params)
        
        # SMILES
        params['structure_col'] = 'SMILES'
        params['structure_type'] = 'smiles'
        
        bleach = Bleach(df, params, self.default_options)
        bleach.drop_na()
        
        expected = pd.DataFrame({
            'SMILES': ['CCC']
        }, index=[4])
        
        self.assertEqual(
            True,
            bleach.df.equals(expected)
        )
        
        # ROMol
        params['structure_col'] = 'ROMol'
        params['structure_type'] = 'mol'
        
        bleach = Bleach(df, params, self.default_options)
        bleach.drop_na()
        
        expected = pd.DataFrame({
            'ROMol': ['C']
        }, index=[3])
        
        self.assertEqual(
            True,
            bleach.df.equals(expected)
        )
        
        # ALL_NA_STRUCTURES warning
        params['structure_col'] = 'drop_empty'
        
        bleach = Bleach(df, params, self.default_options)
        
        try:
            bleach.drop_na()
        except RuntimeWarning as w:
            self.assertEqual(
                w.args[0],
                'ALL_NA_STRUCTURES: All structures in specified column were NA, all rows dropped'
            )
            
        # ALL_NA_TARGETS warning
        params['structure_col'] = 'ROMol'
        params['target_col'] = 'drop_empty'
        
        # Not set to drop NA targets
        bleach = Bleach(df, params, self.default_options)
        bleach.drop_na()
        
        expected = pd.DataFrame({
            'ROMol': ['C']
        }, index=[3])
        
        self.assertEqual(
            True,
            bleach.df.equals(expected)
        )
        
        # Set to drop NA targets
        options = copy.deepcopy(self.default_options)
        options['file_settings']['remove_na_targets']['run'] = True
        
        bleach = Bleach(df, params, options)
        
        try:
            bleach.drop_na()
        except RuntimeWarning as w:
            self.assertEqual(
                w.args[0],
                'ALL_NA_TARGETS: All targets in specified column were NA, all rows dropped'
            )
            
    def test_init_structure_compute(self):
        # From SMILES
        params = copy.deepcopy(self.default_params)
        params['structure_col'] = 'SMILES'
        params['structure_type'] = 'smiles'
        
        bleach = Bleach(self.smiles_df, params, self.default_options)
        bleach.drop_na()
        bleach.init_structure_compute()

        [self.assertIsInstance(m, Chem.rdchem.Mol) for m in bleach.df['ROMol']]
        
        # From SMILES w/ bad Mols
        mol_df = bleach.df.copy()
        
        new_row = pd.DataFrame({
            'SMILES': ['CC', 'C', 'CCC'],
            'ROMol': ['test', 1, 'CCC']
        })
        
        mol_df = pd.concat((mol_df, new_row))
        bleach = Bleach(mol_df, params, self.default_options)
        bleach.drop_na()
        bleach.init_structure_compute()
        
        self.assertEqual(
            bleach.df['SMILES'].tolist(),
            mol_df['SMILES'].tolist()
        )
        [self.assertIsInstance(m, Chem.rdchem.Mol) for m in bleach.df['ROMol']]
        
        # From Mols w/ bad Mols
        params['structure_col'] = 'ROMol'
        params['structure_type'] = 'mol'
        bleach = Bleach(mol_df, params, self.default_options)
        bleach.drop_na()
        bleach.init_structure_compute()
        
        [self.assertIsInstance(m, Chem.rdchem.Mol) for m in bleach.df['ROMol']]
        self.assertEqual(
            len(bleach.df),
            6
        )
        self.assertEqual(
            bleach.df['SMILES'].tolist(),
            self.smiles_df['SMILES'].tolist()[:-1]  # Excluding blank
        )
        
    def test_mol_cleanup(self):
        params = copy.deepcopy(self.default_params)
        options = copy.deepcopy(self.default_options)
        params['structure_col'] = 'SMILES'
        params['structure_type'] = 'smiles'
        
        # Salts only
        options['molecule_settings']['remove_fragments']['filter_method'] = 'none'
        
        bleach = Bleach(self.smiles_df, params, options)
        bleach.drop_na()
        bleach.init_structure_compute()
        bleach.mol_cleanup()
        
        expected = pd.DataFrame({
            'SMILES': ['Cc1cc(/C=C/C#N)cc(C)c1Nc1nc(Nc2ccc(C#N)cc2)ncc1N',
                       'Cc1cc(/C=C/C#N)cc(C)c1Nc1ncc(N)c(Nc2c(C)cc(/C=C/C#N)cc2C)n1',
                       'CCC',
                       'CCC',
                       'C.NO.S']
        })
        
        self.assertTrue(
            bleach.df['SMILES'].equals(expected['SMILES'])
        )
        
        # Filters only
        options['molecule_settings']['remove_fragments']['salts'] = ''
        
        # Carbon count
        options['molecule_settings']['remove_fragments']['filter_method'] = 'carbon_count'
        bleach = Bleach(self.smiles_df, params, options)
        bleach.drop_na()
        bleach.init_structure_compute()
        bleach.mol_cleanup()
        self.assertEqual(
            bleach.df['SMILES'].iloc[-2],
            'C'
        )
        
        # Atom count
        options['molecule_settings']['remove_fragments']['filter_method'] = 'atom_count'
        bleach = Bleach(self.smiles_df, params, options)
        bleach.drop_na()
        bleach.init_structure_compute()
        bleach.mol_cleanup()
        self.assertEqual(
            bleach.df['SMILES'].iloc[-2],
            'NO'
        )
        
        # Molecular weight
        options['molecule_settings']['remove_fragments']['filter_method'] = 'mw'
        bleach = Bleach(self.smiles_df, params, options)
        bleach.drop_na()
        bleach.init_structure_compute()
        bleach.mol_cleanup()
        self.assertEqual(
            bleach.df['SMILES'].iloc[-2],
            'S'
        )
        
        # Salts + filter together
        options['molecule_settings']['remove_fragments']['salts'] = 'Cl,Br'
        options['molecule_settings']['remove_fragments']['filter_method'] = 'carbon_count'
        
        bleach = Bleach(self.smiles_df, params, options)
        bleach.drop_na()
        bleach.init_structure_compute()
        bleach.mol_cleanup()
        
        expected = pd.DataFrame({
            'SMILES': ['Cc1cc(/C=C/C#N)cc(C)c1Nc1nc(Nc2ccc(C#N)cc2)ncc1N',
                       'Cc1cc(/C=C/C#N)cc(C)c1Nc1ncc(N)c(Nc2c(C)cc(/C=C/C#N)cc2C)n1',
                       'CCC',
                       'CCC',
                       'C']
        })
        
        self.assertTrue(
            bleach.df['SMILES'].equals(expected['SMILES'])
        )


if __name__ == '__main__':
    unittest.main()

from multiprocessing.sharedctypes import Value
import unittest
import json
import pandas as pd
import numpy as np
from naclo import Bleach
import warnings


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
                       'CCC',
                       'C.NO.S',
                       '']
        })
        
        warnings.filterwarnings('error')  # Catch in try except
        
        return super().setUpClass()
    
    def test_init_error_checker(self):
        params = self.default_params.copy()
        options = self.default_options.copy()
        
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
        params = self.default_params.copy()
        options = self.default_options.copy()
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
        
        params = self.default_params.copy()
        
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
        options = self.default_options.copy()
        options['file_settings']['remove_na_targets']['run'] = True
        
        bleach = Bleach(df, params, options)
        
        try:
            bleach.drop_na()
        except RuntimeWarning as w:
            self.assertEqual(
                w.args[0],
                'ALL_NA_TARGETS: All targets in specified column were NA, all rows dropped'
            )
   
        
        
        
    
    # def test_remove_fragments(self):
    #     params = self.default_params.copy()
    #     options = self.default_options.copy()
    #     params['smiles_col'] = 'SMILES'
        
    #     # Carbon count
    #     options['molecule_settings']['remove_fragments']['filter_method'] = 'carbon_count'
    #     bleach = Bleach(self.smiles_df, params, options)
    #     df = bleach.main()
    #     self.assertEqual(
    #         df['SMILES'].iloc[-1],
    #         'C'
    #     )
        
    #     # Atom count
    #     options['molecule_settings']['remove_fragments']['filter_method'] = 'atom_count'
    #     bleach = Bleach(self.smiles_df, params, options)
    #     df = bleach.main()
    #     self.assertEqual(
    #         df['SMILES'].iloc[-1],
    #         'NO'
    #     )
        
    #     # Molecular weight
    #     options['molecule_settings']['remove_fragments']['filter_method'] = 'mw'
    #     bleach = Bleach(self.smiles_df, params, options)
    #     df = bleach.main()
    #     self.assertEqual(
    #         df['SMILES'].iloc[-1],
    #         'S'
    #     )


if __name__ == '__main__':
    unittest.main()

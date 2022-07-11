import unittest
import pandas as pd
import numpy as np
from copy import deepcopy

from naclo import binarize_default_params, binarize_default_options
from naclo import Binarize


class TestBinarize(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.default_params = binarize_default_params
        cls.default_options = binarize_default_options
        
        cls.test_df = pd.DataFrame({
            'smiles': [
                'CCC',
                'C',
                'CN=C=O',
                'CN(C)C.Cl',
                None
            ],
            'target': [
                55,
                4,
                7,
                100,
                2000
            ],
            'units': [
                'ug•ml-1',
                'mg/l',
                'unrecognized',
                np.nan,
                'pm'
            ],
            'qualifiers': [
                '>',
                '<',
                '=',
                '>',
                '<'
            ]
        })
        
        cls.default_params['structure_col'] = 'smiles'
        cls.default_params['structure_type'] = 'smiles'
        cls.default_params['target_col'] = 'target'
        cls.default_params['decision_boundary'] =  7  # neg log molar
        
        return super().setUpClass()
    
    def test_convert_units(self):
        options = deepcopy(self.default_options)
        options['convert_units']['units_col'] = 'units'
        binarize = Binarize(self.test_df, params=self.default_params, options=options)
        molar = binarize.convert_units()
        self.fail()
    
    def test_binarize(self):
        options = deepcopy(self.default_options)
        
        # No qualifiers, active boundary
        options['qualifiers']['run'] = False
        options['active_boundary_cond'] = True
        binarize = Binarize(self.test_df, params=self.default_params, options=options)
        self.assertTrue(
            np.array_equal(
                binarize.binarize(),
                np.array([0, 1, 1, 0])
            )
        )
        
        # No qualifiers, inactive boundary
        options['qualifiers']['run'] = False
        options['active_boundary_cond'] = False
        binarize = Binarize(self.test_df, params=self.default_params, options=options)
        self.assertTrue(
            np.array_equal(
                binarize.binarize(),
                np.array([0, 1, 0, 0])
            )
        )
        
        # Qualifiers, active boundary
        options['qualifiers']['run'] = True
        options['qualifiers']['qualifier_col'] = 'qualifiers'
        options['active_boundary_cond'] = True
        binarize = Binarize(self.test_df, params=self.default_params, options=options)
        self.assertTrue(
            np.array_equal(
                binarize.binarize(),
                np.array([0, 1, 1, 0])
            )
        )
        
        # Qualifiers, inactive boundary
        options['qualifiers']['run'] = True
        options['qualifiers']['qualifier_col'] = 'qualifiers'
        options['active_boundary_cond'] = False
        binarize = Binarize(self.test_df, params=self.default_params, options=options)
        self.assertTrue(
            np.array_equal(
                binarize.binarize(),
                np.array([0, 1, 0, 0])
            )
        )
        
    def test_main(self):
        options = deepcopy(self.default_options)
        options['qualifiers']['run'] = True
        options['qualifiers']['qualifier_col'] = 'qualifiers'
        options['convert_units']['units_col'] = 'units'
        binarize = Binarize(self.test_df, params=self.default_params, options=options)
        x = binarize.main()
        

if __name__ == '__main__':
    unittest.main()
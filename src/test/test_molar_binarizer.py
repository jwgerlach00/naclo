import unittest
import numpy as np

from naclo import MolarBinarizer


class TestBinarize(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.test_values = [
            55,
            4,
            7,
            100,
            2000
        ]
        cls.molar_qualifiers = [
            '<',
            '>',
            '<',
            '=',
            '='
        ]
        
        cls.molar_boundary = 10
        
        cls.qualifier_binarizer = MolarBinarizer(molar_vals=cls.test_values, molar_boundary=cls.molar_boundary,
                                                 molar_qualifiers=cls.molar_qualifiers)
        cls.no_qualifier_binarizer = MolarBinarizer(molar_vals=cls.test_values, molar_boundary=cls.molar_boundary)
        return super().setUpClass()
    
    def test_binarize(self):
        qualifier_out = self.qualifier_binarizer.binarize()
        no_qualifier_out = self.no_qualifier_binarizer.binarize()
        
        expected_qualifier = np.array([
            np.nan,
            np.nan,
            1,
            0,
            0
        ])
        expected_no_qualifier = np.array([
            0,
            1,
            1,
            0,
            0
        ])
        
        self.assertTrue(
            np.allclose(
                qualifier_out,
                expected_qualifier,
                equal_nan=True
            )
        )
        self.assertTrue(
            np.array_equal(
                no_qualifier_out,
                expected_no_qualifier
            )
        )
        
if __name__ == '__main__':
    unittest.main()

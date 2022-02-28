from test_skeleton import BasicTest
from mol_stats import *
import mol_conversion


class MolStatsTests(BasicTest):
    def __init__(self) -> None:
        self.test_smile = 'CN=C=O'
        self.test_mols = mol_conversion.smiles_to_mols([self.test_smile, 'C'])
    
    def carbon_num_test(self):
        assert carbon_num(self.test_smile) == 2
    
    def mol_weights_test(self):
        assert list(map(lambda x: round(x, 3), mol_weights(self.test_mols))) == [57.021, 16.031]


if __name__ == '__main__':
    mol_stats_tests = MolStatsTests()
    mol_stats_tests.run_all_tests()

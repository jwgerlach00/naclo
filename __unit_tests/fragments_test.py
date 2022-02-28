import fragments
from mol_conversion import smiles_to_mols, smiles_from_mols


frag_rm = fragments.FragRemover()

def frag_remover_test():
    # Test case
    test_smile = 'CCC.C.C'
    
    # Test all frag removing methods
    assert frag_rm.mw(test_smile) == 'CCC'
    assert frag_rm.atom_count(test_smile) == 'CCC'
    assert frag_rm.carbon_count(test_smile) == 'CCC'
    
def remove_salts_test():
    # Test case
    test_mols = smiles_to_mols(['CN(C)C.Cl', 'CN(C)C.N'])
    
    # Test default salts [Cl,Br]
    assert smiles_from_mols(fragments.remove_salts(test_mols)) == ['CN(C)C', 'CN(C)C.N']
    
    # Test custom salts [N]
    assert smiles_from_mols(fragments.remove_salts(test_mols, salts='[N]')) == ['CN(C)C.Cl', 'CN(C)C']


if __name__ == '__main__':
    frag_remover_test()
    remove_salts_test()

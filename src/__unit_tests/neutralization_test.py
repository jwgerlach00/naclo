from neutralization import *
import mol_conversion


def init_neutralization_rxns_test():
    rxns = init_neutralization_rxns()
    rxns = list(rxns.values()) +  list(rxns.keys())
    assert mol_conversion.mols_2_smiles(rxns) == ['O', '*']  # ????
    

def neutralize_charges_test():
    # Doesn't do anything yet
    pass


init_neutralization_rxns_test()
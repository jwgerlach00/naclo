from rdkit.Chem import SaltRemover


def remove_salts(mols):
    """Removes salts from pandas dataframe.

    Args:
        mols (list): Contains RDKit Mol objects.

    Returns:
        list: Contains RDKit Mol objects with salts removed.
    """
    remover = SaltRemover()
    
    new_mols = []
    for mol in mols:
        m = remover.StripMol(mol)
        new_mols.append(m)
        
    return new_mols

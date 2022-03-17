from rdkit import Chem


def mols_2_smiles(mols):  # *
    """Generates SMILES strings from list of rdkit Mol objects.

    Args:
        mols (iter[rdkit Mol]): Contains RDKit Mols.

    Returns:
        list[str]: Contains SMILES strings.
    """
    return [Chem.MolToSmiles(mol) for mol in mols]

def smiles_2_mols(smiles):  # *
    """Generates rdkit Mol objects from SMILES strings.

    Args:
        smiles (iter[str]): Contains SMILES strings.

    Returns:
        list[rdkit Mol]: Contains RDKit Mols.
    """
    return [Chem.MolFromSmiles(smile) for smile in smiles]

def mols_2_inchi_keys(mols):  # *
    """Generates InChI key strings from rdkit Mol objects.

    Args:
        mols (iter[rdkit Mols]): Contains rdkit Mols.

    Returns:
        list[str]: Contains InChI key strings.
    """
    return [Chem.MolToInchiKey(mol) for mol in mols]

def smiles_2_inchi_keys(smiles):  # *
    """Generates InChI key strings from SMILES strings.

    Args:
        smiles (iter[str]): Contains SMILES strings.

    Returns:
        list[str]: Contains InChI key strings.
    """
    mols = smiles_2_mols(smiles)
    return mols_2_inchi_keys(mols)

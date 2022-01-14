from rdkit import Chem


def smiles_from_mols(mols):
    """Generates list of SMILES from list of Mols.

    Args:
        mols (list): Contains RDKit Mol objects.

    Returns:
        list: Contains SMILES strings.
    """
    return [Chem.MolToSmiles(mol) for mol in mols]

def smiles_to_mols(smiles):
    """Generates list of Mols from list of SMILES.

    Args:
        smiles (list): Contains SMILES strings.

    Returns:
        list: Contains RDKit Mol objects.
    """
    return [Chem.MolFromSmiles(smile) for smile in smiles]

def inchis_from_mols(mols):
    """Generates list of InChI strings from list of Mols.

    Args:
        mols (list): Contains rdkit Mols.

    Returns:
        list: Contains InChI strings.
    """
    return [Chem.rdinchi.MolToInchiKey(mol) for mol in mols]

def inchis_to_mols(inchis):
    """Generates list of Mols from list of InChI strings.

    Args:
        inchis (list): Contains InChI strings.

    Returns:
        list: Contains rdkit Mols.
    """
    return [Chem.inchi.MolFromInchi(inchi) for inchi in inchis]

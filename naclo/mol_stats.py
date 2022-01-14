from rdkit.Chem.Descriptors import ExactMolWt


def find_molecular_weights(mols):
    """Finds molecular weights from list of Mols.

    Args:
        mols (list): Contains RDKit Mol objects.

    Returns:
        list: Molecular weights. Numeric.
    """
    return [ExactMolWt(mol) for mol in mols]

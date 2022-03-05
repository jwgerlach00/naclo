from rdkit import Chem


def drop_val(df, col_name, value):  # *
    """Drops row from df where column entry is equal to value.

    Args:
        df (pandas DataFrame): DataFrame containing rows to drop.
        col_name (str): Name of column to investigate.
        value (pandas object): Value to drop instances of.

    Returns:
        pandas DataFrame: DataFrame with rows dropped.
    """
    return pull_not_val(df, col_name, value).reset_index(drop=True)

def pull_val(df, col_name, value):  # *
    """Retrieves rows from df where column entry is equal to value.

    Args:
        df (pandas DataFrame): DataFrame containing rows to pull.
        col_name (str): Name of column to investigate.
        value (pandas object): Value to pull instances of.

    Returns:
        pandas DataFrame: DataFrame where value is found.
    """
    return df.where(df[col_name] == value).dropna()

def pull_not_val(df, col_name, value):  # *
    """Retrieves rows from df where column entry is not equal to value.

    Args:
        df (pandas DataFrame): DataFrame containing rows to pull.
        col_name (str): Name of column to investigate.
        value (pandas object): Value to not pull instances of.

    Returns:
        pandas DataFrame: DataFrame where value is not found.
    """
    return df.where(df[col_name] != value).dropna()

def smiles_to_mols(df, smiles_name, mol_name, dropna=True):  # *
    """Adds rdkit Mol column to df using SMILES column as reference.

    Args:
        df (pandas DataFrame): DataFrame to add Mol column to.
        smiles_name (str): Name of SMILES column in df.
        molecule_name (str): Name of Mol column in df.
        dropna (bool, optional): Drop NA Mols. Defaults to True.

    Returns:
        pandas DataFrame: DataFrame with Mol column appended.
    """
    df[mol_name] = df[smiles_name].map(Chem.MolFromSmiles, na_action='ignore')  # Results in NA el if fails
    return df.dropna(subset=[mol_name]) if dropna else df

def mols_to_inchis(df, mol_name, inchi_name, dropna=True):  # *
    """Adds InChi Key column to df using Mol column as reference.

    Args:
        df (pandas DataFrame): DataFrame to add InChi column to.
        mol_name (str): Name of InChi column in df.
        inchi_name (str): Name of InChi column in df.
        dropna (bool, optional): Drop NA InChis. Defaults to True.

    Returns:
        pandas DataFrame: DataFrame with InChi column appended.
    """
    df[inchi_name] = df[mol_name].map(Chem.MolToInchiKey, na_action='ignore')  # Results in NA el if fails
    return df.dropna(subset=[inchi_name]) if dropna else df
    
def mols_to_smiles(df, mol_name, smiles_name, dropna=True):  # *
    """Adds SMILES Key column to df using Mol column as reference.

    Args:
        df (pandas DataFrame): DataFrame to add SMILES column to.
        mol_name (str): Name of Mol column in df.
        inchi_name (str): Name of InChi column in df.
        dropna (bool, optional): Drop NA SMILES. Defaults to True.

    Returns:
        pandas DataFrame: DataFrame with SMILES column appended.
    """
    df[smiles_name] = df[mol_name].map(Chem.MolToSmiles, na_action='ignore')  # Results in NA el if fails
    return df.dropna(subset=[smiles_name]) if dropna else df
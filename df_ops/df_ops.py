import numpy as np
from rdkit.Chem import rdinchi


def case_insensitive_column_names(df, names):
    """Reformats df column names to case-match column in df.

    Args:
        df (pandas DataFrame): DataFrame to perform column search on.
        names (list): Names of columns to search for.

    Returns:
        list: DataFrame column names matching names.
    """
    df_cols = df.columns
    low_df_cols = list(map(str.lower, df_cols))  # Convert df column names to lowercase for comparison

    for i, name in enumerate(names):
        low_name = name.lower()  # Convert input name to lowercase for comparison

        # Find existing column name
        if low_name in low_df_cols:
            column_name_index = low_df_cols.index(low_name)
            column_name = df_cols[column_name_index]
        else:
            column_name = low_name

        # Reassign names
        names[i] = column_name

    return names

def remove_na_columns(df, na=('', 'nan', 'None')):
    """Removes columns that are ENTIRELY blank or NaN.

    Args:
        df (pandas DataFrame): DataFrame containing columns to drop.
        na (tuple, optional): Values to drop in addition to np.nan. Defaults to ('', 'nan', 'None').

    Returns:
        pandas DataFrame: DataFrame with columns dropped.
    """
  
    [df.replace(n, np.nan, inplace=True) for n in na]
    df.dropna(axis=1, how="all", inplace=True)  # drop entire column if ALL NaN
    return df

def drop_value_from_column(df, col_name, value):
    """Drops row from df where column entry is equal to value.

    Args:
        df (pandas DataFrame): DataFrame containing rows to drop.
        col_name (str): Name of column to investigate.
        value (pandas object): Value to drop instances of.

    Returns:
        pandas DataFrame: DataFrame with rows dropped.
    """
    return df.where(df[col_name] != value).dropna()

def pull_value_from_column(df, col_name, value):
    """Retrieves row from df where column entry is not equal to value.

    Args:
        df (pandas DataFrame): DataFrame containing rows to pull.
        col_name (str): Name of column to investigate.
        value (pandas object): Value to pull instances of.

    Returns:
        pandas DataFrame: DataFrame where value is found.
    """
    return df.where(df[col_name] == value).dropna()

def remove_na_rows(df, value_column_names, na=('', 'nan', 'None')):
    """Removes entrie rows from df where value is blank or NaN.

    Args:
        df ([type]): [description]
        value_column_name ([type]): [description]

    Returns:
        [type]: [description]
    """
    [df.replace(n, np.nan, inplace=True) for n in na]
    return df.dropna(axis=0, subset=value_column_names)

def remove_header_chars(df, chars):
    """Removes chars in string from df header.

    Args:
        df (pandas DataFrame): DataFrame containing header with chars to remove.
        chars (str): Continuous string of chars to individually remove.

    Returns:
        pandas DataFrame: DataFrame with header chars removed.
    """
    for column_name in df:
        translate_remove = column_name.maketrans('', '', chars)
        new_column_name = column_name.translate(translate_remove)
        df.rename(columns={column_name: new_column_name}, inplace=True)

    return df

def add_molecule_column(df, smiles_name, molecule_name):
    """Adds rdkit Mol column to df using SMILES column as reference.

    Args:
        df (pandas DataFrame): DataFrame to add Mol column to.
        smiles_name (str): Name of SMILES column in df.
        molecule_name (str): Name of Mol column in df.

    Returns:
        pandas DataFrame: DataFrame with Mol column appended.
    """
    mols = []

    for i, smiles in enumerate(df[smiles_name]):
        try:
            mols.append(Chem.MolFromSmiles(smiles))
        except TypeError:
            df.drop(index=i, inplace=True)

    df[molecule_name] = mols

    return df

"""
    Re-loads mols from SMILES (converts to canonical), strips salts, adds inchi column
    :param df:
    :param column_name:
    :param rm_salts: bool, remove salts or not
    :return: list of molecules, list of smiles, list of inchis, list of drop indices
    """

def mol_2_smiles_2_inchi(df, mol_col_name):
    """"""

    drop_indices = []
    smiles = []
    molecules = []
    inchis = []

    for i, mol in enumerate(df[mol_col_name]):
        # Standardize to remove tautomers REMOVED BC TAKES TOO LONG (2 VERSIONS)
        # smiles = MolStandardize.canonicalize_tautomer_smiles(smiles)
        # tautomer_enumerator = MolStandardize.rdMolStandardize.TautomerEnumerator()
        try:
            # Append SMILES, inchi, mol columns
            smiles.append(Chem.MolToSmiles(mol))
            inchis.append(rdinchi.MolToInchiKey(mol))
            molecules.append(mol)
        except:
            drop_indices.append(i)
    return molecules, smiles, inchis, drop_indices

from rdkit import Chem
from rdkit.Chem import PandasTools


def write_sdf(df, out_path, mol_col_name, id_column_name='ID'):
    """Writes dataframe to SDF file. Includes ID name in ID is valid.

    Args:
        df (pandas DataFrame): Data to write.
        out_path (str or file-like): Path to save SDF to.
        mol_col_name (str): Name of molecule column in dataframe.
        id_column_name (str, optional): Name of id column. Defaults to 'ID'.
    """
    try:
        PandasTools.WriteSDF(df, out_path, molColName=mol_col_name, properties=df.columns, idName=id_column_name)
    except KeyError:
        PandasTools.WriteSDF(df, out_path, molColName=mol_col_name, properties=df.columns)

def add_smiles_column_to_sdf(df, smiles_name, molecule_name):
    """Adds SMILES column to df.

    Args:
        df (pandas DataFrame): DataFrame to add SMILES column to.
        smiles_name (str): Name of new SMILES column.
        molecule_name (str): Name of Mol column in df.

    Returns:
        pandas DataFrame: DataFrame with SMILES appended.
    """

    new_smiles_col = []
    for mol in df[molecule_name]:
        smiles = Chem.MolToSmiles(mol)
        new_smiles_col.append(smiles)

    df[smiles_name] = new_smiles_col

    return df

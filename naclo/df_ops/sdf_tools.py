from warnings import warn
from rdkit.Chem import PandasTools


def write_sdf(df, out_path, mol_col_name, id_column_name='RowID'):  # *
    """Writes dataframe to SDF file. Includes ID name if ID is valid.

    Args:
        df (DataFrame): Data to write.
        out_path (str or file-like): Path to save SDF to.
        mol_col_name (str): Name of molecule column in dataframe.
        id_column_name (str, optional): Name of id column. Defaults to 'ID'.
    """
    try:
        PandasTools.WriteSDF(df, out_path, molColName=mol_col_name, properties=df.columns, idName=id_column_name)
    except KeyError:
        PandasTools.WriteSDF(df, out_path, molColName=mol_col_name, properties=df.columns, idName='RowID')
        warn(f'write_sdf \'{id_column_name}\' ID name invalid', UserWarning)

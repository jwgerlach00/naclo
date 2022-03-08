import numpy as np
import re


"""
    File name: nan_ops.py
    Author: Jacob Gerlach
    Description: Assortment of basic pandas DataFrame cleaning operations to handle null values.
    Notes:
        * = Function has associated unit test.
"""


def convert_to_nan(df, na=('', 'nan', 'none')):  # *
    """Converts all instances found in na to np.nan. Case insensitive.

    Args:
        df (pandas DataFrame): Input DataFrame
        na (iterable, optional): Contains patterns to convert. Defaults to ('', 'nan', 'none').

    Returns:
        pandas DataFrame: DataFrame replaced to nan.
    """
    na = [s.lower() for s in na]
    return df.applymap(lambda s: np.nan if str(s).lower() in na else s)

def non_num_to_nan(df, columns):  # *
    """Removes all non numeric (0-9) characters from columns specified.

    Args:
        df (pandas DataFrame): Contains non-numeric values.

    Returns:
        pandas DataFrame: Non-numerics replaced w np.nan.
    """
    
    # Cast to dataframe in case series (1 column)
    try:
        df = df[columns].applymap(lambda s: re.sub('[^0-9]', '', str(s))).replace('', np.nan)
    except AttributeError: 
        df[columns] = df[columns].map(lambda s: re.sub('[^0-9]', '', s)).replace('', np.nan)
        
    return df

def nan_col_indices(df, column_name):  # *
    # replace_blanks_w_na(df)
    return df[df[column_name].isnull()].index.tolist()

def remove_nan_cols(df):  # *
    """Removes columns that are ENTIRELY blank or NaN.

    Args:
        df (pandas DataFrame): DataFrame containing columns to drop.

    Returns:
        pandas DataFrame: DataFrame with columns dropped.
    """
    return convert_to_nan(df).dropna(axis=1, how="all")

def remove_nan_rows(df, value_column_names):  # *
    """Removes entire rows from df where value columns are blank or NaN.

    Args:
        df (pandas DataFrame): DataFrame containing column with NaN values.
        value_column_name (iterable): Contains column names to look for.

    Returns:
        pandas DataFrame: DataFrame with rows dropped where value is NaN.
    """
    return df.dropna(axis=0, subset=value_column_names).reset_index(drop=True)

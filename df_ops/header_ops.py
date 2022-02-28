from difflib import get_close_matches


"""
    File name: header_ops.py
    Author: Jacob Gerlach
    Description: Assortment of basic pandas DataFrame cleaning operations.
    Notes:
        * = Function has associated unit test.
"""


def id_nearest_col(df, name, case_insensitive=True):  # *
    """Sends closest column to input column name provided that it is at least 50% similar. CASE INSENSITIVE.

    Args:
        df (pandas DataFrame): DataFrame containing header with columns of interest.
        name (str): Column name to search for.
        case_insensitive (bool, optional): If True, converts both to lower before comparing. Defaults to True.

    Returns:
        str: Most similar column name.
    """
    similarity_ratio = 0.5  # 50%
    
    # Convert both to lowercase
    if case_insensitive:
        df_head = [col.lower() for col in df]
        name = name.lower()
    else:
        df_head = df.columns
    
    # Find closest match
    matches = get_close_matches(name, df_head, n=1, cutoff=similarity_ratio)
    if len(matches) <= 0:
        return None
    else:
        match_index = df_head.index(matches[0])
        return df.columns[match_index]

def remove_header_chars(df, chars, case_insensitive=True):  # *
    """Removes chars in string from df header.

    Args:
        df (pandas DataFrame): DataFrame containing header with chars to remove.
        chars (str): Continuous string of chars to individually remove.
        case_insensitive (bool, optional): If True, converts both to lower before comparing. Defaults to True.

    Returns:
        pandas DataFrame: DataFrame with header chars removed.
    """
    df = df.copy()
    
    # Convert both to lowercase
    if case_insensitive:
        chars += chars.upper() + chars.lower()
    
    for column_name in df.columns:
        translate_remove = column_name.maketrans('', '', chars)
        new_column_name = column_name.translate(translate_remove)
        df.rename(columns={column_name: new_column_name}, inplace=True)

    return df
import pandas as pd


def filter_duplicates(df, subsets, return_duplicates=False):
    """Filters duplicates to either return DataFrame subset of duplicates or non-duplicates.

    Args:
        df (pandas DataFrame): DataFrame to perform duplicate search on.
        subsets (list): List of column names to compare in duplicate search. \
            All must match to be considered a duplicate.
        return_duplicates (bool, optional): If True, returns non-duplicate df else returns duplicate df. \
            Defaults to False.

    Returns:
        pandas DataFrame: Filtered df.
    """

    # Filter duplicates by inchi, break into 'drop' and 'keep'
    duplicates = pd.DataFrame.duplicated(df, subset=subsets)
    
    return df[duplicates] if return_duplicates else df[~duplicates]

def remove_duplicates(df, subsets):
    """Removes duplicates from df.

    Args:
        df (pandas DataFrame): DataFrame to perform dupicate search on.
        subsets (list): List of column names to compare in duplicate search. \
            All must match to be considered a duplicate.

    Returns:
        pandas DataFrame: DataFrame with duplicates removed.
    """
    return filter_duplicates(df, subsets=subsets)

def average_duplicates(df, subsets, average_by):
    """Averages duplicates from df.

    Args:
        df (pandas DataFrame): DataFrame to perform duplicate search on.
        subsets (list): List of column names to compare in duplicate search. \
            All must match to be considered a duplicate.
        average_by (str): Column name to average along.

    Returns:
        pandas DataFrame: DataFrame with duplicates averaged.
    """
    out = filter_duplicates(df, subsets, return_duplicates=True)
        
    # Ensure average_by exists in the dataframe
    assert(average_by in df.columns)
    
    # Cast value column to float (for avg computation)
    df[subsets] = df[subsets].astype(float)

    # Group duplicate rows and aggregate by mean (adds mean to end of value column)
    avg_value = df.groupby(by=subsets).agg({average_by: ['mean']})

    # Drop original value column ('value' --> 'valuemean')
    out.drop(columns=[average_by], inplace=True)  # drop original values

    # Add averaged values to dataframe
    out = pd.merge(out, avg_value, on=subsets)
    
    return out

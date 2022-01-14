import pandas as pd


def filter_duplicates(df, subset, return_duplicates=False):
    """
    Handles duplicate rows according to selection. Averages, removes, or keeps. If 'None' is selected then averaging
    is not allowed.
    """

    # Filter duplicates by inchi, break into 'drop' and 'keep'
    duplicates = pd.DataFrame.duplicated(df, subset=subset)
    
    return df[duplicates] if return_duplicates else df[~duplicates]

def remove_duplicates(df, subset):
    return filter_duplicates(df, subset)

def average_duplicates(df, subset, average_by):
    out = filter_duplicates(df, subset, return_duplicates=True)
        
    # Ensure average_by exists in the dataframe
    assert(average_by in df.columns)
    
    # Cast value column to float (for avg computation)
    df[subset] = df[subset].astype(float)

    # Group duplicate rows and aggregate by mean (adds mean to end of value column)
    avg_value = df.groupby(by=subset).agg({average_by: ['mean']})

    # Drop original value column ('value' --> 'valuemean')
    out.drop(columns=[average_by], inplace=True)  # drop original values

    # Add averaged values to dataframe
    out = pd.merge(out, avg_value, on=subset)
    
    return out
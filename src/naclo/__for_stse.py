from codecs import xmlcharrefreplace_errors
from typing import List, Tuple, Union, Iterable
import pandas as pd
from numpy import iterable, ma
import numpy as np


def sync_na_drop(df:Union[pd.DataFrame, pd.Series], col_s:Union[str, list], *iterables:Iterable,
                 all_na:bool=True) -> Tuple[pd.DataFrame, List[Iterable]]:
    """Syncs NaN values in a dataframe with a list of iterables.

    Args:
        df (pd.DataFrame): Data to drop from.
        col_s (Union[str, list]): Name of single column or list of columns to subset for drop.
        all_na (bool): True, drops row only if all values in row are NaN. Else drops row if any values is NaN.
        Defaults to True.

    Returns:
        Tuple[List[Iterable], pd.DataFrame]: Tuple of list of masked iterables and NA dropped DataFrame.
    """
    f = df[col_s].isnull
    if isinstance(df[col_s], pd.DataFrame):
        is_na_mask = f().all(axis=1) if all_na else f().any(axis=1).to_numpy()
    elif isinstance(df[col_s], pd.Series):
        is_na_mask = f().to_numpy()
        
    if len(iterables) == 1:
        out = np.array(iterables[0])[~is_na_mask]
    else:
        out = [np.array(it)[~is_na_mask] for it in iterables]
        
    return (df.dropna(subset=[col_s], how=('all' if all_na else 'any')),
            out)
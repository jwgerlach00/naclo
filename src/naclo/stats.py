import numpy as np
from rdkit import DataStructs


def sim_matrix(row_keys, col_keys, key_type='ecfp'):
    """Constructs a tanimoto matrix from ecfp or maccs fingerprints.

    :param row_keys: Collection of fingerprints
    :type row_keys: iter[binary]]
    :param col_keys: Collection of fingerprints 2
    :type col_keys: iter[binary]
    :param key_type: Fingerprint type, defaults to 'ecfp'
    :type key_type: str, optional
    :return: Tanimoto array
    :rtype: np.array()
    """
    
    # Initialize zero array
    full_array = np.zeros((len(row_keys), len(col_keys)))

    # Run through rows and cols of matrix according to key_type
    for i, row in enumerate(row_keys):
        for j, col in enumerate(col_keys):
            if key_type == 'ecfp':
                similarity = DataStructs.TanimotoSimilarity(row, col)
            elif key_type == 'maccs':
                similarity = DataStructs.FingerprintSimilarity(row, col)
            else:
                raise('Key type must be either "ecfp" or "maccs".')
                
            full_array[i][j] = similarity
    
    return full_array

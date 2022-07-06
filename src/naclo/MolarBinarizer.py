import numpy as np
from typing import Iterable, Union, Optional
from copy import copy

from naclo import UnitConverter


class MolarBinarizer:
    def __init__(self, molar_vals:Iterable, molar_boundary:Union[float, int, np.number],
                 molar_qualifiers:Optional[Iterable]=None, active_boundary_cond:bool=True) -> None:
        self.__molar_vals = copy(molar_vals)
        self.__molar_boundary = molar_boundary
        self.__molar_qualifiers = copy(molar_qualifiers)
        self.__active_boundary_cond = active_boundary_cond
    
    @staticmethod
    @np.vectorize
    def __bin_molar_qualifiers(molar_val:Union[int, float, np.number], molar_boundary:Union[int, float, np.number],
                               relation:str, active_boundary_cond:bool=True) -> Union[int, float]:
        """Categorizes a molarity as active or inactive using qualifiers. Molarity is considered active when it is below
        a threshold suggesting that less compound is required to obtain a response.

        Args:
            molar_val (Union[int, float, np.number]): 
            relation (str): Qualifier.
            molar_boundary (Union[int, float, np.number]): Binary decision boundary.
            active_boundary_cond (bool, optional): Sets equal qualifiers to 1 (True) or 0 (False). Defaults to True.

        Returns:
            int, float: 1 (active), 0 (inactive), or np.nan for unload relation or qualifier conflict.
        """
        if molar_val == molar_boundary:  # Boundary case
            return 1 if active_boundary_cond else 0
        
        elif relation in ['>', '≥']:
            return 0 if molar_val > molar_boundary else np.nan  # NaN if less than boundary (qualifier conflict)

        elif relation in ['<', '≤']:
            return 1 if molar_val < molar_boundary else np.nan  # NaN if greater than boundary (qualifier conflict) \
                # else inactive
                
        elif relation == '=':
            return 1 if molar_val < molar_boundary else 0

        else:  # Relation unknown
            return np.nan
    
    @staticmethod
    @np.vectorize
    def __bin_molar_no_qualifiers(molar_val, molar_boundary, active_boundary_cond):
        if molar_val == molar_boundary:
            return 1 if active_boundary_cond else 0
        else:
            return 1 if molar_val < molar_boundary else 0
    
    def binarize(self) -> np.array:
        if self.__molar_qualifiers:
            return self.__bin_molar_qualifiers(molar_val=self.__molar_vals, molar_boundary=self.__molar_boundary,
                                               relation=self.__molar_qualifiers,
                                               active_boundary_cond=self.__active_boundary_cond)
        else:
            return self.__bin_molar_no_qualifiers(molar_val=self.__molar_vals, molar_boundary=self.__molar_boundary,
                                                  active_boundary_cond=self.__active_boundary_cond)

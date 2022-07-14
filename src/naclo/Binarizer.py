from cmath import nan
import numpy as np
from typing import Iterable, Union, Optional
from copy import copy

class Binarizer:
    __less_than = ['<', '≤']
    __greater_than = ['>', '≥']
    
    def __init__(self, values:Iterable, boundary:Union[float, int, np.number], active_operator:str,
                 qualifiers:Optional[Iterable]=None) -> None:
        self.__values = copy(values)
        self.__boundary = boundary
        self.__active_operator = active_operator
        self.__qualifiers = copy(qualifiers)
    
    @staticmethod
    def __equals_edge_case(relation:str, active_operator:str, nan_greater_equal:bool=True) -> Union[int, float]:
        """Binarization handler for when value == threshold.

        Args:
            relation (str): Value qualifier. Set to one of: =, <, >, ≤, or ≥.
            active_operator (str): Condition to set values to active based on threshold. Set to one of: <, >, <=, or >=.
            nan_greater_equal (bool): Conditional to set any '≥' relation to np.nan. Defaults to True.

        Returns:
            Union[int, float]: 1 (active), 0 (inactive), or np.nan for '≥' special case.
        """
        if relation == '≥' and nan_greater_equal:  # Special case --> throw out bc ambiguous
            return np.nan
        return 1 if '=' in active_operator else 0
    
    @staticmethod
    def __active_operator_handler(is_greater_than:bool, active_operator:str) -> int:
        """Binarizes according to operator and (value > threshold) conditional.

        Args:
            is_greater_than (bool): Conditional for value > threshold.
            active_operator (str): Condition to set values to active based on threshold. Set to one of: <, >, <=, or >=.

        Returns:
            int: 1 (active) or 0 (inactive).
            
        Raises:
            ValueError: "active_operator" does not contain ">" or "<"
        """
        if '>' in active_operator:
            return 1 if is_greater_than else 0
        elif '<' in active_operator:
            return 0 if is_greater_than else 1
        else:
            raise(ValueError, '"active_operator" does not contain ">" or "<"')
    
    @staticmethod
    @np.vectorize
    def bin_molar_qualifiers(value:Union[int, float, np.number], boundary:Union[int, float, np.number], relation:str,
                             active_operator:str) -> Union[int, float]:
        """Categorizes value as active (1), inactive (0), or conflict (np.nan) using qualifiers.

        Args:
            value (Union[int, float, np.number]): Number to evaluate threshold relationship.
            boundary (Union[int, float, np.number]): Binary decision boundary.
            relation (str): Value qualifier. Set to one of: =, <, >, ≤, or ≥.
            active_operator: Condition to set values to active based on threshold. Set to one of: <, >, <=, or >=.

        Returns:
            Union[int, float]: 1 (active), 0 (inactive), or np.nan for qualifier conflict / unrecognized relation.
        """
        
        greater_than_bound = None
        qualifier_conflict = False
        
        while not qualifier_conflict:
            if value == boundary:
                return Binarizer.__equals_edge_case(nan_greater_equal=True)
            
            elif relation in Binarizer.__less_than:
                if value < boundary:
                    greater_than_bound = False
                else:
                    qualifier_conflict = True
                    
            elif relation in Binarizer.__greater_than:
                if value > boundary:
                    greater_than_bound = True
                else:
                    qualifier_conflict = True
                
            elif relation == '=':
                greater_than_bound = True if value > boundary else False
                
            else:
                return np.nan  # Unrecognized relation
        
            return Binarizer.__active_operator_handler(greater_than_bound, active_operator)

        return np.nan  # Qualifier conflict
    
    @staticmethod
    @np.vectorize
    def bin_molar_no_qualifiers(value:Union[int, float, np.number], boundary:Union[int, float, np.number],
                                active_operator:str) -> Union[int, float]:
        """Categorizes value as active (1) or inactive (0).

        Args:
            value (Union[int, float, np.number]): Number to evaluate threshold relationship.
            boundary (Union[int, float, np.number]): Binary decision boundary.
            active_operator (str): Condition to set values to active based on threshold. Set to one of: <, >, <=, or >=.

        Returns:
            Union[int, float]: 1 (active), 0 (inactive), or np.nan special greater than equal edge case.
        """
        if value == boundary:
            return Binarizer.__equals_edge_case(nan_greater_equal=True)
        else:
            greater_than_bound = (value > boundary)
            return Binarizer.__active_operator_handler(greater_than_bound, active_operator)
    
    def binarize(self) -> np.array:
        """Main binarize function.

        Returns:
            np.array: Binarized values.
        """
        if self.__qualifiers:
            return self.bin_molar_qualifiers(value=self.__values, boundary=self.__boundary, relation=self.__qualifiers,
                                             active_operator=self.__active_operator)
        else:
            return self.bin_molar_no_qualifiers(value=self.__values, boundary=self.__boundary,
                                                active_operator=self.__active_operator)

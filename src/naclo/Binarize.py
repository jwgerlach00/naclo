from cmath import log10
from typing import List
import pandas as pd
from copy import copy
import numpy as np

import naclo
import stse


class Binarize:
    def __init__(self, df:pd.DataFrame, params:dict, options:dict) -> None:
        self.original_df = df.copy()
        self.df = df.copy()
        
        self.__structure_col = params['structure_col']
        self.__structure_type = params['structure_type']
        self.__target_col = params['target_col']  # Set to standard_value col if doing unit conversion
        self.__decision_boundary = params['decision_boundary']
        
        self.__options = copy(options)
        
        # Drop NA structures and targets
        self.df = stse.dataframes.convert_to_nan(self.df)
        self.df.dropna(subset=[self.__structure_col], inplace=True)
        self.df.dropna(subset=[self.__target_col], inplace=True)
    
    def __mol_weights(self) -> List[float]:
        if self.__structure_type == 'smiles':
            mols = naclo.smiles_2_mols(self.df[self.__structure_col])
        elif self.__structure_type == 'mol':
            mols = self.df[self.__structure_col]
        return naclo.mol_stats.mol_weights(mols)
    
    def __append_unit_convert(self, molar_values) -> None:
        output_units = self.__options['convert_units']['output_units']
        if output_units == 'neg_log_molar':
            self.df[f'neg_log_molar{self.__target_col}'] = [-log10(m) for m in molar_values]
        elif output_units == 'molar':
            self.df[f'molar_{self.__target_col}'] = molar_values
        
    def convert_units(self) -> pd.DataFrame:
        mws = self.__mol_weights()
        unit_converter = naclo.UnitConverter(values=self.df[self.__target_col],  # target_col == standard_val
                                             units=self.df[self.__options['convert_units']['units_col']],
                                             mol_weights=mws)
        return unit_converter.to_molar()
        
    def binarize(self, values) -> np.array:
        if self.__options['qualifiers']['run']:
            col = self.__options['qualifiers']['qualifier_col']
            self.df.dropna(subset=[col], inplace=True)  # Drop NA
            qualifiers = self.df[col].tolist()
        else:
            qualifiers = None
            
        molar_binarizer = naclo.MolarBinarizer(molar_vals=values,
                                               molar_boundary=self.__decision_boundary,
                                               molar_qualifiers=qualifiers,
                                               active_boundary_cond=self.__options['active_boundary_cond'])
        return molar_binarizer.binarize()
    
    def main(self) -> pd.DataFrame:
        if self.__options['convert_units']['run']:
            molar_values = self.convert_units().tolist()
            self.df[f'binarized_{self.__target_col}'] = self.binarize()
            self.__append_unit_convert()
        else:
            self.df[f'binarized_{self.__target_col}'] = self.binarize(self.df[self.__target_col].tolist())
            
        return self.df

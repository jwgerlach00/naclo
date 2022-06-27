from distutils import errors
from logging import warning
from unittest import TextTestRunner
import naclo
import stse
from rdkit import Chem
import json
import pandas as pd
import warnings


class Bleach:
    def __init__(self, df:pd.DataFrame, params:dict, options:dict) -> None:  # *
        self.__recognized_structures = ['smiles', 'mol']
        self.__default_cols = {
            'smiles': 'SMILES',
            'mol': 'ROMol',
            'inchi_key': 'InchiKey',
            'mw': 'MW'
        }
        
        # Save input data
        self.original_df = df.copy()
        self.df = df.copy()
        
        # Load file parameters
        self.structure_col = params['structure_col']
        self.structure_type = params['structure_type']
        self.__structure_checker()
        self.target_col = params['target_col']
        
        self.mol_col = None
        self.smiles_col = None
        self.__set_structure_cols()  # Assign mol and SMILES cols using input + defaults
        self.inchi_key_col = None
        
        # Load user options
        self.mol_settings = options['molecule_settings']
        self.file_settings = options['file_settings']
        
        self.__recognized_options_checker()
        
    def __build_mols(self):
        self.df = naclo.dataframes.df_smiles_2_mols(self.df, self.smiles_col, self.mol_col)
        
    def __build_smiles(self):
        self.df = naclo.dataframes.df_mols_2_smiles(self.df, self.mol_col, self.smiles_col)

    def __structure_checker(self) -> None:  # *
        """Checks for errors related to initialized parameters.

        Raises:
            ValueError: NO_STRUCTURE_COLUMN
            ValueError: STRUCTURE_COLUMN_NOT_FOUND
            ValueError: INVALID_STRUCTURE_TYPE
        """
        if not self.structure_col:
            raise ValueError('NO_STRUCTURE_COLUMN', 'Must specify the name of the structure column in params')
            
        if self.structure_col not in self.df.columns:
            raise ValueError('STRUCTURE_COLUMN_NOT_FOUND', f'The structure column: "{self.structure_col}"" is not \
                present in the data: "{list(self.df.columns)}"')
            
        if self.structure_type not in self.__recognized_structures:
            raise ValueError('INVALID_STRUCTURE_TYPE', f'Structure type: "{self.structure_type}"" is not one of: \
                {self.__recognized_structures}')

    def __recognized_options_checker(self) -> None:  # *
        """Checks for errors related to unrecognized options.

        Raises:
            ValueError: BAD_OPTION(S)
        """
        with open('naclo/assets/recognized_bleach_options.json') as f:
            recognized_options = json.load(f)
            
        input = stse.dictionaries.branches({
            'molecule_settings': self.mol_settings,
            'file_settings': self.file_settings
        })
        recognized = stse.dictionaries.branches(recognized_options)
        
        errors = {}
        for key, value in recognized.items():
            if isinstance(value, list):
                if not input[key] in recognized[key]:
                    errors[f'BAD_OPTION{key.upper()}'] = f'"{input[key]}" is not an accepted value for "{key}", set \
                        to one of: "{recognized[key]}"'
            else:
                if not type(input[key]) == type(recognized[key]):
                    errors[f'BAD_OPTION{key.upper()}'] = f'{type(input[key])} is not an accepted type for {key}, \
                        input a {type(recognized[key])}'
    
        if errors:
            raise ValueError(errors)
        
    def __drop_na_structures(self) -> None:  # *
        """Drops NA along a feature column depending on whether Mols or SMILES were inputted.
        """
        self.df.dropna(subset=[self.structure_col], inplace=True)
        if not len(self.df):
            warnings.warn('ALL_NA_STRUCTURES: All structures in specified column were NA, all rows dropped',
                          RuntimeWarning)
    
    def __drop_na_targets(self) -> None:  # *
        run_na_targets = self.file_settings['remove_na_targets']['run']

        if self.target_col and run_na_targets and len(self.df):  # If run and TARGET COLUMN DECLARED
            self.df.dropna(subset=[self.target_col], inplace=True)
            if not len(self.df):
                warnings.warn('ALL_NA_TARGETS: All targets in specified column were NA, all rows dropped',
                              RuntimeWarning)

        elif run_na_targets:
            warnings.warn('NA_TARGETS: options.file_settings.remove_na_targets was set to run but no activity column \
                was specified', RuntimeWarning)

    def __remove_fragments(self):
        """Removes salts if specified, then removes other fragments by appropriate method if specified.
        """
        option = self.mol_settings['remove_fragments']
        
        # Remove salts
        if option['salts']:
            self.df[self.mol_col] = naclo.fragments.remove_salts(self.df[self.mol_col], salts='[{0}]'.format(
                option['salts'].replace(' ', '')))
            self.__build_smiles()
            self.df = stse.dataframes.convert_to_nan(self.df, na=[''])  # Convert bc NA is just empty string
            self.df.dropna(subset=[self.smiles_col], inplace=True)  # Drop NA bc may include molecule that is ONLY salts
        
        # Filter
        if option['filter_method'] and option['filter_method'] != 'none':
            self.df[self.smiles_col] = self.df[self.smiles_col].apply(
                self.__filter_fragments_factory(option['filter_method']))
            self.__build_mols()
    
    def __neutralize_charges(self):
        """Neutralizes Mols and rebuilds SMILES.
        """
        self.df[self.mol_col] = naclo.neutralize.neutralize_charges(self.df[self.mol_col])
        self.__build_smiles
        
    def __set_structure_cols(self) -> None:
        self.mol_col = self.structure_col if self.structure_type == 'mol' else self.__default_cols['mol']
        self.smiles_col = self.structure_col if self.structure_type == 'smiles' else self.__default_cols['smiles']
        
    

    
    def main(self) -> pd.DataFrame:
        """Main bleach loop.

        Returns:
            pandas DataFrame: Cleaned df
        """
        self.drop_na()  # Before init_structure bc need NA
        self.init_structure_compute()
        self.mol_cleanup()  # Clean Mols and SMILES
        self.handle_duplicates()  # Drop/average/keep duplicates
        # self.append_columns()  # Add or remove columns from final output
        # self.remove_header_chars()  # Remove characters from headers
    
        return self.df
        
    # Step 1
    def drop_na(self) -> None:  # *
        """Converts blanks to NA, removes entire NA columns, drops NA Mols or SMILES, handles NA targets.
        """
        
        # Convert all df blanks and 'none' to NA
        self.df = stse.dataframes.convert_to_nan(self.df)

        # Drop rows
        self.__drop_na_structures()
        self.__drop_na_targets()
        
        # Drop cols
        self.df = stse.dataframes.remove_nan_cols(self.df)  # After dropping rows because columns may BECOME empty
    
    # Step 2 
    def init_structure_compute(self) -> None:  # *
        """Computes Mols if the input doesn't contain a Mol column already, refreshes SMILES by regenerating from Mols.
        """
        
        # BUILDING WILL ALSO DELETE NA ROWS FROM DF
        if self.structure_type == 'mol':
            self.__build_smiles()
            # Rebuilding Mols not necessary
            
        elif self.structure_type == 'smiles':
            self.__build_mols()
            # Canonicalize SMILES
            self.__build_smiles()  # Canonicalize SMILES
        
    # Step 3
    def mol_cleanup(self):  # *
        """Cleans Mols and SMILES.
        """
        
        # Step 1: Deal with fragments (includes salt step -- may include a molecule that is ONLY salts)
        # --> this is before dropping NA
        self.__remove_fragments()
        
        # Step 2: Neutralize mols
        if self.mol_settings['neutralize_charges']['run']:
            self.__neutralize_charges()
    
    def __compute_inchi_keys(self):
        self.inchi_key_col = self.__default_cols['inchi_key']
        self.df = naclo.dataframes.df_mols_2_inchi_keys(self.df, self.mol_col, self.inchi_key_col)
    
    # Step 4
    def handle_duplicates(self):  # *
        """Averages, removes, or keeps duplicates. ONLY BY INCHI KEY FOR NOW.

        Args:
            df (pandas DataFrame): Data to transform
        """
        self.__compute_inchi_keys()
        
        dup = self.file_settings['duplicate_compounds']
        
        if dup['selected'] == 'average' and self.target_col:
            self.df = stse.duplicates.average(self.df, subsets=[self.inchi_key_col], average_by=self.target_col)
        elif dup['selected'] == 'remove' or (dup['selected'] == 'average' and not self.target_col):
            self.df = stse.duplicates.remove(self.df, subsets=[self.inchi_key_col])
    
    # Step 6
    def append_columns(self):
        """Drops and adds columns depending on what the user wants returned.

        Args:
            df (pandas DataFrame): Data to transform
        """
        self.__drop_columns()
        self.__add_columns()
    
    # Step 7
    def remove_header_chars(self):
        """Removes any chars listed in a string of chars from the df column headers.
        """
        chars = self.file_settings['remove_header_chars']['chars']
        self.df = stse.dataframes.remove_header_chars(self.df, chars)
    
    def __drop_columns(self):  # Note: everything already has a Molecule column added
        """Removes columns that the user does not want in the final output.
        """
        option = self.file_settings['append_columns']
        
        # Drop added columns from built if not requested
        if not option['MolFile']:
            self.df.drop(self.mol_col, inplace=True)
        if not option['InchiKey']:
            self.df.drop(self.inchi_keys_col, inplace=True)
        if not option['SMILES']:
            self.df.drop(self.smiles_col, inplace=True)

    def __add_columns(self) -> None:
        """Add columns that the user wants in the final output.

        Args:
            df (pandas DataFrame): Data to transform

        Returns:
            pandas DataFrame: Transformed data
        """
        option = self.file_settings['append_columns']
        
        # Add MW column
        if option['MW']:
            self.df.assign(MW = naclo.mol_weights(self.df[self.mol_col]))
    
    @staticmethod
    def __filter_fragments_factory(method):
        if method == 'carbon_count':
            return naclo.fragments.carbon_count
        elif method == 'mw':
            return naclo.fragments.mw
        elif method == 'atom_count':
            return naclo.fragments.atom_count
        else:
            raise ValueError('Filter method is not allowed')

from distutils import errors
from logging import warning
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
            'inchi_key': 'InchiKey'
        }
        
        # Reference point to ensure none are lost (saved after blanks are dropped)
        self.__reference_mols = []
        self.__reference_smiles = []
        
        # Modified molecular properties
        self.__smiles = []
        self.__mols = []
        self.__inchi_keys = []
        
        # Save input data
        self.original_df = df.copy()
        self.df = df.copy()
        
        # Load file parameters
        self.structure_col = params['structure_col']
        self.structure_type = params['structure_type']
        self.target_col = params['target_col']
        
        # Load user options
        self.mol_settings = options['molecule_settings']
        self.file_settings = options['file_settings']
        
        self.__init_error_checker()
        self.__recognized_options_checker()
            
    def __init_error_checker(self) -> None:  # *
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
    
    def main(self) -> pd.DataFrame:
        """Main bleach loop.

        Returns:
            pandas DataFrame: Cleaned df
        """
            
        self.drop_na()  # Blanks to NA, NA columns, NA features (Mols/SMILES), NA targets
        self.compute_mols()  # Create and sanitize, refresh SMILES
        
        self.__set_reference_point()  # Reference point after dropping NAs and computing
        
        self.mol_cleanup()  # Clean Mols and SMILES
        
        # Ensure that cleaning dropped no Mols or SMILES
        if not self.__no_mol_drop():
            self.errors['STRUCTURES_DROPPED_AT_CLEAN'] = 'Cleaning dropped Mols of SMILES'
        
        built_df = self.build_df()  # Build dataframe from cleaned Mols, cleaned SMILES, and computed Inchi keys
        self.handle_duplicates(built_df)  # Drop/average/keep duplicates
        self.append_columns()  # Add or remove columns from final output
        self.remove_header_chars()  # Remove characters from headers
    
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
    def compute_mols(self):
        """Computes Mols if the input doesn't contain a Mol column already, refreshes SMILES by regenerating from Mols.
        """
        try:
            self.__mols = self.df[self.mol_col]  # Test if mol column already exists
            assert all([isinstance(x, Chem.rdchem.Mol) for x in self.__mols])  # Ensure the entire column is Mol
        except KeyError:
            self.df[self.smiles_col].dropna(inplace=True)
            self.df = naclo.dataframes.df_smiles_2_mols(self.df, self.smiles_col, self.mol_col)
            
            self.__mols = self.df[self.mol_col].tolist()
        
        self.__smiles = naclo.mols_2_smiles(self.__mols)
        
    # Step 3
    def mol_cleanup(self):
        """Cleans Mols and SMILES.
        """
        
        # Step 1: Deal with fragments (includes salt step -- may include a molecule that is ONLY salts)
        # --> this is before dropping NA
        self.__remove_fragments()
        
        # Step 2: Neutralize mols
        if self.mol_settings['neutralize_charges']['run']:
            self.__neutralize_charges()
    
    # Step 4
    def build_df(self):
        """Builds df prior to handling of duplicates.

        Returns:
            pandas DataFrame: self.df with Mols, SMILES, and Inchi Keys added
        """
        
        out = self.df.copy()
        out[self.mol_col if self.mol_col else self.__default_cols['mol']] = self.__mols
        out[self.smiles_col if self.smiles_col else self.__default_cols['smiles']] = self.__smiles
        out[self.inchi_keys_col if self.inchi_keys_col else self.__default_cols['inchi_key']] = \
            naclo.mols_2_inchi_keys(self.__mols)
        return out
    
    # Step 5
    def handle_duplicates(self, df:pd.DataFrame):
        """Averages, removes, or keeps duplicates.

        Args:
            df (pandas DataFrame): Data to transform
        """
        
        dup = self.file_settings['duplicate_compounds']
        
        if dup['selected'] == 'average' and self.target_col:
            df = stse.duplicates.average(df, subsets=[self.inchi_keys_col], average_by=self.target_col)
        elif dup['selected'] == 'remove' or (dup['selected'] == 'average' and not self.target_col):
            df = stse.duplicates.remove(df, subsets=[self.inchi_keys_col])
        
        self.df = df
    
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
        
    def __set_reference_point(self):
        """Sets reference Mols and SMILES.
        """
        self.__reference_smiles = self.__smiles
        self.__reference_mols = self.__mols
    
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
            
    def __remove_fragments(self):
        """Removes salts if specified, then removes other fragments by appropriate method if specified.
        """
        option = self.mol_settings['remove_fragments']
        
        # Remove salts (separate option -- bc may include a molecule that is ONLY salts --> NA drop is later)
        if option['salts']:
            self.__mols = naclo.fragments.remove_salts(self.__mols, salts='[{0}]'.format(option['salts'].replace(' ', '')))
        
        # Controlled by method of removing fragments
        if not option['filter_method']:  # Break if no method
            return
        elif option['filter_method'] == 'carbon_count':
            self.__smiles = [naclo.fragments.carbon_count(s) for s in self.__smiles]
        elif option['filter_method'] == 'mw':
            self.__smiles = [naclo.fragments.mw(s) for s in self.__smiles]
        elif option['filter_method'] == 'atom_count':
            self.__smiles = [naclo.fragments.atom_count(s) for s in self.__smiles]
        
        self.__mols = naclo.smiles_2_mols(self.__smiles)  # Build Mols from new SMILES
        
    def __neutralize_charges(self):
        """Neutralizes Mols and rebuilds SMILES.
        """
        self.__mols = naclo.neutralize.neutralize_charges(self.__mols)
        self.__smiles = naclo.mols_2_smiles(self.__mols)  # Build SMILES from new Mols
    
    def __no_mol_drop(self):
        """Evaluates truth of equality between lengths of Mols, SMILES, reference Mols, reference SMILES.

        Returns:
            bool: Truth of equality
        """
        return len(self.__mols) == len(self.__smiles) == len(self.__reference_mols) == len(self.__reference_smiles)


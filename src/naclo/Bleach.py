import json
import pandas as pd
import warnings
from typing import Callable

# sourced from github.com/jwgerlach00
import naclo
import stse
from naclo.__asset_loader import recognized_bleach_options as recognized_options


class Bleach:
    def __init__(self, df:pd.DataFrame, params:dict, options:dict) -> None:  # *
        # Load user options
        self.mol_settings = options['molecule_settings']
        self.file_settings = options['file_settings']
        self.__recognized_options_checker()

        self.__recognized_structures = ['smiles', 'mol']
        self.__default_cols = {
            'smiles': 'SMILES',
            'mol': 'ROMol',
            'inchi_key': 'InchiKey',
            'mw': 'MW'
        }

        # Save user input data
        self.original_df = df.copy()
        self.df = df.copy()

        # Load file parameters
        self.structure_col = params['structure_col']
        self.structure_type = params['structure_type']
        self.target_col = params['target_col']
        self.__param_checker()

        self.mol_col = None
        self.smiles_col = None
        self.__set_structure_cols()  # Assign mol and SMILES cols using input + defaults
        self.inchi_key_col = None


# -------------------------------------------------- ERROR CHECKING -------------------------------------------------- #
    def __recognized_options_checker(self) -> None:  # *
        """Checks for errors related to unrecognized option assignments.

        Raises:
            ValueError: BAD_OPTION(S)
        """
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

    def __param_checker(self) -> None:  # *
        """Checks for errors related to declared parameters.

        Raises:
            ValueError: NO_STRUCTURE_COLUMN
            ValueError: STRUCTURE_COLUMN_NOT_FOUND
            ValueError: INVALID_STRUCTURE_TYPE
            ValueError: TARGET_COLUMN_NOT_FOUND
        """
        if not self.structure_col:
            raise ValueError('NO_STRUCTURE_COLUMN', 'Must specify the name of the structure column in params')

        if self.structure_col not in self.df.columns:
            raise ValueError('STRUCTURE_COLUMN_NOT_FOUND', f'The structure column: "{self.structure_col}"" is not \
                present in the data: "{list(self.df.columns)}"')

        if self.structure_type not in self.__recognized_structures:
            raise ValueError('INVALID_STRUCTURE_TYPE', f'Structure type: "{self.structure_type}"" is not one of: \
                {self.__recognized_structures}')

        if self.target_col and self.target_col not in self.df.columns:
            raise ValueError('TARGET_COLUMN_NOT_FOUND', f'The target column: "{self.target_col}"" is not present in \
                the data: "{list(self.df.columns)}"')


# -------------------------------------------------- STATIC METHODS -------------------------------------------------- #
    @staticmethod
    def __filter_fragments_factory(filter:str) -> Callable:
        """Returns a callable SMILES fragment filter function using a key.

        Args:
            filter (str): Function key.

        Raises:
            ValueError

        Returns:
            Callable: Filter function.
        """
        if filter == 'carbon_count':
            return naclo.fragments.carbon_count
        elif filter == 'mw':
            return naclo.fragments.mw
        elif filter == 'atom_count':
            return naclo.fragments.atom_count
        else:
            raise ValueError('Filter method is not allowed')


# -------------------------------------------------- PRIVATE METHODS ------------------------------------------------- #
    def __set_structure_cols(self) -> None:
        """Sets Mol and SMILES columns using declared structure type."""
        self.mol_col = self.structure_col if self.structure_type == 'mol' else self.__default_cols['mol']
        self.smiles_col = self.structure_col if self.structure_type == 'smiles' else self.__default_cols['smiles']

    def __drop_na_structures(self) -> None:
        """Drops NA along declared structure column."""
        self.df.dropna(subset=[self.structure_col], inplace=True)
        if not len(self.df):
            warnings.warn('ALL_NA_STRUCTURES: All structures in specified column were NA, all rows dropped',
                          RuntimeWarning)

    def __drop_na_targets(self) -> None:
        """Drops NA along declared target column"""
        run_na_targets = self.file_settings['remove_na_targets']['run']

        if self.target_col and run_na_targets and len(self.df):  # If run and TARGET COLUMN DECLARED
            self.df.dropna(subset=[self.target_col], inplace=True)
            if not len(self.df):
                warnings.warn('ALL_NA_TARGETS: All targets in specified column were NA, all rows dropped',
                              RuntimeWarning)

        elif run_na_targets:  # If run but not declared target
            warnings.warn('NA_TARGETS: options.file_settings.remove_na_targets was set to run but no activity column \
                was specified', RuntimeWarning)

    def __build_smiles(self) -> None:
        """Creates a SMILES column in the dataset using dataset MolFile column. DROPS NA."""
        self.df = naclo.dataframes.df_mols_2_smiles(self.df, self.mol_col, self.smiles_col)

    def __build_mols(self) -> None:
        """Creates MolFile column in the dataset using dataset SMILES column. DROPS NA."""
        self.df = naclo.dataframes.df_smiles_2_mols(self.df, self.smiles_col, self.mol_col)

    def __remove_fragments(self) -> None:
        """Removes salts if specified. Drops NA as a result of salt removal. Filters out other fragments by specified
        method."""
        option = self.mol_settings['remove_fragments']

        # Remove salts
        if option['salts']:
            self.df[self.mol_col] = naclo.fragments.remove_salts(self.df[self.mol_col], salts='[{0}]'.format(
                option['salts'].replace(' ', '')))
            self.__build_smiles()

            # Drop NA (blank string after salts)
            self.df = stse.dataframes.convert_to_nan(self.df, na=[''])  # Convert bc NA is just empty string
            self.df.dropna(subset=[self.smiles_col], inplace=True)  # Drop NA bc may include molecule that is ONLY salts

        # Filter
        if option['filter_method'] and option['filter_method'] != 'none':
            self.df[self.smiles_col] = self.df[self.smiles_col].apply(
                self.__filter_fragments_factory(option['filter_method']))
            self.__build_mols()

    def __neutralize_charges(self) -> None:
        """Neutralizes Mols. Rebuilds SMILES."""
        self.df[self.mol_col] = naclo.neutralize.neutralize_charges(self.df[self.mol_col])
        self.__build_smiles

    def __compute_inchi_keys(self) -> None:
        """Declares inchi key column name using default. Appends inchi keys to dataset."""
        self.inchi_key_col = self.__default_cols['inchi_key']
        self.df = naclo.dataframes.df_mols_2_inchi_keys(self.df, self.mol_col, self.inchi_key_col)

    def __drop_columns(self) -> None:
        """Removes columns that the user does not want in the final output."""
        option = self.file_settings['append_columns']

        # Drop added columns from built if not requested
        if not option['mol']:
            self.df.drop(self.mol_col, inplace=True, axis=1)
        if not option['inchi_key']:
            self.df.drop(self.inchi_key_col, inplace=True, axis=1)
        if not option['smiles']:
            self.df.drop(self.smiles_col, inplace=True, axis=1)

    def __add_columns(self) -> None:
        """Add columns that the user wants in the final output."""
        option = self.file_settings['append_columns']
        
        # Add MW column
        if option['mw']:
            self.df = self.df.assign(MW = naclo.mol_weights(self.df[self.mol_col]))


# ------------------------------------------------- PUBLIC FUNCTIONS ------------------------------------------------- #
    # Step 1
    def drop_na(self) -> None:  # *
        """Converts blanks to NA. Drops NA Mols or SMILES. Handles NA targets. Removes entire NA columns"""

        # Convert all df blanks and 'none' to NA
        self.df = stse.dataframes.convert_to_nan(self.df)

        # Drop rows
        self.__drop_na_structures()
        self.__drop_na_targets()

        # Drop cols
        self.df = stse.dataframes.remove_nan_cols(self.df)  # After dropping rows because columns may BECOME empty

    # Step 2 
    def init_structure_compute(self) -> None:  # *
        """Builds (or rebuilds from Mols) SMILES. Builds Mols if not present in dataset."""
        if self.structure_type == 'mol':
            self.__build_smiles()
            # Rebuilding Mols not necessary

        elif self.structure_type == 'smiles':
            self.__build_mols()
            self.__build_smiles()  # Canonicalize SMILES

    # Step 3
    def mol_cleanup(self):  # *
        """Cleans Mols and SMILES."""

        # Step 1: Deal with fragments (includes salt step -- may include a molecule that is ONLY salts (NA dropped))
        self.__remove_fragments()

        # Step 2: Neutralize mols
        if self.mol_settings['neutralize_charges']['run']:
            self.__neutralize_charges()

    # Step 4
    def handle_duplicates(self):  # *
        """Computes inchi keys. Averages, removes, or keeps duplicates. ONLY BY INCHI KEY FOR NOW."""
        self.__compute_inchi_keys()

        dup = self.file_settings['duplicate_compounds']

        if dup['selected'] == 'average' and self.target_col:
            self.df = stse.duplicates.average(self.df, subsets=[self.inchi_key_col], average_by=self.target_col)
        elif dup['selected'] == 'remove' or (dup['selected'] == 'average' and not self.target_col):
            self.df = stse.duplicates.remove(self.df, subsets=[self.inchi_key_col])

    # Step 5
    def append_columns(self) -> None:  # *
        """Drops and adds columns depending on what the user wants returned.

        Args:
            df (pandas DataFrame): Data to transform
        """
        self.__drop_columns()
        self.__add_columns()

    # Step 6
    def remove_header_chars(self) -> None:  # *
        """Removes any chars listed in a string of chars from the df column headers.
        """
        chars = self.file_settings['remove_header_chars']['chars']
        self.df = stse.dataframes.remove_header_chars(self.df, chars)


# ----------------------------------------------------- MAIN LOOP ---------------------------------------------------- #
    def main(self) -> pd.DataFrame:
        """Main bleach loop.

        Returns:
            pandas DataFrame: Cleaned df
        """
        self.drop_na()  # Before init_structure bc need NA
        self.init_structure_compute()
        self.mol_cleanup()
        self.handle_duplicates()  # Includes InchiKey generation
        self.append_columns()
        self.remove_header_chars()
        
        return self.df

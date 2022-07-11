import warnings


def drop_na_structures(self) -> None:
        """Drops NA along declared structure column."""
        self.df.dropna(subset=[self.structure_col], inplace=True)
        if not len(self.df):
            warnings.warn('ALL_NA_STRUCTURES: All structures in specified column were NA, all rows dropped',
                          RuntimeWarning)

def drop_na_targets(self) -> None:
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
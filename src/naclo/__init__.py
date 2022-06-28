__doc__ == """
**naclo** is a Python cleaning toolset for small molecule drug discovery datasets.
"""


from turtle import clear
from naclo.mol_stats import *
from naclo.mol_conversion import *
from naclo import database
from naclo import dataframes
from naclo import fragments
from naclo import neutralize
from naclo import rdpickle
from naclo import visualization
from naclo.Bleach import Bleach

import json

with open('naclo/assets/default_options.json') as f:
    bleach_default_options = json.load(f)
    
with open('naclo/assets/default_params.json') as f:

    bleach_default_params = json.load(f)

from importlib import resources
import json


with resources.open_text('naclo.assets', 'default_params.json') as f:
    bleach_default_params = json.load(f)

with resources.open_text('naclo.assets', 'default_options.json') as f:
    bleach_default_options = json.load(f)

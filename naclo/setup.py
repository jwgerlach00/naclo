from setuptools import setup, find_packages

VERSION = '0.0.2'
DESCRIPTION = 'Collection of compound data cleaning operations.'
LONG_DESCRIPTION = 'Collection of compound data cleaning operations.'

setup(
    name='naclo',
    version=VERSION,
    author='Jacob Gerlach',
    author_email='jwgerlach00@gmail.com',
    description=DESCRIPTION,
    long_description=LONG_DESCRIPTION,
    packages=find_packages(),
    install_requires=[],
    keywords=['python', 'molecule', 'smiles', 'inchi', 'rdkit', 'chem'],
    classifiers= [
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Education",
        "Programming Language :: Python :: 2",
        "Programming Language :: Python :: 3",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: Microsoft :: Windows",
    ]
)
from setuptools import setup, find_packages


setup(
    name='naclo',
    version='0.1.77',
    license='MIT',
    author='Jacob Gerlach',
    author_email='jwgerlach00@gmail.com',
    url='https://github.com/jwgerlach00/naclo',
    description='Cleaning toolset for small molecule drug discovery datasets',
    packages=find_packages('src'),
    package_dir={'': 'src'},
    package_data={'naclo': ['assets/default_params.json',
                            'assets/default_options.json',
                            'assets/recognized_bleach_options.json']},
    python_requires='>=3.6',
    install_requires=[
        'numpy',
        'pandas',
        'rdkit',
        'rdkit_pypi',
        'setuptools',
        'stse'
    ],
)

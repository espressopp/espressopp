from setuptools import setup, Extension

print(setuptools.find_packages(where="."))

setup(
    name='espressopp',
    version='@PROJECT_VERSION@',
    packages=setuptools.find_packages(where="."),
    package_dir={'': '.', 'espressopp': 'espressopp'},
    package_data={'': ['_espressopp.so']},
)

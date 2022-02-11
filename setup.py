from setuptools import setup, find_packages
import sys

sys.path.append('./PolyMolParser')

setup(name='PolyMolParser',
        version='2022.2.11',
        description='PolyMolParser',
        long_description="README",
        author='Kan Hatakeyama',
        license="MIT",
        packages = find_packages(),
    )
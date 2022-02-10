from setuptools import setup, find_packages
import sys

sys.path.append('./polysmiles')

setup(name='PolyMolParser',
        version='2022.2.10',
        description='PolyMolParser',
        long_description="README",
        author='Kan Hatakeyama',
        license=license,
        packages = find_packages(),
    )
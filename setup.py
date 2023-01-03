from setuptools import setup, find_packages

VERSION = '0.0.0'
DESCRIPTION = 'Parser for Quantum Espresso'
LONG_DESCRIPTION = 'My python package for Quantum Espresso, made in ANL'

# Setting up
setup(
        name="parserQE",
        version=VERSION,
        author="RayHEsparza",
        author_email="rayhe88@gmail.com",
        description=DESCRIPTION,
        long_description=LONG_DESCRIPTION,
        packages=find_packages(),
        install_requires=[],
        keywords=['QuantumEspresso','Parser'],
        classifiers=[
            "Programming Language :: Python :: 3",
            "Operating System :: Linux"
            ] 
    )

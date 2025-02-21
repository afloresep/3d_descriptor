from setuptools import setup, find_packages

setup(
    name="three_desc",
    version="0.1.0",
    description="A package for computing 3D descriptors from molecular structures.",
    author="Alejandro Flores",
    packages=find_packages(),  # This will automatically find the "my_3d_descriptors" package folder
    install_requires=[
        "rdkit",
        "numpy",
        "pandas",
        "matplotlib",
        "tqdm", 
        "pyarrow",
    ],
    extras_require={
        "dev": [
            "pytest",
            "flake8",
            "black",
            "mypy",
            "pytest-cov",
        ],
    },
    python_requires=">=3.9",
)
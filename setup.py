from setuptools import setup, find_packages

setup(
    name="wsigma",
    version="1.1.0",
    description="evaluation of sound speed at given pressure, temperature and chemical composition",
    author="Ivan Selyakov, Bertram Boehrer",
    packages=find_packages(),
    python_requires='>=3.8',
)
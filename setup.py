```python
from setuptools import setup, find_packages

setup(
    name="manual-visium-loader",
    version="0.1",
    description="Manual loader for 10x Visium spatial transcriptomics data",
    author="Your Name",
    author_email="your@email.com",
    packages=find_packages(),
    install_requires=[
        "scanpy",
        "anndata",
        "pandas",
        "numpy",
        "scipy",
        "Pillow"
    ],
    python_requires=">=3.7",
)

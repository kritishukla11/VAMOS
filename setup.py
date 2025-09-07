from setuptools import setup, find_packages

setup(
    name="VAMOS",
    version="0.1.0",
    packages=find_packages(),
    install_requires=[
        "numpy",
        "pandas",
        "scikit-learn",
        "matplotlib",
        "seaborn",
        "scipy",
        "jupyter",
        "notebook",
        "nbconvert",
        "biopython"
    ],
    author="Your Name",
    author_email="your.email@example.com",
    description="Variant Mapping and Oncogenic Signatures (VAMOS) analysis toolkit",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/yourusername/VAMOS",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
                python_requires=">=3.8",
extras_require={
    "dev": ["pytest", "flake8", "black", "isort"]
},

    entry_points={
        "console_scripts": [
            "vamos=vamos.cli:main"
        ]
    }
)


import os
from setuptools import setup, find_packages

# Setup base directory
base_dir = os.path.abspath(os.path.dirname(__file__))

# Versioning
VERSION = '1.0.0'
version_path = os.path.join(base_dir, 'alignit', 'version.py')
with open(version_path, 'w') as version_file:
    version_content = "\n".join([
        "# THIS FILE IS GENERATED FROM SETUP.PY",
        f"version = '{VERSION}'",
        "__version__ = version"
    ])
    version_file.write(version_content)

# Package setup
setup(
    name='alignit',
    version=VERSION,
    description='Alignment tool for genomic data',  # Updated to be more descriptive
    long_description=open(os.path.join(base_dir, 'README.md')).read(),  
    long_description_content_type='text/markdown',  
    author='Tiffany Zhang, Jenny Nguyen, Corey Nguyen',
    author_email='con002@ucsd.edu',
    url='https://github.com/cnidyllic/alignit',  # Assuming a GitHub repo URL
    packages=find_packages(),
    include_package_data=True,
    entry_points={
        "console_scripts": [
            "alignit=alignit.main:main"  # Updated to reflect new tool name
        ],
    },
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: MIT License',  # Assuming an MIT license
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
    ],
    python_requires='>=3.7',  # Specifies the Python version requirement
)

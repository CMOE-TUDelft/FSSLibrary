from setuptools import setup, find_packages

setup(
    name='FSSLibrary',  # Replace with your package name
    version='0.1.1',
    description='Set of python library made specifically for Floating and Submerged Structures course.',
    author='Shagun Agarwal, Oriol Colomés',
    author_email='shagun.1995@gmail.com, j.o.colomesgene@tudelft.nl',
    url='https://github.com/CMOE-TUDelft/FSSLibrary',
    packages=find_packages(),  # Automatically finds all packages/modules
    install_requires=[
        # List dependencies here, e.g.:
        'numpy',
        'scipy',
        'typing'
    ],
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',  # Or your license
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.7',
)

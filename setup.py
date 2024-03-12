# include information about your package
from setuptools import setup, find_packages

setup(
    name='HTHubID',
    version='0.0.1',
    author='Chi-Yun Wu',
    author_email='yun820627@gmail.com',
    description='A brief description of your package',
    long_description=open('README.md').read(),
    url='https://yourpackagehomepage.com',
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    license='LICENSE',
    python_requires='>=3.6',
    install_requires=[
        'numpy >= 1.19.2',
        'gensim >= 4.3.2',
        'pandas',
        'h5py',
        'scanpy',
        'scipy',
        'matplotlib',
        'scikit-learn',
        'argparse',
        'seaborn',
    ]
)

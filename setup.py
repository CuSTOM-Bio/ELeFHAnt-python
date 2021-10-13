from setuptools import find_packages, setup
setup(
    name='elefhant',
    packages=find_packages(include=['elefhant']),
    version='1.0.0',
    description='A supervised machine learning approach for label harmonization and annotation of single cell RNA-seq data',
    author='Konrad Thorner',
    license='MIT',
    download_url = '',
    install_requires=['numpy==1.20.3','pandas==1.3.3','scanpy==1.8.1','anndata==0.7.6','bbknn==1.5.1','matplotlib==3.4.3','seaborn==0.11.2','gseapy==0.10.5'],
    setup_requires=[]
)
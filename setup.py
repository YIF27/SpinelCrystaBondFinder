from setuptools import setup, find_packages

setup(
    name='SpinelCrystalBondFinder',
    version='0.1',
    packages=find_packages(),
    install_requires=[
        # dependencies, e.g., 'numpy>=1.20.0'
    ],
    author='Ying Fang',
    author_email='yif27@pitt.edu',
    description='spinel crystal bond analysis',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    license='PITT',
    keywords='some keywords',
    url='https://github.com/YIF27/SpinelCrystaBondFinder'
)

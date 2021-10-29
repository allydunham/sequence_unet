from setuptools import setup, find_packages

setup(name='sequence_unet',
      version="1.0.0",
      description="Make protein predictions with Sequence UNET and train new models",
      url='https://github.com/allydunham/sequence_unet',
      author='Alistair Dunham',
      author_email='alistair.dunham@ebi.ac.uk',
      license='Apache 2.0',
      packages=find_packages(include="sequence_unet"),
	  include_package_data=False,
      install_requires=['numpy', 'pandas', 'biopython', 'tensorflow>=2.0', 'proteinnetpy>=0.5.2'],
      extras_require={},
      entry_points = {
        'console_scripts': ['sequence_unet=sequence_unet.scripts.make_preds:main',
                            'split_fasta=sequence_unet.scripts.split_fasta:main',
                            'filter_fasta=sequence_unet.scripts.filter_fasta:main']
      },
      zip_safe=True)

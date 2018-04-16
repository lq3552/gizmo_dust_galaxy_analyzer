from distutils.core import setup
from dust_gal_analyzer import __version__ as version

setup(name='DGalA',
		version=version,
		description = 'A tool to extract and analyze galaxies from GIZMO simulations',
		author = 'Qi Li',
		author_email=['lq3552@gmail.com'], 
		url='https://bitbucket.org/lq3552/dust_gal_analyzer',
		requires=['numpy(>=1.12.0)', 'matplotlib(>=2.0.0)','yt(>=3.4.1)','caesar'],
		package_dir={'dust_gal_analyzer': 'dust_gal_analyzer'}, # the present directory maps to src 
		packages = ['dust_gal_analyzer'],
	 )

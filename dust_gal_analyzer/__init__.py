__version__ = '1.4'
__all__ = ['DustyGalaxy','DustyGalaxyExtractor','constants']

# Import the classes and functions from the files of the same
# name. This allows users to access the classes and functions
# directly, as dust_gal_analyzer.dusty_galaxy, rather than 
#having to use the cumbersome notation 
#dust_gal_analyzer.dusty_galaxy.dusty_galaxy or similar.

from .DustyGalaxy import DustyGalaxy
from .DustyGalaxyExtractor import DustyGalaxyExtractor

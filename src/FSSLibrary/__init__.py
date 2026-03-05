import warnings as _warn
from FSSLibrary import Constants
from FSSLibrary import BeamMatrices
from FSSLibrary import LinearWave
from FSSLibrary import MoorLib
from FSSLibrary import FFTBasic
from FSSLibrary import WaveSpectrum

# Ignore future warnings
_warn.simplefilter(action="ignore", category=FutureWarning)

__version__ = "v0.2.0"

__copyright__ = """
Shagun Agarwal and Oriol Colomés."""

__license__ = "Revised BSD License"

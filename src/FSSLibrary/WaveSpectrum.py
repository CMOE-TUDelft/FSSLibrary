import numpy as np
from scipy import optimize
from FSSLibrary import LinearWave


class SpecStruct:
    
    def __init__(
        self, Hs: float, Tp: float, d: float, 
        omega: np.ndarray, S: np.ndarray, A: np.ndarray, 
        phRad: np.ndarray, *, k: np.ndarray|None=None
    ):
        """
        Wave spectrum structure containing spectral parameters
        
        ## Parameters:
        - Hs (float): Significant wave height [m]
        - Tp (float): Peak wave period [s]
        - d (float): Water depth [m]
        - omega (array): Angular frequency array [rad/s]
        - S (array): Power spectral density array [m^2*s]
        - A (array): Wave amplitude array [m]        
        - phRad (array): Phase array [rad]                
        - k (array | None, optional): Wave number array [rad/m]
        """
        
        N = len(omega)
        if (k is not None) and len(k) != N:
            raise ValueError("Input array k must have the same length as omega.")
        if not (len(S) == N and len(A) == N and len(phRad) == N):
            raise ValueError("Input arrays omega, S, A, and phRad must have the same length.")        
        
        self.Hs = Hs
        self.Tp = Tp
        self.d = d

        self.numOmega = N
        self.omega = np.array(omega)
        self.f = self.omega / (2 * np.pi)
        self.S = np.array(S)
        self.A = np.array(A)        
        self.phRad = np.array(phRad)        

        if k is not None:
            self.k = np.array(k)
            wvAll = [ 
                LinearWave.LinearWave2D(
                    self.d, 1/fi, 2*Ai, phi=phi, L=2*np.pi/ki, msg=False ) 
                for fi,Ai,phi,ki in zip(self.f, self.A, self.phRad, self.k) 
            ]        
            self.wvAll = wvAll
        else:
            wvAll = [ 
                LinearWave.LinearWave2D(
                    self.d, 1/fi, 2*Ai, phi=phi, msg=False ) 
                for fi,Ai,phi in zip(self.f, self.A, self.phRad) 
            ]        
            self.wvAll = wvAll
            self.k = np.array([ wv.k for wv in self.wvAll ])


    
    
    def __str__(self):
        return (
            f"WaveSpectrum.SpecStruct:\n"
            f"  Hs: {self.Hs}\n"
            f"  Tp: {self.Tp}\n"
            f"  d: {self.d}\n"
            f"  numOmega: {self.numOmega}\n"
            f"  omega: array of length {len(self.omega)}\n"
            f"  f: array of length {len(self.f)}\n"
            f"  S: array of length {len(self.S)}\n"
            f"  A: array of length {len(self.A)}\n"
            f"  k: array of length {len(self.k)}\n"
            f"  phRad: array of length {len(self.phRad)}\n"
            f"  wvAll: list of length {len(self.wvAll)}"
        )
    


    def waveElevation(self, t: float|np.ndarray, x: float) -> float|np.ndarray:
        """
        Calculate wave elevation at location x and time(s) t by summing all wave components.
        
        ## Parameters:
        - t (float or array): Time or time array [s]
        - x (float): Horizontal position [m]
        
        ## Returns:
        - float or array: Wave elevation [m]
        """
        
        if isinstance(t, float):
            eta = 0.0
        else:
            eta = np.zeros_like(t, dtype=float)
        
        for wv in self.wvAll:
            eta += wv.waveElevation(t, x)

        return eta

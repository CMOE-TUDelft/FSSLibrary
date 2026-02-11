import numpy as np
from scipy import optimize
from FSSLibrary.Constants import g, rhoW

def dispersionRelSol(sol, *data):
    g, d, T = data
    L = sol
    return (L-(g/2/np.pi*T*T*np.tanh(2*np.pi/L*d)))

def getWaveLen(g, d, T):
    """
    Calculate the wave length using the dispersion relation.

    ## Parameters:
    - g (float): Gravitational acceleration
    - d (float): Water depth
    - T (float): Wave period
    
    ## Returns:
    - float: Wave length
    
    ## Notes: 
    - Remember that `fsolve` is used to find the roots of a function, and it requires an initial guess. The initial guess is crucial for the convergence of the solution. In this case, we can use the deep water wave length as the initial guess, which is given by L = g/2/pi * T^2.
    - `fsolve` (to `optimize.brentq`) are both, ROOT FINDING METHODS. `fsolve` is a more general method that can find roots of any function, while `brentq` is specifically designed for finding roots of continuous functions. In this case, since we are dealing with a continuous function (the dispersion relation), `brentq` is more appropriate and efficient for finding the wave length.
    - For minimisation you can use `optimize.minimize_scalar` or `optimize.minimize`, but in this case, since we are interested in finding the wave length that satisfies the dispersion relation (i.e., where the function equals zero), root finding is the correct approach. Minimisation would not be appropriate here because we are not trying to find the minimum value of a function, but rather the specific value of L that makes the function equal to zero.
    - ROOT FINDING IS NOT THE SAME AS MINIMISATION. 
        - Root finding is about finding the input values that make a function equal to zero.
        - Minimisation is about finding the input values that make a function reach its minimum value.
    """
    data = (g, d, T)
    root = optimize.fsolve(dispersionRelSol, g/2/np.pi*T*T ,args=data)
    return(root[0])

class LinearWave2D:    
    
    def __init__(            
            self, d:float, T:float, H:float, *,
            phi:float=0.0, x0:float=0.0, L:float|None=None, msg:bool=True):
        """
        Initialize a LinearWave object.

        ## Parameters:
        - d (float): Water depth
        - T (float): Wave period
        - H (float): Wave height
        - phi (float, optional): Phase angle. Defaults to 0.0.
        - x0 (float, optional): Initial position. Defaults to 0.0.
        - L (float | None, optional): Wave length. If None, calculated from dispersion relation. Defaults to None.
        - msg (bool, optional): Print wave parameters. Defaults to True.
        """
        
        self.rhoW = rhoW
        self.g = g
        self.d = d
        self.x0 = x0
        self.phi = phi                
        self.T = T        
        self.H = H 
        
        if L is not None:
            self.L = L
        else:               
            self.L = getWaveLen(self.g, d, T)        

        self.k = 2*np.pi/self.L
        self.w = 2*np.pi/self.T

        if(msg):
            print('Wave-Length L = ', self.L)
            print('d/L = ', self.d / self.L)
        
        
    def wavePhase(self, t: float, x: float) -> float:
        return self.k*(x-self.x0) - self.w*t + self.phi
        
    
    def waveElevation(self, t: float, x: float) -> float:
        et = self.H / 2 * np.cos( self.wavePhase(t, x) )
        return et
    
    
    def pressureTotPoi(self, t, x, z):
        et = self.waveElevation(t, x)
        if(z>0):
            if(z>et):
                pdyn = 0.0*et
            else:
                pdyn = et-z
        else:
            if(z<et):
                Kp = np.cosh( self.k* (self.d + z) ) / np.cosh( self.k* self.d )
                pdyn = Kp*et - z
            else:
                pdyn=0.0*et
        return self.rhoW * self.g * pdyn
    
    def pressureTot(self, t, x, z):
        PdynRes = [self.pressureTotPoi(t,ix,iz) for ix, iz in zip(x,z)]        
        return PdynRes
    
    
    def pressureDynPoi(self, t, x, z):
        et = self.waveElevation(t, x)
        if(z>0):
            if(z>et):
                pdyn = 0.0*et
            else:
                pdyn = et-z
        else:
            if(z<et):
                Kp = np.cosh( self.k* (self.d + z) ) / np.cosh( self.k* self.d )
                pdyn = Kp*et                
            else:
                pdyn=0.0*et
        return self.rhoW * self.g * pdyn
    
    def pressureDyn(self, t, x, z):
        PdynRes = [self.pressureDynPoi(t,ix,iz) for ix, iz in zip(x,z)]        
        return PdynRes
    
    
    def particleVelPoi(self, t, x, zin):        
        et = self.waveElevation(t, x)   
        # # Wheeler stretching
        # z = self.d * ( self/d + zin ) / ( self/d + et ) - self.d     
        # Not applying Wheeler Stretching
        z = zin
        if(z>0):
            vx = 0.0 * et
            vz = 0.0 * et
        else:
            mag = self.H/2 * self.w
            vx =  mag * np.cosh( self.k*(self.d + z) ) / np.sinh( self.k * self.d )
            vx = vx * np.cos( self.wavePhase(t, x) )
            vz =  mag * np.sinh( self.k*(self.d + z) ) / np.sinh( self.k * self.d )
            vz = vz * np.sin( self.wavePhase(t, x) )
        return vx, vz
    
    
    def particleAccPoi(self, t, x, zin):
        et = self.waveElevation(t, x)
        # # Wheeler stretching
        # z = self.d * ( self/d + zin ) / ( self/d + et ) - self.d     
        # Not applying Wheeler Stretching
        z = zin
        if(z>0):
            vx = 0.0 * et
            vz = 0.0 * et
        else:
            mag = self.H/2 * self.w**2
            vx =  mag * np.cosh( self.k*(self.d + z) ) / np.sinh( self.k * self.d )
            vx = vx * np.sin( self.wavePhase(t, x) )
            vz =  mag * np.sinh( self.k*(self.d + z) ) / np.sinh( self.k * self.d )
            vz = -vz * np.cos( self.wavePhase(t, x) )
        return vx, vz

    
    def particleVelMax(self, x, zin):
        et = self.H  / 2
        # # Wheeler stretching
        # z = self.d * ( self/d + zin ) / ( self/d + et ) - self.d     
        # Not applying Wheeler Stretching
        z = zin
        mag = self.H/2 * self.w
        vx =  mag * np.cosh( self.k*(self.d + z) ) / np.sinh( self.k * self.d )            
        vz =  mag * np.sinh( self.k*(self.d + z) ) / np.sinh( self.k * self.d )
        if(z>0):
            vx = vx * 0.0
            vz = vz * 0.0
        return vx, vz


    def particleAccMax(self, x, zin):        
        et = self.H  / 2
        # # Wheeler stretching
        # z = self.d * ( self/d + zin ) / ( self/d + et ) - self.d     
        # Not applying Wheeler Stretching
        z = zin
        mag = self.H/2 * self.w**2
        vx =  mag * np.cosh( self.k*(self.d + z) ) / np.sinh( self.k * self.d )
        vz =  mag * np.sinh( self.k*(self.d + z) ) / np.sinh( self.k * self.d )
        if(z>0):
            vx = vx * 0.0
            vz = vz * 0.0
        return vx, vz
            
    
# LinearWave2D(rhoW, g, d, T, H, phi=0, x0=0, msg=True)
# rhoW, gravity, water-depth, time-period, wave-height, wave-phase, x0, printMsg




class LinearWaveDeep2D:
    # LinearWaveDeep2D(rhoW, g, T, H, phi=0, x0=0, msg=True)
    # rhoW, gravity, time-period, wave-height, wave-phase, x0, printMsg
    def __init__(self, rhoW, g, T, H, phi=0.0, x0=0.0, msg=True):
        self.rhoW = rhoW
        self.g = g
        self.T = T
        self.H = H
        self.phi = phi
        self.x0 = x0
        self.L = g/2.0/np.pi * T**2
        self.k = 2*np.pi/self.L
        self.w = 2*np.pi/self.T
        if(msg):
            print('Wave-Length L = ', self.L)
        
    def wavePhase(self, t, x):
        return self.k*(x-self.x0) - self.w*t + self.phi
        
    
    def waveElevation(self, t, x):
        et = self.H / 2 * np.cos( self.wavePhase(t, x) )
        return et
    
    
    def pressureTotPoi(self, t, x, z):
        et = self.waveElevation(t, x)
        if(z>0):
            if(z>et):
                pdyn = 0.0*et
            else:
                pdyn = et-z
        else:
            if(z<et):
                Kp = np.exp( self.k * z )
                pdyn = Kp*et - z
            else:
                pdyn=0.0*et
        return self.rhoW * self.g * pdyn
    
    def pressureTot(self, t, x, z):
        PdynRes = [self.pressureTotPoi(t,ix,iz) for ix, iz in zip(x,z)]        
        return PdynRes
    
    
    def pressureDynPoi(self, t, x, z):
        et = self.waveElevation(t, x)
        if(z>0):
            if(z>et):
                pdyn = 0.0
            else:
                pdyn = et-z
        else:
            if(z<et):
                Kp = np.exp( self.k * z )
                pdyn = Kp*et
            else:
                pdyn=0.0*et
        return self.rhoW * self.g * pdyn
    
    def pressureDyn(self, t, x, z):
        PdynRes = [self.pressureDynPoi(t,ix,iz) for ix, iz in zip(x,z)]        
        return PdynRes
    
    
    def particleVelPoi(self, t, x, z):
        et = self.waveElevation(t, x)
        if(z>0):
            vx = 0.0 * et
            vz = 0.0 * et
        else:
            mag = self.H/2 * self.w
            vx =  mag * np.exp( self.k * z )
            vx = vx * np.cos( self.wavePhase(t, x) )
            vz =  mag * np.exp( self.k * z )
            vz = vz * np.sin( self.wavePhase(t, x) )
        return vx, vz
    
    
    def particleAccPoi(self, t, x, z):
        et = self.waveElevation(t, x)
        if(z>0):
            vx = 0.0 * et
            vz = 0.0 * et
        else:
            mag = self.H/2 * self.w**2
            vx =  mag * np.exp( self.k * z )
            vx = vx * np.sin( self.wavePhase(t, x) )
            vz =  mag * np.exp( self.k * z )
            vz = -vz * np.cos( self.wavePhase(t, x) )
        return vx, vz

    
    def particleVelMax(self, x, z):
        mag = self.H/2 * self.w
        vx =  mag * np.exp( self.k * z )
        vz =  mag * np.exp( self.k * z )
        if(z>0):
            vx = vx * 0.0
            vz = vz * 0.0
        return vx, vz


    def particleAccMax(self, x, z):        
        mag = self.H/2 * self.w**2
        vx =  mag * np.exp( self.k * z )
        vz =  mag * np.exp( self.k * z )
        if(z>0):
            vx = vx * 0.0
            vz = vz * 0.0
        return vx, vz  
    
# LinearWaveDeep2D(rhoW, g, T, H, phi=0, x0=0, msg=True)
# rhoW, gravity, time-period, wave-height, wave-phase, x0, printMsg





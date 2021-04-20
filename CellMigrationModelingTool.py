import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import levy

class point:
    def __init__(self, xyz = (None, None, None), spherical = (None, None, None), t = None):
        self.__xyz = xyz
        self.__spherical = spherical
    
    @property
    def x(self):
        return self.__xyz[0]
    
    @property
    def y(self):
        return self.__xyz[1]
    
    @property
    def z(self):
        return self.__xyz[2]
    
    @property
    def xyz(self):
        return self.__xyz
    
    @property
    def mag(self):
        return self.__spherical[0]
    
    @property
    def polar_angle(self):
        return self.__spherical[1]
    
    @property
    def azimuthal_angle(self):
        return self.__spherical[2]
    
    @property
    def spherical(self):
        return self.__spherical
    
    @x.setter
    def x(self, x):
        xyz = np.array(self.__xyz)
        xyz[xyz==None] = 0
        xyz[0] = x
        self.__xyz = tuple(xyz)
        self.__shperical = spherical(self.__xyz)
        
    @y.setter
    def y(self, y):
        xyz = np.array(self.__xyz)
        xyz[xyz==None] = 0
        xyz[1] = y
        self.__xyz = tuple(xyz)
        self.__shperical = spherical(self.__xyz)
        
    @z.setter
    def x(self, z):
        xyz = np.array(self.__xyz)
        xyz[xyz==None] = 0
        xyz[2] = z
        self.__xyz = tuple(xyz)
        self.__shperical = spherical(self.__xyz)
    
    @mag.setter
    def mag(self, mag):
        
        
    @spherical.setter
    def spherical(self, xyz):
        x, y, z = xyz[0], xyz[1], xyz[2]
        magnitude = np.sqrt(x**2 + y**2 + z**2)
        if magnitude == 0:
            polar_angle = None
            azimuthal_angle = None
        elif x == 0:
            azimuthal_angle = None
            polar_angle = np.arccos(z/magnitude)
        else:
            polar_angle = np.arccos(z/magnitude)
            azimuthal_angle = np.arctan(y/x)
        self.__spherical = (magnitude, polar_angle, azimuthal_angle)
    
    @xyz.setter
    def xyz(self, spherical):
        magnitude = spherical[0]
        polar_angle = spherical[1]
        azimuthal_angle = spherical[2]
        x = magnitude * np.sin(polar_angle) * np.cos(azimuthal_angle)
        y = magnitude * np.sin(polar_angle) * np.sin(azimuthal_angle)
        z = magnitude * np.cos(polar_angle)
        self.__xyz = (x, y, z)
        return self.__xyz
            

            
class vector:
    def __init__(self, mag, polarangle, azimuthalangle, start, end):

        
class trajectory:
    def __init__(points = [], vectors = [], start = point(x = 0, y = 0), timestep = 1):
        pass
    
    

class LevyFlight:
    def __init__(self, loc = 0, scale = 1, n_tracks = 1, len_tracks = 100, 
                 timestep = 1):
        self.loc = loc
        self.scale = scale
        self.n_tracks = n_tracks
        self.len_tracks = len_tracks
        self.tracks = []
    def run(self):
        start_pos = 
        direction = np.random.random_sample(size = self.len_tracks) * 2* np.pi
        stepsize = levy.rvs(self.loc,self.scale,size = self.len_tracks)
        
LF = LevyFlight(10)
print(LF.repeats)
print(LF.tracks)
x = np.linspace(0, 20, 100) 
y1 = levy .pdf(x, 1, 3) 
y2 = levy .pdf(x, 0, 4) 
plt.plot(x, y1, "*", x, y2, "r--") 
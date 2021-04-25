import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import datetime
from scipy.stats import levy

class point:
    def __init__(self, xyz = (None, None, None), rpa = (None, None, None),
                 t = datetime.datetime.now(), n = None):
        self._xyz = xyz
        self._rpa = rpa
        self._t = t
        self._n = n
    
    @property
    def x(self):
        return self._xyz[0]
    
    @property
    def y(self):
        return self._xyz[1]
    
    @property
    def z(self):
        return self._xyz[2]
    
    @property
    def xyz(self):
        return self._xyz
    
    @property
    def radius(self):
        return self._rpa[0]
    
    @property
    def polar_angle(self):
        return self._rpa[1]
    
    @property
    def azimuthal_angle(self):
        return self._rpa[2]
    
    @property
    def rpa(self):
        return self._rpa
    
    @property
    def t(self):
        return self._t
    
    @property 
    def n(self):
        return self._n
    
    @x.setter
    def x(self, x):
        xyz = np.array(self._xyz)
        xyz[xyz==None] = 0
        xyz[0] = x
        self._xyz = tuple(xyz)
        self._rpa = self.cartesian2spherical(self._xyz)
        
    @y.setter
    def y(self, y):
        xyz = np.array(self._xyz)
        xyz[xyz==None] = 0
        xyz[1] = y
        self._xyz = tuple(xyz)
        self._rpa = self.cartesian2spherical(self._xyz)
        
    @z.setter
    def z(self, z):
        xyz = np.array(self._xyz)
        xyz[xyz==None] = 0
        xyz[2] = z
        self._xyz = tuple(xyz)
        self._rpa = self.cartesian2spherical(self._xyz)
    
    @xyz.setter
    def xyz(self, xyz):
        self._xyz = xyz
        self._rpa = self.cartesian2spherical(self._xyz)
    
    @radius.setter
    def radius(self, radius):
        rpa = np.array(self._rpa)
        rpa[rpa==None] = 0
        rpa[0] = radius
        self._rpa = tuple(rpa)
        self._xyz = self.spherical2cartesian(self._rpa)
        
    @polar_angle.setter
    def polar_angle(self, polar_angle):
        rpa = np.array(self._rpa)
        rpa[rpa==None] = 0
        rpa[1] = polar_angle
        self._rpa = tuple(rpa)
        self._xyz = self.spherical2cartesian(self._rpa)
        
    @azimuthal_angle.setter
    def azimuthal_angle(self, azimuthal_angle):
        rpa = np.array(self._rpa)
        rpa[rpa==None] = 0
        rpa[2] = azimuthal_angle
        self._rpa = tuple(rpa)
        self._xyz = self.spherical2cartesian(self._rpa)

    @rpa.setter
    def rpa(self, rpa):
        self._rpa = rpa
        self._xyz = self.spherical2cartesian(self._rpa)
    
    @t.setter
    def t(self, timestamp):
        self._t = timestamp
    
    @n.setter
    def n(self, point_number):
        self._n = point_number
        
    def cartesian2spherical(self, xyz):
        x, y, z = xyz[0], xyz[1], xyz[2]
        radius = np.sqrt(x**2 + y**2 + z**2)
        if radius == 0:
            return (0,0,0)
        polar_angle = np.arccos(z/radius)
        if polar_angle == 0:
            return (radius,0,0)
        if y >=0:
            azimuthal_angle = np.arctan2(y,x)
        else:
            azimuthal_angle = np.arctan2(y,x) + 2 * np.pi
        #azimuthal_angle = np.arcsin(y/(radius*np.sin(polar_angle)))
        #azimuthal_angle = np.arctan(y/x)
        return (radius, polar_angle, azimuthal_angle)
    
    def spherical2cartesian(self, rpa):
        radius = rpa[0]
        polar_angle = rpa[1]
        azimuthal_angle = rpa[2]
        x = radius * np.sin(polar_angle) * np.cos(azimuthal_angle)
        y = radius * np.sin(polar_angle) * np.sin(azimuthal_angle)
        z = radius * np.cos(polar_angle)
        return (x, y ,z)

            
class vector:
    def __init__(self, mpa = (None, None, None), start = point(), end = point()):
        self._mpa = mpa
        self._start = start
        self._end = end
        self._dt = self._end.t - self._start.t
        if all(self._start.xyz) and all(self._end.xyz):
            self._mpa = self.start_end2mpa(self._start, self._end)
        elif(all(self._start.xyz) and all (self._mpa)):
            self._end.xyz = self.start_mpa2end(self._start, self._mpa)
        elif(all(self._end.xyz) and all(self._mpa)):
            self.start.xyz = self.end_mpa2start(self._end, self._mpa)

    @property
    def mpa(self):
        return self._mpa
    @property
    def magnitude(self):
        return self._mpa[0]
    @property
    def polar_angle(self):
        return self._mpa[1]
    @property
    def azimuthal_angle(self):
        return self._mpa[2]
    @property
    def start(self):
        return self._start
    @property
    def end(self):
        return self._end
    @property
    def dx(self):
        return self._end.x - self._start.x
    @property
    def dy(self):
        return self._end.y - self._start.y
    @property
    def dz(self):
        return self._end.z - self._start.z
    @property
    def dt(self):
        return self._dt
    @property
    def velocity(self):
        return self._mpa[0]/self._dt.seconds
    
    @start.setter
    def start(self, start):
        self._start = start
        if all(self._end.xyz):
            self._mpa = self.start_end2mpa(self._start, self._end)
        else:
            self._end.xyz = (0,0,0)
            self._mpa = self.start_end2mpa(self._start, self._end)

    @end.setter
    def end(self, end):
        self._end = end
        if all(self._start.xyz):
            self._mpa = self.start_end2mpa(self._start, self._end)
        else:
            self._start.xyz = (0, 0, 0)
            self._mpa = self.start_end2mpa(self._start, self._end)
            
    @mpa.setter
    def mpa(self, mpa):
        self._mpa = mpa
        if all(self._start.xyz):
            self._end.xyz = self.start_mpa2end(self._start, self._mpa)
        elif all(self._end.xyz):
            self._start = self.end_mpa2start(self._end, self._mpa)
        else:
            self._start.xyz = (0, 0, 0)
            self._end.xyz = self.start_mpa2end(self._start, self._mpa)            
                
    def start_mpa2end(self, start, mpa):
        return tuple(np.array(start.xyz) + 
                              np.array(point().spherical2cartesian(mpa)))
    
    def end_mpa2start(self, end, mpa):
        return tuple(np.array(end.xyz) - 
                              np.array(point().spherical2cartesian(mpa)))
    
    def start_end2mpa(self, start, end):
        dx = end.x - start.x
        dy = end.y - start.y
        dz = end.z - start.z
        return point().cartesian2spherical((dx, dy, dz))


class trajectory:
    def __init__(points = [], vectors = [], start = "hallo", timestep = 1):
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
        start_pos = 1
        direction = np.random.random_sample(size = self.len_tracks) * 2* np.pi
        stepsize = levy.rvs(self.loc,self.scale,size = self.len_tracks)
        


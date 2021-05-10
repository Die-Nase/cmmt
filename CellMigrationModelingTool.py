import numpy as np
import pandas as pd
import datetime
import pint

from scipy.stats import levy
import matplotlib.pyplot as plt
import pint_pandas

ureg = pint.UnitRegistry()
Q_ = ureg.Quantity

                 
class point:
    
    def __init__(self,  xyz = (None, None, None), rpa = (None, None, None), 
                 t = datetime.datetime.now(), n = None, 
                 spaceUnit = 'meter', angle_measure = 'radian'):
        self._xyz= Q_(xyz, spaceUnit)
        self._t = t
        self._n = n
        if self.chk4input(rpa) and self.chk4input(xyz):
            raise Exception("Two inputs where one was expected. rpa will be ignored")
        elif self.chk4input(rpa):
            self._xyz = self.spherical2cartesian(rpa, spaceUnit, angle_measure)
            
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
        return self.rpa[0]
    
    @property
    def polar_angle(self):
        return self.rpa[1][0]
    
    @property
    def azimuthal_angle(self):
        return self.rpa[1][1]
    
    @property
    def rpa(self):
        return self.cartesian2spherical(self._xyz.magnitude, self._xyz.units)
    
    @property
    def t(self):
        return self._t
    
    @property 
    def n(self):
        return self._n
    
    @property 
    def spaceUnit(self):
        return self._xyz.units
    
    # @x.setter
    # def x(self, x):
        
    #     if len(x) == 2:
    #         x, unit = x
    #         Q_(x, unit).
    #     xyz = self._xyz.magnitude
    #     unit = self._xyz.units
    #     xyz[0] = x
    #     self._xyz = Q_(xyz, unit)
        
    # @y.setter
    # def y(self, y):
    #     xyz = self._xyz.magnitude
    #     unit = self._xyz.units
    #     xyz[1] = y
    #     self._xyz = Q_(xyz, unit)
        
    # @z.setter
    # def z(self, z):
    #     xyz = self._xyz.magnitude
    #     unit = self._xyz.units
    #     xyz[2] = z
    #     self._xyz = Q_(xyz, unit)
    
    # @xyz.setter
    # def xyz(self, xyz):
    #     unit = self._xyz.units
    #     self._xyz = Q_(xyz, unit)
    
    # @radius.setter
    # def radius(self, radius):
    #     rpa = list(self._rpa)
    #     rpa[0] = radius
    #     self._rpa = tuple(self.quality_control(rpa, [0]*3, 
    #                     [self._spaceUnit,ureg('radians'),ureg('radians')]))
    #     self._xyz = self.spherical2cartesian(self._rpa)
        
    # @polar_angle.setter
    # def polar_angle(self, polar_angle):
    #     rpa = list(self._rpa)
    #     rpa[1] = polar_angle
    #     self._rpa = tuple(self.quality_control(rpa, [0]*3, 
    #                     [self._spaceUnit,ureg('radians'),ureg('radians')]))
    #     self._xyz = self.spherical2cartesian(self._rpa)
        
    # @azimuthal_angle.setter
    # def azimuthal_angle(self, azimuthal_angle):
    #     rpa = list(self._rpa)
    #     rpa[2] = azimuthal_angle
    #     self._rpa = tuple(self.quality_control(rpa, [0]*3, 
    #                     [self._spaceUnit,ureg('radians'),ureg('radians')]))
    #     self._xyz = self.spherical2cartesian(self._rpa)

    # @rpa.setter
    # def rpa(self, rpa):
    #     self._rpa = tuple(self.quality_control(list(rpa), [0]*3, 
    #             [self._spaceUnit,ureg('radians'),ureg('radians')]))
    #     self._xyz = self.spherical2cartesian(self._rpa)
    
    # @t.setter
    # def t(self, timestamp):
    #     self._t = timestamp
    
    # @n.setter
    # def n(self, point_number):
    #     self._n = point_number
        
    # @spaceUnit.setter
    # def spaceUnit(self, spaceUnit):
    #     self._xyz._units = ureg.Unit(spaceUnit)
        
    def cartesian2spherical(self, xyz, spaceUnit):
        x, y, z = Q_(xyz, spaceUnit)
        radius = np.sqrt(x**2 + y**2 + z**2)
        if radius == 0:
            return (radius, Q_([0,0],'radian'))
        polar_angle = np.arccos(z/radius)
        if polar_angle == 0:
            return (radius, Q_([0,0],'radian'))
        if y >=0:
            azimuthal_angle = np.arctan2(y,x)
        else:
            azimuthal_angle = np.arctan2(y,x) + 2 * np.pi * ureg('radians')
        #azimuthal_angle = np.arcsin(y/(radius*np.sin(polar_angle)))
        #azimuthal_angle = np.arctan(y/x)
        return (radius, Q_([polar_angle.magnitude, azimuthal_angle.magnitude], 'radian'))
    
    def spherical2cartesian(self, rpa, spaceUnit, angle_measure):
        radius = Q_(rpa[0], spaceUnit)
        polar_angle, azimuthal_angle = Q_([rpa[1],rpa[2]], angle_measure)
        x = radius * np.sin(polar_angle) * np.cos(azimuthal_angle)
        y = radius * np.sin(polar_angle) * np.sin(azimuthal_angle)
        z = radius * np.cos(polar_angle)
        return Q_([x.magnitude, y.magnitude ,z.magnitude],spaceUnit)

    # def quality_control(self, input_values, default_values, default_Units):
    #     for i in range(len(input_values)):
    #         if input_values[i] == None:
    #             input_values[i] = default_values[i]* default_Units[i]
    #         elif Q_(input_values[i]).dimensionless:
    #             input_values[i] = input_values[i] * default_Units[i]
    #     return input_values
    
    def chk4input(self, tuple_in):
        return not any(map(lambda x: x is None, tuple_in))


class vector:
    def __init__(self, start = point(), end = point()):
        self._start = start
        self._end = end
    
    @property
    def mpa(self):
        return self.start_end2mpa(self._start,self._end)
    @property
    def magnitude(self):
        mpa = self.start_end2mpa(self._start,self._end)
        return mpa[0]
    @property
    def polar_angle(self):
        mpa = self.start_end2mpa(self._start,self._end)
        return mpa[1][0]
    @property
    def azimuthal_angle(self):
        mpa = self.start_end2mpa(self._start,self._end)
        return mpa[1][1]
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
        return (self._end.t - self._start.t).total_seconds() * ureg('seconds')
    @property
    def velocity(self):
        return self.magnitude/self.dt
    
    # @start.setter
    # def start(self, start):
    #     self._start = start
    #     if all(self._end.xyz):
    #         self._mpa = self.start_end2mpa(self._start, self._end)
    #     else:
    #         self._end.xyz = tuple((0,0,0) * self._end.spaceUnit)
    #         self._mpa = self.start_end2mpa(self._start, self._end)
    
    # @end.setter
    # def end(self, end):
    #     self._end = end
    #     if all(self._start.xyz):
    #         self._mpa = self.start_end2mpa(self._start, self._end)
    #     else:
    #         self._start.xyz = tuple((0,0,0) * self._start.spaceUnit)
    #         self._mpa = self.start_end2mpa(self._start, self._end)
        
    # @mpa.setter
    # def mpa(self, mpa):
    #     self._mpa = mpa
    #     if all(self._start.xyz):
    #         self._end.xyz = self.start_mpa2end(self._start, self._mpa)
    #     elif all(self._end.xyz):
    #         self._start = self.end_mpa2start(self._end, self._mpa)
    #     else:
    #         self._start.xyz = tuple((0,0,0) * self._start.spaceUnit)
    #         self._end.xyz = self.start_mpa2end(self._start, self._mpa)            
                
    # def start_mpa2end(self, start, mpa):
    #     return tuple(ureg.Quantity.from_list(list(start.xyz)) +
    #                  ureg.Quantity.from_list(list(point().spherical2cartesian(mpa))))
    
    # def end_mpa2start(self, end, mpa):
    #      return tuple(ureg.Quantity.from_list(list(end.xyz)) +
    #                  ureg.Quantity.from_list(list(point().spherical2cartesian(mpa))))
    
    def start_end2mpa(self, start, end):
        xyz = end.xyz - start.xyz
        return point().cartesian2spherical(xyz.magnitude, xyz.units)


class trajectory:
    def __init__(self, fxyz = np.empty((1,5)), generation_datetime = datetime.datetime.now(),
                 freq = '1 S', spaceUnit = 'm'):
        t = pd.Series(pd.date_range(generation_datetime, freq=freq, periods=len(fxyz)))
        data ={'fn': pd.Series(fxyz[:,0], dtype='Int64'),
               't': pd.Series(t, dtype = 'datetime64'),
               # 'x': pd.Series(fxyz[:,1], dtype="pint["+spaceUnit+"]"),
               # 'y': pd.Series(fxyz[:,2], dtype="pint["+spaceUnit+"]"),
               # 'z': pd.Series(fxyz[:,3], dtype="pint["+spaceUnit+"]" )}
               'x': pd.Series(fxyz[:,1]),
               'y': pd.Series(fxyz[:,2]),
               'z': pd.Series(fxyz[:,3])}
        self._txyz = pd.DataFrame(data)
        self._points = []
        self._vectors = []
    
    
    @property
    def txyz(self):
        return self._txyz
    
    @txyz.setter
    def txyz(self, txyz):
        self._txyz = txyz
        
    @property
    def points(self):
        if not self._points:
            self._points = self.get_points(self._txyz)
        return self._points
    
    def get_points(self, txyz):
        points = [None] * len(txyz)
        for row in txyz.itertuples():
            points[row.Index] = point(xyz = (row.x, row.y, row.z),
                                          t = row.t, n = row.fn)
        return points
    
    @property
    def vectors(self):
        if not self._vectors:
            self._vectors = self.get_vectors(self.points)
        return self._vectors
    
    def get_vectors(self, points):
        vectors = [None] * (len(points)-1)
        for i in range(0, len(vectors)):
            vectors[i] = vector(start = points[i], end = points[i+1])
        return vectors
            
    
    @property
    def velocities(self):
        velocity_list = [None] * len(self._vectors)
        for i in range(0, len(velocity_list)):
            velocity_list[i] = self._vectors[i].velocity
        return velocity_list
            
    
    @property
    def displacements(self):
        displacement_list = [None] * len(self._vectors)
        for i in range(0, len(displacement_list)):
            displacement_list[i] = self._vectors[i].magnitude
        return displacement_list
    
            
        
    def trackmate_xml2DataFrame(self, trackmate_xml):
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
        

fxyz = np.array([[181, 372.44840892346417, 1392.9970222092763, 0.0],
               [182, 373.16006520865074, 1393.1116296971932, 0.0],
               [183, 372.5888414573269, 1392.087922094057, 0.0],
               [184, 371.47162594500435, 1391.5351846322474, 0.0],
               [185, 369.24421493144075, 1390.6588627745375, 0.0],
               [186, 369.6450338418909, 1390.725567955707, 0.0],
               [187, 370.624446126946, 1391.4296428129408, 0.0],
               [188, 371.50277040471093, 1391.6409472786852, 0.0],
               [189, 372.94501560620296, 1392.0921455475623, 0.0],
               [190, 373.18290170287753, 1391.7978370942435, 0.0],
               [191, 371.2552424409938, 1391.4683611969374, 0.0],
               [192, 369.7177104170941, 1387.0306214492473, 0.0],
               [193, 368.7203413928448, 1390.0195814760107, 0.0],
               [194, 370.00135913784425, 1386.7504562823817, 0.0],
               [195, 370.0178789337686, 1386.7260815938664, 0.0]])

track = trajectory(fxyz = fxyz, freq = '30 S', spaceUnit = "um")
track.points
print(track.vectors)

import pandas as pd
import pint


df = pd.DataFrame({
    "torque": pd.Series([1, 2, 2, 3], dtype="pint[lbf ft]"),
    "angular_velocity": pd.Series([1, 2, 2, 3], dtype="pint[rpm]") })

df['power'] = df['torque'] * df['angular_velocity']
df.dtypes
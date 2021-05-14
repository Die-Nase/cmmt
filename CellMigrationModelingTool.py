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
    
    def __init__(self,  xyz = [None, None, None], rpa = (None, None, None), 
                 t = datetime.datetime.now(), n = None, 
                 spaceUnit = 'meter', angle_measure = 'radian'):
        self._xyz= Q_(xyz, spaceUnit)
        self._t = t
        self._n = n
        if self.chk4input(rpa) and self.chk4input(xyz):
            raise Exception("Two inputs where one was expected. rpa will be ignored")
        elif self.chk4input(rpa):
            self._xyz = self.spherical2cartesian(rpa, spaceUnit, angle_measure)
            
    def __eq__(self, other):
        if type(other) == point:    
            if all(self.xyz == other.xyz) and self.t == other.t and self.n == other.n:
                return True
            else:
                return False
        else: 
            raise Exception("cannot compare {other_type} to {self_type}".format(
                other_type = type(other), self_type = type(self)))
    
    def __add__(self, other):
        if  isinstance(other, point):
            return vector(start = self, end = other)
        elif isinstance(other, vector):
            xyz = self.xyz + (other.end.xyz - other.start.xyz)
            return point(xyz = xyz.magnitude, spaceUnit = xyz.units)
        # add dxyz and dmpa
        else:
            raise Exception("cannot add {other_type} to {self_type}".format(
                other_type = type(other), self_type = type(self)))
    
    def __sub__(self, other):
        if  type(other) == point:
            return vector(start = other, end = self)
        elif type(other) == vector:
            xyz = self.xyz - (other.end.xyz - other.start.xyz)
            return point(xyz = xyz.magnitude, spaceUnit = xyz.units)
        #  add  dxyz and dmpa
        else:
            raise Exception("cannot subtract {other_type} from {self_type}".format(
                other_type = type(other), self_type = type(self)))
    
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
 
    def chk4input(self, tuple_in):
        return not any(map(lambda x: x is None, tuple_in))


class vector:
    def __init__(self, start = point(), end = point()):
        self._start = start
        self._end = end
        
    def __eq__(self, other):
        if  type(other) == vector:
            if self.start == other.start and self.end == other.end:
                return True
            else:
                return False
        else:
            raise Exception("cannot compare {other_type} to {self_type}".format(
                other_type = type(other), self_type = type(self)))
    
    def __add__(self, other):
        if type(other) == vector:
            end = self.end + (other.end - other.start)
            return vector(start = self.start, end = end)
        else:
            raise Exception("cannot add {other_type} to {self_type}".format(
                other_type = type(other), self_type = type(self)))
        
    def __sub__(self, other):
        if type(other) == vector:
            end = self.end - (other.end - other.start)
            return vector(start = self.start, end = end)
        else:
            raise Exception("cannot substract {other_type} from {self_type}".format(
                other_type = type(other), self_type = type(self)))
            
    def dot_product(self, other):
        return (self.dx * other.dx + self.dy * other.dy + self.dz * other.dz)
    
    def cross_product(self, other):
        dx = self.dy * other.dz - self.dz * other.dy
        dy = self.dz * other.dx - self.dx * other.dz
        dz = self.dx * other.dy - self.dy * other.dx
        dxdydz = Q_.from_list([dx, dy, dz])
        return self.dxdydz2vector(dxdydz)
    
    def normalize(self):
        dxdydz = self.dxdydz/self.magnitude
        return self.dxdydz2vector(dxdydz, start = self.start)
    
    def translate(self, start_point):
        return self.dxdydz2vector(self.dxdydz, start_point = start_point)
    
    def translate2origin(self):
        return self.translate(point(xyz = [0,0,0]))
    
    def dxdydz2vector(self,dxdydz, start_point = point(xyz = [0,0,0])):
        return vector(start = start_point, 
                      end = point(xyz = dxdydz.magnitude,
                      spaceUnit = dxdydz.units))
    
    def mpa2vector(self, mpa, start_point = point(xyz = [0,0,0])):
        rpa = (mpa[0].magnitude, mpa[1][0].magnitude, mpa[1][1].magnitude)
        end_point = point(rpa = rpa, spaceUnit = mpa[0].units)
        return vector(start = start_point, end = end_point)

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
    def dxdydz(self):
        return self.end.xyz - self.start.xyz
    @property
    def dt(self):
        return (self._end.t - self._start.t).total_seconds() * ureg('seconds')
    @property
    def velocity(self):
        return self.magnitude/self.dt
    
    def start_end2mpa(self, start, end):
        xyz = end.xyz - start.xyz
        return point().cartesian2spherical(xyz.magnitude, xyz.units)


class trajectory:
    def __init__(self, fxyz = np.empty((1,5)), generation_datetime = datetime.datetime.now(),
                 freq = '1 S', spaceUnit = 'm'):
        t = pd.Series(pd.date_range(generation_datetime, freq=freq, periods=len(fxyz)))
        data ={'fn': pd.Series(fxyz[:,0], dtype='Int64'),
               't': pd.Series(t, dtype = 'datetime64'),
                'x': pd.Series(fxyz[:,1], dtype="pint["+spaceUnit+"]"),
                'y': pd.Series(fxyz[:,2], dtype="pint["+spaceUnit+"]"),
                'z': pd.Series(fxyz[:,3], dtype="pint["+spaceUnit+"]" )}

        self._txyz = pd.DataFrame(data)
        self._points = []
        self._vectors = []
    
    
    @property
    def txyz(self):
        return self._txyz
    
    # @txyz.setter
    # def txyz(self, txyz):
    #     self._txyz = txyz
        
    @property
    def points(self):
        if not self._points:
            self._points = self.get_points(self._txyz)
        return self._points

    @property
    def vectors(self):
        if not self._vectors:
            self._vectors = self.get_vectors(self.points)
        return self._vectors
            
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

    def get_points(self, txyz):
        points = [None] * len(txyz)
        for row in txyz.itertuples():
            xyz_list = [row.x, row.y, row.z]
            xyz = Q_.from_list(xyz_list)
            points[row.Index] = point(xyz = xyz.magnitude, spaceUnit = xyz.units,
                                      t = row.t, n = row.fn)
        return points
        
    def get_vectors(self, points):
        vectors = [None] * (len(points)-1)
        for i in range(0, len(vectors)):
            vectors[i] = vector(start = points[i], end = points[i+1])
        return vectors
        
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

track = trajectory(fxyz = fxyz, freq = '30 S', spaceUnit = "um", generation_datetime = datetime.datetime.now())
track.points

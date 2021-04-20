import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import levy

class point:
    def __init__(self, x = 0, y = 0, z = 0):
        self.x = x
        self.y = y
        self.z = z
    
    def polar(self):
        pass

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
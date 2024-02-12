#!/usr/bin/env python
#
#         python find_xyz.py 2010-11-15T17:26:45.808
#
#           reads the "orbital_elements" file
#
# uses Pyephem
#    see http://www.rhodesmill.org/brandon/projects/pyephem.html
#

from ephem import *
import time, math, subprocess, sys, os
import datetime
import matplotlib.pyplot as plt
import numpy as np

"""def increment_date(year, month, day, time):
    if month in ['09', '04', '06', '11']:
        max_days = 30
    if month """

# date format required    2019-02-19T18:00:00

start_date = sys.argv[1]
no_of_days = int(sys.argv[2])

year = start_date[0:4]
month = start_date[5:7]
day = start_date[8:10]
time = start_date[11:]
start_date = datetime.datetime(int(year),int(month), int(day), int(time[0:2]), int(time[3:5])*60 + int(time[6:]))

dates = []
suns = []
minors = []
diffs = []

elements = open('orbital_elements')
minor = EllipticalBody()
minor._epoch = '2000/01/01.5'   # epoch date for J2000.0
while 1:
    line = elements.readline()
    if line == "": break
    word = line.split()
    if line.find('Object:') >= 0: obj_name = ' '.join(word[2:])
    if line.find('Mean anomaly') >= 0: minor._M = float(word[3])
    if line.find('Argument of perihelion') >= 0: minor._om = float(word[4])
    if line.find('Long. of ascending node') >= 0: minor._Om = float(word[5])
    if line.find('Inclination') >= 0: minor._inc = float(word[2])
    if line.find('Eccentricity') >= 0: minor._e = float(word[2])
    if line.find('Semimajor axis') >= 0: minor._a = float(word[3])
    if line.find('Epoch of osculation') >= 0:
        minor._epoch_M = word[4]+'/'+word[5]+'/'+word[6]
    if line.find('Absolute magnitude') >= 0: minor._H = float(word[4])
    if line.find('Slope parameter') >= 0: minor._G = float(word[4])   

durham = Observer()
durham.epoch = '2000'
durham.long = '-1:34:24'
durham.lat = '54:46:01'
durham.temp = 10
durham.elev = 119.5

for i in range(0, no_of_days):
    durham.date =start_date + datetime.timedelta(days=i) # year+'/'+month+'/'+day+ ' '+obstime

    dates.append(durham.date)

    const = 180./math.pi
    minor.compute(durham)
    pos = Equatorial(minor.a_ra,minor.a_dec, epoch='2000')
    ecl = Ecliptic(pos)
    ecl_long = float(ecl.lon)
    ecl_lat  = float(ecl.lat)

    d_minor = minor.earth_distance
    x_minor = d_minor*math.cos(ecl_lat)*math.cos(ecl_long)
    y_minor = d_minor*math.cos(ecl_lat)*math.sin(ecl_long)
    z_minor = d_minor*math.sin(ecl_lat)
    print ("minor", x_minor, y_minor, z_minor)
    print (" ")
    minors.append([x_minor, y_minor, z_minor])

    xxx = Sun()
    xxx.compute(durham)
    pos = Equatorial(xxx.a_ra,xxx.a_dec, epoch='2000')
    ecl = Ecliptic(pos)
    ecl_long = float(ecl.lon)
    ecl_lat  = float(ecl.lat)
    #print ecl_long*180./math.pi , ecl_lat*180./math.pi 
    d_sun = xxx.earth_distance
    x_sun = d_sun*math.cos(ecl_lat)*math.cos(ecl_long)
    y_sun = d_sun*math.cos(ecl_lat)*math.sin(ecl_long)
    z_sun = d_sun*math.sin(ecl_lat)
    print ("Sun  ",x_sun, y_sun, z_sun)
    print (" ")
    print ("DIFF ",x_minor-x_sun, y_minor-y_sun, z_minor-z_sun)
    print(" ")
    suns.append([x_sun, y_sun, z_sun])
    diffs.append([x_minor-x_sun, y_minor-y_sun, z_minor-z_sun])

#return dates, minors, suns, diffs
    
earths = -1*np.array(suns)
    
fig = plt.figure()
ax = plt.axes(projection='3d')
ax.plot3D([i[0] for i in diffs], [i[1] for i in diffs], [i[2] for i in diffs], color = 'red', label = f'{obj_name}')
ax.plot3D([i[0] for i in earths], [i[1] for i in earths], [i[2] for i in earths], color = 'blue', label = 'Earth')
ax.plot3D([0], [0], [0], color = 'gold', label = 'Sun', marker = 'o', markersize = 12, linestyle = 'None')
plt.legend()
ax.set_xlabel('x (AU)'); ax.set_ylabel('y (AU)'); ax.set_zlabel('z (AU)')
ax.set_xlim(-2, 2)
ax.set_ylim(-2, 2)
plt.savefig(f'plots/{obj_name}_orbit')
plt.show()
    
print(dates)
    

#print (sys.argv[1], "%7.5f" % d_minor, "%7.5f" % (180./math.pi*ecl_long_minor), "%7.5f" % (180./math.pi*ecl_lat_minor), "%7.5f" % x_minor, "%7.5f" %  y_minor, "%7.5f" % z_minor,  "%7.5f" % d_sun, "%7.5f" % (180./math.pi*ecl_long_sun), "%7.5f" % (180./math.pi*ecl_lat_sun), "%7.5f" % x_sun, "%7.5f" %  y_sun, "%7.5f" % z_sun)

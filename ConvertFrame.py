from Tools import *
import math
from math import sin, cos

#Number of minuutes it takes for the earth to complete one sidereal rotation
SIDEREAL_PERIOD = 23.9344696*60

#The time of the 2000 spring equinox (minutes since epoch)
EQ_SINCE_EPOCH = 26937636.0

def convertFrame(frame1, frame2, u, t=0.0):
    """
    Convert from one reference frame to another, based on the time
    ECI: Earth Centered Inertial Frame(time independant)
    ECEF: Earth Centered Earth Inertial Frame (time dependant)
    coord: Polar coordinates: [Longitude, Latitude, Radius] (time dependant)
    """
    # theta is angle of earth from ECI. 0 is noon on vernal equinox
    theta = 2*math.pi*(t-EQ_SINCE_EPOCH)/SIDEREAL_PERIOD
    if frame1 == "ECI" and frame2 == "ECEF":
        R = [
            [cos(-theta), -sin(-theta), 0],
            [sin(-theta), cos(-theta), 0],
            [0, 0, 1]
        ]
        return mat_multi(R, u)
    elif frame1 == "ECEF" and frame2 == "coord":
        long = math.degrees(math.atan2(u.j, u.i))
        if long > 180:
            long -= 360
        lat  = math.degrees(math.atan2(u.k, (u.i**2+u.j**2)**.5))
        if lat > 90:
            lat -= 360
        h = mag(u)
        return Vec_3(long, lat, h)
    elif frame1 == "ECI" and frame2 == "coord":
        u2 = convertFrame("ECI", "ECEF", u, t)
        return convertFrame("ECEF", "coord", u2, t)
    print("Unknown frame conversion.")
    return u
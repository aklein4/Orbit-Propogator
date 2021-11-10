
from Tools import *
from Orbit import Orbit
from ConvertFrame import convertFrame, EQ_SINCE_EPOCH

import math
from matplotlib import pyplot as plt
import datetime
import time
import numpy as np

#gravity constant of the earth
MU_EARTH_S = 3.986004418 * (10**5) #(km^3)/(s^2)
MU_EARTH_MIN = 3.986004418 * (10**5) * 3600 #(km^3)/(min^2)

#Gravity constant of the sun
MU_SUN_S = 1.32712440018 * (10**11) #(km^3)/(s^2)
MU_SUN_MIN = 1.32712440018 * (10**11) * 3600 #(km^3)/(min^2)

#Radius of the Earth in km
RAD_EARTH = 6378

#time conversions
MIN_to_HOUR = 60
MIN_to_DAY = 60*24

#image for orbit plotting
IMG = plt.imread("grid_earth.jpg")

#Tracks the orbit of the sun
SUN = Orbit(MU_SUN_MIN)
SUN.set_params(.0167, 149600000, math.radians(23.5), math.radians(0), math.radians(288.1),
               math.radians(360-288.1), EQ_SINCE_EPOCH)

def main():

    #initialize ground map
    plt.figure(1)
    plt.xlim([-180, 180])
    plt.ylim([-90, 90])
    plt.imshow(IMG, extent=[-180, 180, -90, 90])
    plt.pause(.001)

    #initialize 3D graph
    plt.figure(2)
    space = plt.axes(projection="3d")
    space.set_xlim3d([-8000, 8000])
    space.set_ylim3d([-8000, 8000])
    space.set_zlim3d([-8000, 8000])
    u, v = np.mgrid[0:2 * np.pi:20j, 0:np.pi:10j]
    x = 6000 * np.cos(u) * np.sin(v)
    y = 6000 * np.sin(u) * np.sin(v)
    z = 6000 * np.cos(v)
    space.plot_surface(x, y, z, color="b")
    plt.pause(.001)

    #initialize sattelite dictionary
    sats = {
        "hubble": Orbit(MU_EARTH_MIN),
        "sun": SUN
    }
    #Hubble telescope automatically loaded
    sats["hubble"].set_params(.0002493, 6914, math.radians(28.4709), math.radians(323.96),
                  0, 0, 0, ma=False)


    """ Terminal Function: """
    while True:
        command = input("command<< ").lower().strip()
        if command == "exit":
            break
        command = command.split(",")
        for i in range(len(command)):
            command[i] = command[i].strip().lower()
        if command[0] == "range" and len(command) >= 3:
            lower = int("0"+command[1])
            upper = int("0"+command[2])
            step = 1
            if len(command) > 3:
                step = int("0"+command[3])
            command = input("ranged command<< ").lower().strip()
            command = command.split(",")
            for i in range(len(command)):
                command[i] = command[i].strip()
            for t in range(lower, upper, step):
                cop = command.copy()
                cop.append(str(t))
                print("--"+str(t)+"--")
                com(sats, cop, space)
        else:
            com(sats, command, space)
    return

def com(sats, command, space):
    if len(command) > 0:
        if command[0] == "new" or command[0] == "set" and len(command) == 1:
            name = input("name: ")
            sats[name] = Orbit(MU_EARTH_MIN)
            init = "x"
            while init != "parameters" and init != "state":
                init = input("parameters or state? ")
            if init == "parameters":
                e = float("0"+input("eccentricity: "))
                a = float("0"+input("semi-major axis: "))
                i = float("0"+input("inclination: "))
                omega = float("0"+input("angle of node: "))
                w = float("0"+input("argument of paragee: "))
                v0 = float("0"+input("anomoly: "))
                t0 = float("0"+input("t0: "))
                MA = input("Mean anomoly? ")
                if MA[0].lower() == "y":
                    MA = True
                else:
                    MA = False
                sats[name].set_params(e, a, i, omega, w, v0, t0, ma=MA)
            else:
                ri = float("0"+input("ri: "))
                rk = float("0"+input("rj: "))
                rj = float("0"+input("rk: "))
                veli = float("0"+input("veli: "))
                velj = float("0"+input("velj: "))
                velk = float("0"+input("velk: "))
                t0 = float("0"+input("t0: "))
                sats[name].set_state(Vec_3(ri, rj, rk), Vec_3(veli, velj, velk), t0)
            print("Satellite created: '"+name+"'")
        elif command[0] == "get" and len(command) >= 3:
            name = command[1]
            if not name in sats.keys():
                print("Invalid Satellite:", name)
            else:
                n = 0
                if len(command) > 3:
                    n = float(command[3])
                val = sats[name].get(command[2], t=n)
                if type(val).__name__ == "Vec_3":
                    val.show()
                else:
                    print(val)
        elif command[0] == "check_eclipse" and len(command) == 3:
            name = command[1]
            if not name in sats.keys():
                print("Invalid Satellite:", name)
            else:
                t = float("0"+command[2])
                check = check_eclipse(sats[name], SUN, t)
                if check:
                    print("In shade.")
                else:
                    print("In sun.")
        elif command[0] == "plot" and len(command) == 3:
            name = command[1]
            if not name in sats.keys():
                print(sats)
                print("Invalid Satellite:", name)
            else:
                sat = sats[name]
                t = float(command[2])
                r = sat.get("r", t)
                #plot ground map
                plt.figure(1)
                coord = convertFrame("ECI", "coord", r, t)
                long = coord.i
                lat = coord.j
                shade = check_eclipse(sat, SUN, t)
                point = "r"
                if shade:
                    point += "x"
                else:
                    point += "."
                plt.plot(long, lat, point, markersize="3")
                plt.pause(.001)

        elif command[0] == "plot3D" and len(command) == 3:
            name = command[1]
            if not name in sats.keys():
                print(sats)
                print("Invalid Satellite:", name)
            else:
                sat = sats[name]
                t = float(command[2])
                r = sat.get("r", t)
                #plot space
                plt.figure(2)
                space.scatter(r.i, r.j, r.k, color="r")
                plt.pause(.001)

        elif command[0] == "clear" or command[0] == "open":
            #clear ground track
            plt.figure(1)
            plt.clf()
            plt.imshow(IMG, extent=[-180, 180, -90, 90])
            plt.pause(.001)
            #clear space
            plt.figure(2)
            space.clear()
            space.set_xlim3d([-8000, 8000])
            space.set_ylim3d([-8000, 8000])
            space.set_zlim3d([-8000, 8000])
            u, v = np.mgrid[0:2 * np.pi:20j, 0:np.pi:10j]
            x = 6000 * np.cos(u) * np.sin(v)
            y = 6000 * np.sin(u) * np.sin(v)
            z = 6000 * np.cos(v)
            space.plot_surface(x, y, z, color="b")
            plt.pause(.001)
            print("Map cleared.")

        elif command[0] == "eval":
            print(eval(input("expression: ")))

        elif command[0] == "date":
            year = int("0"+input("year: "))
            month = int("0" + input("month: "))
            day = int("0" + input("day: "))
            hour = int("0" + input("hour: "))
            min = int("0" + input("min: "))
            print("Minutes since epoch:")
            print(math.floor(datetime.datetime(year, month, day, hour, min).timestamp()/60))
        elif command[0] == "time":
            print("Current minutes since epoch:")
            print(math.floor(time.time()/60))

        else:
            print("Invalid command.")


def check_eclipse(sat, sun, t):
    """
    input the orbits of a satellite and the sun, along with a time (min since epoch)
    return whether or not the satellite is being eclipsed by the earth
    """
    r_sat = sat.get("r", t)
    r_sun = sun.get("r", t)
    if dot(r_sat, r_sun) >= 0:
        #sat on same side of earth as sun
        return False
    par = norm(r_sun)*dot(norm(r_sun), r_sat)
    perp = r_sat-par
    if mag(perp) < RAD_EARTH:
        #sat is behind earth's cylindical shadow
        return True
    return False

if __name__ == "__main__":
    main()
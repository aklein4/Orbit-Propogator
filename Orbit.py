from Tools import *
import math
from math import sin, cos

class Orbit:
    """
    Contains information about an orbit in constant reference frame.
    Used to track a satellite's movement
    """
    def __init__(self, mu):
        self._mu = mu
        self._e = 0 #eccentricity
        self._a = 0 #semimajor axis
        self._i = 0 #inclination
        self._omega = 0 #right ascension of the ascending node
        self._w = 0 #argument of paragee
        self._t0 = 0 #time when anomaly is 0
        self._initialized = False

    def set_params(self, e, a, i, omega, w, v0, t, ma=False):
        # v0 is the true anomoly at t
        self._e = e
        if(self._e) >= 1:
            print("Error: can not handle hyperbolic orbits.")
            return
        self._a = a
        self._i = i
        self._omega = omega
        self._w = w
        self._initialized = True
        #calculate t0
        M = v0
        if not ma:
            M = self.convertAnom("v", "m", v0)
        n = self.get("n")
        dt = M/n
        self._t0 = t - dt

    def set_state(self, r, vel, t):
        # Store the orbit parameters
        self._e = self._calc("e", r, vel)
        if(self._e) >= 1:
            print("Error: can not handle hyperbolic orbits.")
            return
        self._a = self._calc("a", r, vel)
        self._i = self._calc("i", r, vel)
        self._omega = self._calc("omega", r, vel)
        self._w = self._calc("w", r, vel)
        self._initialized = True
        #calculate t0
        v0 = self._calc("v", r, vel)
        M = self.convertAnom("v", "m", v0)
        n = self.get("n")
        dt = M/n
        self._t0 = t - dt

    def _calc(self, prop, r, vel):
        """
        Calculate orbit properties using state vectors.
        Requires no outside variables other than mu.
        :param prop: property to calculate
        :param r: position vector
        :param vel: velocity vector
        :return: property that was calculated
        """
        prop = prop.lower().strip()
        if prop == "e":
            # eccentricity
            return mag(self._calc("e_vec", r, vel))
        if prop == "a":
            # Semimajor axis
            En = (dot(vel, vel) / 2) - (self._mu / mag(r))
            a = -self._mu / (2 * En)
            return a
        if prop == "i":
            # Inclination
            h = cross(r, vel)
            return math.acos(h.k / mag(h))
        if prop == "omega":
            # right-ascension of ascending node
            node = self._calc("node", r, vel)
            omega = math.acos(node.i / mag(node))
            if node.j < 0:
                omega = 2 * math.pi - omega
            return omega
        if prop == "w":
            # argument of paragee
            node = self._calc("node", r, vel)
            e_vec = self._calc("e_vec", r, vel)
            w = math.acos(dot(node, e_vec) / (mag(node) * mag(e_vec)))
            if e_vec.k < 0:
                w = 2 * math.pi - w
            return w
        if prop == "v":
            # true anomoly
            e_vec = self._calc("e_vec", r, vel)
            v = math.acos(dot(e_vec, r) / (mag(e_vec) * mag(r)))
            if e_vec.k < 0:
                v = 2 * math.pi - v
            return v
        if prop == "node":
            # line of nodes
            h = cross(r, vel)
            node = Vec_3(-h.j, h.i, 0)
            if(mag(node)==0):
                node.i = 1
            return node
        if prop == "e_vec":
            # eccentricity vector
            h = cross(r, vel)
            c1 = cross(vel, h)/self._mu
            c2 = norm(r)
            e_vec = c1-c2
            return e_vec
        print("Error: invalid property calculation call.")
        return 0/0

    def get(self, prop, t=0.0):
        """
        Returns a property of the orbit.
        Requires all orbit parameters to be initialized, other than t0.
        :param prop: property to be returned
        :param t: current time to use in calculation
        :return: property
        """
        if not self._initialized:
            print("Orbit not initialized.")
            return 0
        prop = prop.lower().strip()
        if prop == "e": return self._e
        if prop == "a": return self._a
        if prop == "i": return self._i
        if prop == "omega": return self._omega
        if prop == "w": return self._w
        if prop == "t0": return self._t0
        if prop == "params":
            params = {
                "e": self._e,
                "a": self._a,
                "i": self._i,
                "omega": self._omega,
                "w": self._w,
                "v": self.get("va", t)
            }
            return params
        if prop == "n":
            # mean motion
            return self._mu ** (.5) / self._a ** (3 / 2)
        if prop == "period":
            # period time
            return 2*math.pi/self.get("n", t)
        if prop == "rp":
            return self._a*(1-self._e)
        if prop == "ra":
            return self._a*(1+self._e)
        if prop == "rad":
            # current radius
            E = self.get("EA", t)
            rad = self._a * (1 - self._e * math.cos(E))
            return rad
        if prop == "h":
            h = (self._mu * self._a * (1 - (self._e ** 2))) ** (.5)
            return h
        if prop == "ma":
            # calculate mean anomoly
            n = self.get("n")
            M = (n * (t - self._t0)) % (2 * math.pi)
            return M
        if prop == "va":
            # calculate true anomoly
            M = self.get("Ma", t)
            v = self.convertAnom("m", "v", M)
            return v
        if prop == "ea":
            # calculate eccentric anomoly
            M = self.get("Ma", t)
            E = self.convertAnom("m", "e", M)
            return E
        if prop == "r":
            v = self.get("va", t)
            rad = self.get("rad", t)
            r_o = Vec_3(cos(v), sin(v), 0) * rad
            return self._o_to_ECI(r_o)
        if prop == "vel":
            E = self.get("ea", t)
            rad = self.get("rad", t)
            vel_o = Vec_3(-sin(E), cos(E) * ((1 - self._e ** 2) ** .5), 0) * ((self._mu * self._a) ** .5) / rad
            return self._o_to_ECI(vel_o)
        if prop == "state":
            return {
                "r": self.get("r", t),
                "vel": self.get("vel", t)
            }
        print("Invalid property get.")
        return 0

    def _o_to_ECI(self, u):
        #Convert from orbital plane vector to ECI vector
        w = self._w
        omega = self._omega
        i = self._i
        # Matrix to convert from orbital plane to ECI
        R = [
            [cos(w)*cos(omega)-sin(w)*cos(i)*sin(omega), -sin(w)*cos(omega)-cos(w)*cos(i)*sin(omega), 0],
            [cos(w)*sin(omega)+sin(w)*cos(i)*cos(omega), cos(w)*cos(i)*cos(omega)-sin(w)*sin(omega), 0],
            [sin(w)*sin(i), cos(w)*sin(i), 0]
        ]
        return mat_multi(R, u)

    def convertAnom(self, typeIn, typeOut, anom):
        """
        Convert between mean, eccentric, and true anomalies
        """
        typeIn = typeIn.strip().lower()
        typeOut= typeOut.strip().lower()
        if typeIn == "v" and typeOut == "e":
            E = 2 * math.atan(math.tan(anom / 2) / ((1 + self._e) / (1 - self._e)) ** (0.5))
            return E
        elif typeIn == "v" and typeOut == "m":
            E = self.convertAnom("v", "e", anom)
            M = self.convertAnom("e", "m", E)
            return M
        elif typeIn == "e" and typeOut == "m":
            M = anom - self._e * math.sin(anom)
            return M
        elif typeIn == "e" and typeOut == "v":
            num = sin(anom) * (1 - self._e ** 2) ** .5
            den = cos(anom) - self._e
            v = math.atan2(num, den)
            if v < 0:
                return 2 * math.pi + v
            return v
        elif typeIn == "m" and typeOut == "e":
            # orbit is circular
            if self._e == 0:
                return anom
            # calculate eccentric anomaly with newton's method
            E0 = anom - 1
            E1 = anom
            while abs(E1 - E0) > .01:
                E0 = E1
                f = E0 - self._e * math.sin(E0) - anom
                f_ = 1 - self._e * math.cos(E0)
                if f_ == 0:
                    E1 = E0
                else:
                    E1 = E0 - f / f_
            E = E1
            return E
        elif typeIn == "m" and typeOut == "v":
            E = self.convertAnom("m", "e", anom)
            v = self.convertAnom("e", "v", E)
            return v
        print("Invalid anomoly conversion.")
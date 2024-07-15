import numpy as np
import math


# Problem 1
def calculate_buoyancy(v: float, density_fluid: float):
    """Return: Force of buoyancy
    v: Volume displaced by AUV
    density_fluid: Density of fluid AUV is submerged in"""

    g = 9.81
    buoyancy = density_fluid * v * g

    return buoyancy


# Problem 2
def will_it_float(v: float, mass: float):
    """Return: True if floats, False if sinks, None if nuetral
    v: Volume displaced by AUV
    mass: Mass of AUV"""

    g = 9.81
    if calculate_buoyancy(v, mass / v) > mass * g:
        return True
    elif calculate_buoyancy(v, mass / v) < mass * g:
        return False
    else:
        return None


# Problem 3
def calculate_pressure(depth: float):
    """Return: Pressure Force at a depth
    depth: Depth of AUV"""

    g = 9.81
    pressure = depth * g + 1
    # Add full pressure equation for air

    return pressure


# Problem 4
def calculate_acceleration(F: float, m: float):
    """Return: Acceleration of AUV
    F: Force exerted on AUV
    m: Mass of AUV"""

    acceleration = F / m

    return acceleration


# Problem 5
def calculate_angular_acceleration(tau: float, I: float):
    """Return: Angular acceleration of AUV
    tau: Torque exerted on AUV
    I: Moment of inertia"""

    acceleration = tau / I

    return acceleration


# Problem 6
def calculate_torque(F_magnitude: float, F_direction: float, r: float):
    """Return: torque exerted on AUV
    F_magnitude: Magnitude of force exerted on AUV
    F_direction: Direction of Force exerted on AUV
    r: Distance from axis of rotation to point of force"""

    torque = F_magnitude * r
    if F_direction < 0:
        torque *= -1

    return torque


# Problem 7
def calculate_moment_of_inertia(m: float, r: float):
    """Return: Momment of inertia of AUV
    m: Mass of AUV
    r: Distance from AUV axis of rotation to force applied"""

    moi = m * r**2

    return moi


# Problem 8
def calculate_auv_acceleration(
    F_magnitude: float, F_angle: float, mass=100, volume=0.1, thruster_distance=0.5
):
    """Return: Magnitude of acceleration of AUV
    F_magnitude: Magnitude of force exerted on AUV
    F_angle: Angle of Force applied
    mass: Mass of AUV
    volume: Volume of AUV
    thruster_distance: Distance of thruster to axis of rotation"""

    # acceleration = calculate_acceleration(F_magnitude, mass)
    force = F_magnitude * np.array([np.cos(F_angle), np.sin(F_angle)])
    acceleration = force / mass

    return acceleration


def calculate_auv_angular_acceleration(
    F_magnitude: float, F_angle: float, inertia=1, thruster_distance=0.5
):
    """Return: Magnitude of angular acceleration of AUV
    F_magnitude: Magnitude of force exerted on AUV
    F_angle: Angle of Force applied
    inertia: Moment of inertia of AUV
    thruster_distance: Distance of thruster to axis of rotation"""

    a_acceleration = calculate_angular_acceleration(F_magnitude, inertia)

    return a_acceleration


# Problem 9
def calculate_auv2_acceleration(T: np.ndarray, alpha: float, theta: float, mass=100):
    """Return: Array of magnitude of acceleration of AUV in x and y dimensions
    T: Array of magnitudes of forc applied by each thruster
    alpha: Angle of thrusters in radians
    theta: Angle of AUV in radians
    mass: Mass of AUV"""

    sin = math.sin(alpha * 180 / math.pi)
    x_force1 = T[1] * sin - T[0] * sin
    x_force2 = T[3] * sin - T[2] * sin
    x_net = x_force1 + x_force2

    cos = math.cos(alpha * 180 / math.pi)
    y_force1 = T[0] * cos - T[3] * cos
    y_force2 = T[1] * cos - T[2] * cos
    y_net = y_force1 + y_force2

    force = math.sqrt(x_net**2 + y_net**2)
    F_angle = math.atan(x_net / y_net)

    acceleration = calculate_acceleration(force, mass)

    return np.array([x_net, y_net]) / mass


def calculate_auv2_angular_acceleration(
    T: np.ndarray, alpha: float, L: float, l: float, inertia=100
):
    """Return: Magnitude of acceleration of AUV in x and y dimensions
    T: Array of magnitudes of forc applied by each thruster
    alpha: Angle of thrusters in radians
    L: Distance from the center of mass of the AUV to the thrusters
    l: Distance from the center of mass of the AUV to the thrusters
    inertia: Inertia of AUV"""

    sin = math.sin(alpha * 180 / math.pi)
    x_force1 = T[1] * sin - T[0] * sin
    x_force2 = T[3] * sin - T[2] * sin
    x_angular = (x_force1 - x_force2) * L

    cos = math.cos(alpha * 180 / math.pi)
    y_force1 = T[0] * cos - T[3] * cos
    y_force2 = T[1] * cos - T[2] * cos
    y_angular = (y_force1 - y_force2) * L

    force = y_angular + x_angular
    torque = force * L

    acceleration = calculate_angular_acceleration(torque, inertia)

    return acceleration


# Problem 10
def simulate_auv2_motion(
    T: np.ndarray,
    alpha: float,
    L: float,
    l: float,
    mass=100,
    inertia=100,
    dt=0.1,
    t_final=10,
    x0=0,
    y0=0,
    theta0=0,
):
    """Return: a series of 7 ndarrays for direct use by plot_auv2_motion()
    T: Array of magnitudes of forces applied by each thruster
    alpha: Angle of thrusters in radians
    L: Distance from the center of mass of the AUV to the thrusters
    l: Distance from the center of mass of the AUV to the thrusters
    mass: Mass of AUV
    inertia: Moment of inertia of AUV
    dt: Time step for t array
    t_final: Time to end simulation on
    x0: Initial x position
    y0: Initial y position
    theta0: Initial AUV angle"""

    t = [0]
    x = [x0]
    y = [y0]
    theta = [theta0]
    v = [0]
    omega = [0]
    a = [0]
    x_accel = 0
    y_accel = 0

    for i in range(1, t_final, dt):
        t.append(i)
        x_accel, y_accel = calculate_auv2_acceleration(T, alpha, theta, mass)
        a.append(math.sqrt(x_accel**2 + y_accel**2))
        v.append(a[i])
        x.append(x[i-1] + )
        y.append(x[i-1] + )
        
        calculate_auv2_angular_acceleration(T, alpha, L, l, inertia)

    return (
        np.ndarray(t),
        np.ndarray(x),
        np.ndarray(y),
        np.ndarray(theta),
        np.ndarray(v),
        np.ndarray(omega),
        np.ndarray(a),
    )


def plot_auv2_motion(
    t: np.ndarray,
    x: np.ndarray,
    y: np.ndarray,
    theta: np.ndarray,
    v: np.ndarray,
    omega: np.ndarray,
    a: np.ndarray,
):
    """ """
    return None

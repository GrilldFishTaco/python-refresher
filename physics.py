import numpy as np
import math


# Problem 1
def calculate_buoyancy(v, density_fluid):
    """ """
    g = 9.81
    bouyancy = density_fluid * v * g

    return buoyancy


# Problem 2
def will_it_float(v, mass):
    """ """
    g = 9.81
    if calculate_buoyancy(v, mass / v) > mass * g:
        return True
    elif calculate_buoyancy(v, mass / v) < mass * g:
        return False
    else:
        return None


# Problem 3
def calculate_pressure(depth):
    """ """
    g = 9.81
    pressure = depth * g

    return pressure


# Problem 4
def calculate_acceleration(F, m):
    """ """
    acceleration = F / m

    return acceleration


# Problem 5
def calculate_angular_acceleration(tau, I):
    """ """
    acceleration = tau / I

    return acceleration


# Problem 6
def calculate_torque(F_magnitude, F_direction, r):
    """ """
    torque = F_magnitude * r
    if F_direction < 0:
        torque *= -1

    return torque


# Problem 7
def calculate_moment_of_inertia(m, r):
    """ """
    moi = m * r**2

    return moi


# Problem 8
def calculate_auv_acceleration(
    F_magnitude, F_angle, mass=100, volume=0.1, thruster_distance=0.5
):
    """ """
    acceleration = calculate_acceleration(F_magnitude, mass)

    return acceleration


def calculate_auv_angular_acceleration(
    F_magnitude, F_angle, inertia=1, thruster_distance=0.5
):
    """ """
    a_acceleration = calculate_angular_acceleration(F_magnitude, inertia)

    return a_acceleration


# Problem 9
def calculate_auv2_acceleration(T, alpha, theta, mass=100):
    """ """
    sin = math.sin(alpha)
    x_force1 = T[1] * sin - T[0] * sin
    x_force2 = T[3] * sin - T[2] * sin
    x_net = x_force1 + x_force2

    cos = math.cos(alpha)
    y_force1 = T[0] * cos - T[3] * cos
    y_force2 = T[1] * cos - T[2] * cos
    y_net = y_force1 + y_force2

    force = math.sqrt(x_net**2 + y_net**2)
    F_angle = math.atan(x_net / y_net)

    acceleration = calculate_acceleration(force, mass)

    return acceleration, F_angle


def calculate_auv2_angular_acceleration(T, alpha, L, l, inertia=100):
    """ """
    sin = math.sin(alpha)
    x_force1 = T[1] * sin - T[0] * sin
    x_force2 = T[3] * sin - T[2] * sin
    x_angular = (x_force1 - x_force2) * L

    cos = math.cos(alpha)
    y_force1 = T[0] * cos - T[3] * cos
    y_force2 = T[1] * cos - T[2] * cos
    y_angular = (y_force1 - y_force2) * L

    force = y_angular + x_angular
    torque = force * L

    acceleration = calculate_angular_acceleration(torque, inertia)

    return acceleration


# Problem 10
def simulate_auv2_motion(
    T, alpha, L, l, mass=100, inertia=100, dt=0.1, t_final=10, x0=0, y0=0, theta0=0
):
    """ """
    t = [0]
    x = [x0]
    y = [y0]
    theta = [theta0]
    v = [0]
    omega = [0]
    a = [0]

    for i in range(1, t_final, dt):
        t.append(i)
        calculate_auv2_acceleration(T, alpha, theta, mass)
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


def plot_auv2_motion(t, x, y, theta, v, omega, a):
    """ """
    return None

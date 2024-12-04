import math

def lambda(x):
    d = 1/(1000 * 1000)
    m = 1
    theta = x * math.pi / 180
    theta_unc = 1

    return (d / m) * math.cos(theta) * theta_unc

Web VPython 3.2
        
# https://github.com/yellowstoneparks/sdic/
# Lyapanov exponent for the Lorenz equations
        
num_instances = min(int(input("How many instances would you like to spawn? (max 20): ")), 20)
        
# Defining some initial variables
sigma_range = (9.5, 10.5)
rho_range = (27, 29)
beta_range = (2.6, 2.8)
delta_t = 0.001
num_steps = 100 * 1000
d = 1e-5
lorenz = graph(title="Lorenz Coordinates")

parameters = []

def similarColor(v1, v2, thres=.9):
    return thres > sqrt((v1.x - v2.x)**2 +
                        (v1.y - v2.y)**2 +
                        (v1.z - v2.z)**2)

xs = vector(random(), random(), random())

ys = vector(random(), random(), random())
while not similarColor(ys, xs): ys = vector(random(), random(), random())

zs = vector(random(), random(), random())
while not similarColor(ys, zs) and not similarColor(xs, zs): zs = vector(random(), random(), random())

temp1 = gcurve(graph=lorenz, color=xs, label="X")
temp2 = gcurve(graph=lorenz, color=ys, label="Y")
temp3 = gcurve(graph=lorenz, color=zs, label="Z")
    
# parameters for each instance
for _ in range(num_instances):
    one_instance={
        "sigma": sigma_range[0] + random() * (sigma_range[1] - sigma_range[0]),
        "rho": rho_range[0] + random() * (rho_range[1] - rho_range[0]),
        "beta": beta_range[0] + random() * (beta_range[1] - beta_range[0]),
        "x": 1.0,
        "y": 1.0,
        "z": 1.0,
        "x_prime": 1.0 + d,
        "y_prime": 1.0 + d,
        "z_prime": 1.0 + d,
        "diff_0": sqrt(3 * d**2),
        "lambda": 0.0,
        "curve_x": gcurve(graph=lorenz, color=xs),
        "curve_y": gcurve(graph=lorenz, color=ys),
        "curve_z": gcurve(graph=lorenz, color=zs),
        "path": sphere(pos=vector(1.0, 1.0, 1.0), radius=0.1, color=vector(random(), random(), random()), make_trail=True, trail_radius=0.1)
    }
    parameters.append(one_instance)
    
for i in range(num_steps):
    rate(9000)
    for param in parameters:
        sigma = param["sigma"]
        rho = param["rho"]
        beta = param["beta"]
        x, y, z = param["x"], param["y"], param["z"]
        x_prime, y_prime, z_prime = param["x_prime"], param["y_prime"], param["z_prime"]
        diff_0 = param["diff_0"]

        # Calculating the d_/dt for each coord
        dxdt = sigma * (y - x)
        dydt = x * (rho - z) - y
        dzdt = x * y - beta * z

        dxdt_prime = sigma * (y_prime - x_prime)
        dydt_prime = x_prime * (rho - z_prime) - y_prime
        dzdt_prime = x_prime * y_prime - beta * z_prime

        # updating x y z and prime x y z
        param["x"] += dxdt * delta_t
        param["y"] += dydt * delta_t
        param["z"] += dzdt * delta_t

        param["x_prime"] += dxdt_prime * delta_t
        param["y_prime"] += dydt_prime * delta_t
        param["z_prime"] += dzdt_prime * delta_t

        # Calculating new diff
        diff = sqrt((param["x"] - param["x_prime"])**2 +
                    (param["y"] - param["y_prime"])**2 +
                    (param["z"] - param["z_prime"])**2)
        if diff == 0:
            diff = 1e-10  # avoid log(0)

        param["lambda"] += log(diff / diff_0)

        # rescaling so it all fits
        scaling_factor = d / diff
        param["x_prime"] = param["x"] + (param["x_prime"] - param["x"]) * scaling_factor
        param["y_prime"] = param["y"] + (param["y_prime"] - param["y"]) * scaling_factor
        param["z_prime"] = param["z"] + (param["z_prime"] - param["z"]) * scaling_factor

        param["diff_0"] = d
        
        param["path"].pos = vector(param["x"], param["y"], param["z"])

        # plotting
        param["curve_x"].plot(i, param["x"])
        param["curve_y"].plot(i, param["y"])
        param["curve_z"].plot(i, param["z"])

# getting Lyapunov exponents for each
total_time = num_steps * delta_t
for idx, param in enumerate(parameters):
    param["lambda"] /= total_time
    print("Instance " + str(idx + 1) + " - Sigma: " + str(round(param['sigma'], 2)) +
      ", Rho: " + str(round(param['rho'], 2)) +
      ", Beta: " + str(round(param['beta'], 2)) +
      ", Lyapunov Exponent: " + str(round(param['lambda'], 4)))


### END OF LORENZ

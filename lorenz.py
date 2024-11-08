Web VPython 3.2
        
# Lyapanov exponent for the Lorenz equations
sigma = 10
rho = 28
beta = 8 / 3

x, y, z = 1.0, 1.0, 1.0
delta_t = 0.001
num_steps = 200000

d = 1e-5
x_prime, y_prime, z_prime = x + d, y + d, z + d

diff_0 = sqrt((x - x_prime)**2 + (y - y_prime)**2 + (z - z_prime)**2)

path = sphere(pos=vector(x, y, z), radius=0.1, make_trail=True, trail_radius=0.1)
lorenz = graph(title="Lorenz Coordinates")
lorenz_x_curve = gcurve(graph=lorenz, color=color.green, label="X")
lorenz_y_curve = gcurve(graph=lorenz, color=color.blue, label="Y")
lorenz_z_curve = gcurve(graph=lorenz, color=color.red, label="Z")

# sum
lamb = 0.0

for i in range(num_steps):
    rate(8000)
    
    dxdt = sigma * (y - x)
    dydt = x * (rho - z) - y
    dzdt = x * y - beta * z
    
    x += dxdt * delta_t
    y += dydt * delta_t
    z += dzdt * delta_t
    
    dxdt_prime = sigma * (y_prime - x_prime)
    dydt_prime = x_prime * (rho - z_prime) - y_prime
    dzdt_prime = x_prime * y_prime - beta * z_prime
    
    x_prime += dxdt_prime * delta_t
    y_prime += dydt_prime * delta_t
    z_prime += dzdt_prime * delta_t
    
    diff = sqrt((x - x_prime)**2 + (y - y_prime)**2 + (z - z_prime)**2)
    
    if diff == 0:
        diff = 1e-10  # avoid log(0)
        
    lamb += log(diff / diff_0)
    
    # Rescale
    scaling_factor = d / diff
    x_prime = x + (x_prime - x) * scaling_factor
    y_prime = y + (y_prime - y) * scaling_factor
    z_prime = z + (z_prime - z) * scaling_factor
    
    diff_0 = d  # we've rescaled 'diff' to 'd'
    
    path.pos = vector(x, y, z)
    
    lorenz_x_curve.plot(i, x)
    lorenz_y_curve.plot(i, y)
    lorenz_z_curve.plot(i, z)

total_time = num_steps * delta_t
lamb /= total_time
print("Lyapunov Exponent:", lamb)


### END OF LORENZ

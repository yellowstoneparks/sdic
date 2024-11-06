Web VPython 3.2

# x n+1 = r*xn*(1-xn)
# f'(xn)=r*(1−2*xn)

# write code here
x0 = 0.448
xn = x0
n = 1000
delta_r = 0.001

r_vs_lamb = graph(title="r vs lyapunov exponent", xtitle="r", ytitle="lyapunov exponent", ymin=-2, ymax=1, xmin=0)
exponents_pos = gdots(graph=r_vs_lamb, size=.5, color=color.green)
exponents_neg = gdots(graph=r_vs_lamb, size=.5, color=color.blue)

for r in range(0, 4.001, delta_r):
    def der(xn):
        return r * (1 - 2 * xn)
    
    lamb = 0
    xn = x0
    
    for i in range(0, n):
        xn = r * xn * (1 - xn)
        lamb += log(abs(der(xn))) 
    
    lamb /= n
    
    if lamb > 0:
        exponents_pos.plot(r, lamb)
    else:
        exponents_neg.plot(r, lamb)


# Lyapanov exponent for the Lorenz equations
sigma = 8.5
rho = 24
beta = 8 / 3

x, y, z, delta_t, num_steps = 1, 1, 1, 0.01, 10000

path = sphere(pos=vector(x,y,z), radius=.1, make_trail=True, trail_radius=.1)
lorenz = graph(title="Lorenz Coordinates")
lorenz_x_curve = gcurve(graph=lorenz, color=color.green, label="X")
lorenz_y_curve = gcurve(graph=lorenz, color=color.blue, label="Y")
lorenz_z_curve = gcurve(graph=lorenz, color=color.red, label="Z")

lamb_x, lamb_y, lamb_z = 0, 0, 0

for i in range(num_steps):
    rate(50)
    
    dxdt = sigma * (y - x)
    dydt = x * (rho - z) - y
    dzdt = x * y - beta * z
    
    x += dxdt * delta_t
    y += dydt * delta_t
    z += dzdt * delta_t
    
    lamb_x += log(abs(dxdt))
    lamb_y += log(abs(dydt))
    lamb_z += log(abs(dzdt))
    
    path.pos = vector(x, y, z)
    
    lorenz_x_curve.plot(i, x)
    lorenz_y_curve.plot(i, y)
    lorenz_z_curve.plot(i, z)


lamb_x /= num_steps
lamb_y /= num_steps
lamb_z /= num_steps

print(lamb_x, lamb_y, lamb_z)

### END OF LORENZ

##simple harmonic oscillator
##plan: want to graph seperation (or the equations of the two functions) over time and on the same graph and print the lyapunov exponent


simple_harmonic_oscillator = graph(title = "Lyapunov Exponent of a Simple Harmonic Osicllator", xtitle = "time", ytitle = "xcor")


oscillator_1 = 

##feigenbaum ratio 

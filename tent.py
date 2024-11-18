Web VPython 3.2

# NOTE: I changed mew to r because I didn't want to use special characters

# variables
ITERATION = 100
SKIP = 100 # for transients
k0 = 40000 # base value of r (aka control parameter)
m = 9999999 # divisor used to initalize x0
u = 2 # constant

# graphs
r_vs_lamb = graph(title="R vs Lyapunov exponent", xtitle="R", ytitle="Lyapunov exponent",xmin=0,xmax=4,ymin=-6,ymax=1)
tentmap = graph(title="Tent Map", xtitle="x_t", ytitle="x_t+1",xmin=0,xmax=1,ymin=0,ymax=1)
tent = gdots()
exponents_pos = gdots(graph = r_vs_lamb, size=2, color=color.orange, label = "Positive Lyapunov Exponent")
exponents_neg = gdots(graph = r_vs_lamb, size=2, color=color.blue, label = "Negative Lyapunov Exponent")

# parameter settings
start = 0 
end = 4
delta_r = 0.001

#initial conditions
k = u * k0
x0 = k / m
x = x0

# plot tent map
for t in range(SKIP):
    k = m - abs(m - u * k) # update k with tent map
    x = k / m # normalize k to update x

lyap = 0.0
for t in range(0,ITERATION,delta_r):
    pre_x = x
    k = m - abs(m - u * k)
    x = k / m
    
    tent.plot(pre_x,x)

# compute lyapunov exponents for different values of r
for r in range(start, end, delta_r):
    k = r * k0
    x0 = r / m
    x = x0

    # transients, used to stabilize system
    for t in range(SKIP):
        k = m - abs(m - r * k)
        x = k / m

    lyap = 0.0
    for t in range(ITERATION):
        k = m - abs(m - r * k) # update k with tent map
        x = k / m # update x with normalized k

        temp = r * (1.0 - r * x) / abs(1.0 - r * x) # dererative of tent map
        current_lyap = log(abs(temp)) # apply lyapunov equation
        lyap += current_lyap # sum up lyapunov
        
    lyap /= ITERATION # average lyapunov exponent

    # plot the lyapunov exponent based sign
    if lyap > 0:
        exponents_pos.plot(r, lyap) # positive lyapunov exponent
    else:
        exponents_neg.plot(r, lyap) # negative lyapunov exponent

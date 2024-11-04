Web VPython 3.2

# x n+1 = r*xn*(1-xn)
# f'(xn)=r*(1−2*xn)

# write code here
x0 = 0.3348
xn = x0
n = 1000
delta_r = 0.01

def der(xn):
    return r * (1 - 2 * xn)

r_vs_lamb = graph(title="r vs lyapunov exponent", xtitle="r", ytitle="lyapunov exponent", ymin=-10, ymax=10)
exponents = gcurve(graph=r_vs_lamb)

for r in range(0, 7, delta_r):
    lamb = 0
    
    for i in range(0, n):
        xn = r * xn * (1 - xn)
        lamb += log(abs(der(xn))) 
    
    lamb /= n
    exponents.plot(r, lamb)

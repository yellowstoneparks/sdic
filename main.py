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

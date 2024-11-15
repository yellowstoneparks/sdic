Web VPython 3.2

ITERATION = 100
SKIP = 100

k0 = 43597
m = 9999991
k = k0

r_vs_lamb = graph(title="r vs lyapunov exponent", xtitle="r", ytitle="lyapunov exponent",xmin=0,xmax=4,ymin=-6,ymax=1)
tentmap = graph(title="tent map", xtitle="r", ytitle="lyapunov exponent",xmin=0,xmax=1,ymin=0,ymax=1)
tent = gdots()
exponents_pos = gdots(graph = r_vs_lamb, size=2, color=color.orange, label = "Positive Lyapunov Exponent")
exponents_neg = gdots(graph = r_vs_lamb, size=2, color=color.blue, label = "Negative Lyapunov Exponent")

start = 0
end = 4
delta_r = 0.001

for A in range(start, end, delta_r):
    k = A * k0

    x0 = k / m
    x = x0

    for t in range(SKIP):
        k = m - abs(m - A * k)
        x = k / m

    lyap = 0.0
    for t in range(ITERATION):
        pre_x = x
        k = m - abs(m - A * k)
        x = k / m

        temp = A * (1.0 - A * x) / abs(1.0 - A * x)
        current_lyap = log(abs(temp))
        lyap += current_lyap
        
        tent.plot(A,x)
        
    lyap /= ITERATION

    if lyap > 0:
        exponents_pos.plot(A, lyap)
    else:
        exponents_neg.plot(A, lyap)

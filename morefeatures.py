Web VPython 3.2

# write code here
x0 = 0.448
xn = x0
n = 1000
delta_r = 0.001

startofisland = None
startlist = []
endofisland = None
endlist = []
lastlamb = None

# ymin=-2, ymax=1, 
r_vs_lamb = None
exponents_pos = None
exponents_neg = None
# points of no return, when it switches from non-chaotic to chaotic and vice versa
ponr = None

def calclyapunov(start,end,Ymin,Ymax,iscurve):
    global startlist, endlist, lastlamb, exponents_pos, exponents_neg, r_vs_lamb, ponr
#    r_vs_lamb = graph(title="r vs lyapunov exponent", xtitle="r", ytitle="lyapunov exponent", ymin=Ymin, ymax=Ymax, xmin=start, xmax=end)
    r_vs_lamb = graph(title="r vs lyapunov exponent", xtitle="r", ytitle="lyapunov exponent", fast=False, width=800, ymin=Ymin, ymax=Ymax, xmin=start,xmax=end)
    if iscurve == 0:
        exponents_pos = gdots(graph = r_vs_lamb, size=2, color=color.green)
        exponents_neg = gdots(graph = r_vs_lamb, size=2, color=color.blue)
        ponr = gdots(graph = r_vs_lamb, size = 3, color = color.red, label = "Point of No Return")
    else:
        exponents_pos = gdots(graph = r_vs_lamb, size=5, color=color.green)
        exponents_neg = gdots(graph = r_vs_lamb, size=5, color=color.blue)
        ponr = gdots(graph = r_vs_lamb, size = 7, color = color.red, label = "Point of No Return")
    for r in range(start, end, delta_r):
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
#         island of stability appears when lyapunov expnent becomes negative
        if lastlamb is not None:
            if lamb < 0 and lastlamb > 0:
                endofisland = round(r - delta_r,3)
                startlist.append(endofisland)
                ponr.plot(endofisland,lamb)
            if lamb > 0 and lastlamb < 0:
                startofisland = round(r,3)
                endlist.append(startofisland)
                ponr.plot(startofisland,lamb)
        lastlamb = lamb
    
calclyapunov(0,4,-1,1,0)
endlist = endlist[1:]
print("Starts of Island of Stability: ",startlist)
print("Ends of Island of Stability: ",endlist)

resetButton = button(bind = reset, text = "Click to examine island of stability", pos=scene.title_anchor)
index = 0

def reset(evt):
    global ponr, exponents_pos, exponents_neg, r_vs_lamb, index
    if startlist and endlist and index < len(startlist):
        print("Islands of Stability - Start:", startlist[index], "End:", endlist[index])
        calclyapunov(startlist[index]-0.01, endlist[index]+0.01,-0.5,0.5,1)
    else:
        print("End of Islands of Stability")
    index += 1

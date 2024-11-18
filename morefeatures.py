Web VPython 3.2

# basic variables
x0 = 0.448
n = 1000
delta_r = 0.001

# island of stability variables
startlist = []
endlist = []
lastlamb = None

# graph elements 
r_vs_lamb = None
exponents_pos = None
exponents_neg = None
ponr = None # points of no return, when it switches from non-chaotic to chaotic and vice versa

# func to calc and plot lyapunov exp over a range of r values
def calclyapunov(start,end,Ymin,Ymax,iscurve):
    global startlist, endlist, lastlamb, exponents_pos, exponents_neg, r_vs_lamb, ponr
    # initialize graph + configure based on curve style
    r_vs_lamb = graph(title="r vs lyapunov exponent", xtitle="r", ytitle="lyapunov exponent", fast=False, width=800, ymin=Ymin, ymax=Ymax, xmin=start,xmax=end)
    if iscurve == 0:
        exponents_pos = gdots(graph = r_vs_lamb, size=2, color=color.orange, label = "Positive Lyapunov Exponent")
        exponents_neg = gdots(graph = r_vs_lamb, size=2, color=color.blue, label = "Negative Lyapunov Exponent")
        ponr = gdots(graph = r_vs_lamb, size = 3, color = color.red, label = "Point of No Return")
    else:
        exponents_pos = gdots(graph = r_vs_lamb, size=12, color=color.orange, label = "Positive Lyapunov Exponent")
        exponents_neg = gdots(graph = r_vs_lamb, size=12, color=color.blue, label = "Negative Lyapunov Exponent")
        ponr = gdots(graph = r_vs_lamb, size = 15, color = color.red, label = "Point of No Return")
    # loop to calc lyapunov
    for r in range(start, end, delta_r):
        def der(xn):
            return r * (1 - 2 * xn)
        lamb = 0
        xn = x0
        for i in range(0, n):
            xn = r * xn * (1 - xn)
            lamb += log(abs(der(xn))) 
        lamb /= n
        # plot lyapunov based on sign
        if lamb > 0:
            exponents_pos.plot(r, lamb)
        else:
            exponents_neg.plot(r, lamb)
        # find islands of stability
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

# initial lyapunov calc
calclyapunov(0,4,-1,1,0)

endlist = endlist[1:] # had to remove first since it was faulty
print("Starts of Island of Stability: ",startlist)
print("Ends of Island of Stability: ",endlist)

# button to zoom in on islands of stability
islandButton = button(bind = island, text = "Click to examine island of stability", pos=scene.title_anchor)
index = 0
scene.delete()

def island(evt):
    global ponr, exponents_pos, exponents_neg, r_vs_lamb, index
    if startlist and endlist and index < len(startlist):
        print("Islands of Stability - Start:", startlist[index], "End:", endlist[index])
        calclyapunov(startlist[index]-0.01, endlist[index]+0.01,-0.5,0.5,1) # extra space of 0.01 to allow room to see
    else:
        print("End of Islands of Stability")
    index += 1

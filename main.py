# x n+1 = r*xn*(1-xn)
# f'(xn)=r*(1−2*xn)

# write code here
x0 = 0.3348
xn = x0
n = 1000
r = 3.2

def der(xn):
    return r * (1 - 2 * xn)

lamb = 0

for i in range(0, n):
    xn = r * xn * (1 - xn)
    lamb += log(abs(der(xn))) 

lamb = lamb / n
print(lamb)

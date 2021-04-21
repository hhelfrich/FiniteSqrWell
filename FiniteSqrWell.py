import numpy as np
import matplotlib.pyplot as plt

c = 1E-10 #m
V_0 = 8E-15 #J
m = 9.109E-31
hbar = 1.0546E-34 #J s
z_0 = (c/hbar)*np.sqrt(2*m*V_0)
max_iter = 1000
def find_pairs(f, step, a, b):
    x = a
    pairs = []
    while (x + step < b):
        if (f(x+step)/f(x) < 0):
             pairs.append([x, x+step])
        x += step
    return pairs

def bisection(f, pairs, tolerance, max_iter):
    zeros = []
    for pair in pairs:
        mid = (pair[1]-pair[0])/2 + pair[0]
        iter = 1
        while (abs(f(mid)) > tolerance and iter < max_iter):
            if (f(mid)/f(pair[0]) <0): pair[1] = mid
            else: pair[0] = mid
            mid = (pair[1]-pair[0])/2 + pair[0]
            iter += 1
        if (iter < 1000):
            zeros.append(mid)
    return zeros

def symmetric(x):
    if (x ==0):
        return 1
    else:
        return np.sqrt((z_0/x)**2-1) - np.tan(x)

def antisymmetric(x):
    if (x == 0):
        return 1
    else:
        return np.sqrt((z_0/x)**2-1) + np.cos(x)/np.sin(x)

x = np.linspace(0,10,1000)
y = np.tan(x)
plt.ylim(0,8)
y2 = np.sqrt((z_0/x)**2-1)
plt.plot(x,y)
plt.plot(x,y2)
plt.savefig("symmetric.pdf")


plt.clf()

y3 = np.cos(x)/np.sin(x)

plt.plot(x, y2)
plt.plot(x, y3)
plt.savefig("anti.pdf")

y4 = np.sqrt((z_0/x)**2-1) - np.tan(x)
y5 = np.sqrt((z_0/x)**2-1) - (np.cos(x)/np.sin(x))

plt.clf()
x1 = np.linspace(-8,8,1000)
plt.ylim(-8,8)
plt.plot(x1, y4)


pairs = find_pairs(symmetric, 0.1, 0, 10)
print(pairs)
zeros = bisection(symmetric, pairs, 1E-10, 1000)
print(zeros)

Energies = []
for z in zeros:
    Energies.append(hbar**2/(2*m*c**2)*z**2 - V_0)
print(Energies)

pairs = find_pairs(antisymmetric, 0.1, 0, 10)
print(pairs)
zeros = bisection(antisymmetric, pairs, 1E-10, 1000)
print(zeros)
for z in zeros:
    Energies.append(hbar**2/(2*m*c**2)*z**2 - V_0)
print(Energies)

def Infinite(n):
    return n**2*np.pi**2*hbar**2/(8*m*c**2) - V_0

print("infinite", Infinite(1), Infinite(2), Infinite(3), Infinite(4), Infinite(5))
# Joel Surfleet's Isentropic Toolbox

def bisection(f,a,b,TOL):

    for all in range(100):
        if f(a)/f(b) < TOL:
            break
        elif abs(f(a)) > abs(f(b)):
            a = (a + b) / 2
        else:
            b = (a + b) / 2

    c = (a + b) / 2
    
    return c

class isentropic:
    def __init__(fluid,k):
        fluid.k = k

    def M(fluid,x):
        fluid.Tratio = 1 / ( 1 + (((fluid.k-1)/2) * (x**2)))
        fluid.T(fluid.Tratio)
        
    def P(fluid,x):
        fluid.Pratio   = x
        fluid.rhoratio = x ** (1/fluid.k)
        fluid.Tratio   = x ** ((fluid.k-1)/fluid.k)
        
    def rho(fluid,x):
        fluid.Pratio   = x ** (fluid.k)
        fluid.rhoratio = x
        fluid.Tratio   = x ** (fluid.k-1)
        
    def T(fluid,x):
        fluid.Pratio   = x ** (fluid.k/(fluid.k-1))
        fluid.rhoratio = x ** (1/(fluid.k-1))
        fluid.Tratio   = x

    def static(fluid,P,rho,T):
        fluid.Ps   = P
        fluid.rhos = rho
        fluid.Ts   = T

        fluid.P0   = P   / fluid.Pratio
        fluid.rho0 = rho / fluid.rhoratio
        fluid.T0   = T   / fluid.Tratio

        fluid.a = (fluid.k * P / rho) ** (1/2)

    def stag(fluid,P0,rho0,T0):
        fluid.P0   = P0
        fluid.rho0 = rho0
        fluid.T0   = T0

        fluid.Ps   = P0   * fluid.Pratio
        fluid.rhos = rho0 * fluid.rhoratio
        fluid.Ts   = T0   * fluid.Tratio

class normalShock:
    def __init__(fluid,M,k):
        fluid.k = k
        fluid.rhoratio = ((fluid.k+1) * M ** 2) / (2 + (fluid.k - 1) * M ** 2)
        fluid.Pratio = 1 + (2 * fluid.k) / (fluid.k + 1) * (M ** 2 - 1)
        fluid.Tratio = fluid.Pratio * fluid.rhoratio
    
class rayleigh:
    def __init__(flow,k):
        flow.k = k

    def M(flow,M):
        M = M ** 2

        flow.Pratio = (1 + flow.k) / (1 + flow.k * M)
        flow.Tratio = M * ((1 + flow.k) / (1 + flow.k * M)) ** 2
        flow.rhoratio = 1 / M * ((1 + flow.k * M) / (1 + flow.k))
        flow.T0ratio = (((flow.k + 1) * M) / ((1 + flow.k * M) ** 2)) * (2 + (flow.k - 1) * M)
        flow.P0ratio = ((1 + flow.k) / (1 + flow.k * M)) * (((2 + (flow.k - 1) * M) / (flow.k + 1)) ** (flow.k / (flow.k - 1)))

    def T0(flow,x):
        
        
        
        return flow.M2


# air = isentropic(1.4)
# air.M(0.2)
# air.static(1,1,273)
# print(air.P0,air.T0)

# flow = rayleigh(1.4)
# flow.Mstar(0.2)

# P02 = 1.028

# print()
# print(flow.Pratio)
# print(flow.Tratio)
# print(flow.P0ratio)
# print(flow.T0ratio)

# T02ratio = (1270)/(275.2) * flow.T0ratio

# print(T02ratio)

# Example work from aero 303 HW #2 problem 3
air = isentropic(1.4)
air.M(3)

print()
print("Pressure ratio",air.Pratio)
print("Density Ratio",air.rhoratio)
print("Temperature Ratio",air.Tratio)
print()

air.static(101000,1.4077,250)

print("Stagnation Pressure",air.P0/1000,"kPa")
print("Stagnation Density",air.rho0,"kg/m^3")
print("Stagnation Temperature",air.T0,"K")
print()

print("Speed of Sound:",air.a,"m/s")
print("Air Speed:",air.a * 3,"m/s")
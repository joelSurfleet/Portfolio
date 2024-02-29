# Joel Surfleet's Isentropic Toolbox

from math import log,sin,radians

gamma = 1.4

class isentropic:
    def __init__(fluid,k):
        fluid.k = k

    def M(fluid,M):
        fluid.Tratio = 1 / (1 + (((fluid.k - 1) / 2) * (M ** 2)))
        fluid.T(fluid.Tratio)
        
    def P(fluid,x):
        fluid.Pratio = x
        fluid.Dratio = x ** (1/fluid.k)
        fluid.Tratio = x ** ((fluid.k-1)/fluid.k)
        
    def D(fluid,x):
        fluid.Pratio = x ** (fluid.k)
        fluid.Dratio = x
        fluid.Tratio = x ** (fluid.k-1)
        
    def T(fluid,x):
        fluid.Pratio = x ** (fluid.k/(fluid.k-1))
        fluid.Dratio = x ** (1/(fluid.k-1))
        fluid.Tratio = x

    def getStag(fluid,P,D,T):
        fluid.Ps = P
        fluid.Ds = D
        fluid.Ts = T

        fluid.P0 = P / fluid.Pratio
        fluid.D0 = D / fluid.Dratio
        fluid.T0 = T / fluid.Tratio

        fluid.a = (fluid.k * P / D) ** (1/2)

    def getStatic(fluid,P0,D0,T0):
        fluid.P0 = P0
        fluid.D0 = D0
        fluid.T0 = T0

        fluid.Ps = P0 * fluid.Pratio
        fluid.Ds = D0 * fluid.Dratio
        fluid.Ts = T0 * fluid.Tratio
    
class rayleigh:
    def __init__(flow,k):
        flow.k = k

    def M(flow,M):
        M = M ** 2

        flow.Pratio = (1 + flow.k) / (1 + flow.k * M)
        flow.Tratio = M * ((1 + flow.k) / (1 + flow.k * M)) ** 2
        flow.Dratio = 1 / M * ((1 + flow.k * M) / (1 + flow.k))
        flow.T0ratio = (((flow.k + 1) * M) / ((1 + flow.k * M) ** 2)) * (2 + (flow.k - 1) * M)
        flow.P0ratio = ((1 + flow.k) / (1 + flow.k * M)) * (((2 + (flow.k - 1) * M) / (flow.k + 1)) ** (flow.k / (flow.k - 1)))

    def getStar(flow,P,D,T,P0,T0):
        flow.Pstar = P / flow.Pratio
        flow.Tstar = T / flow.Tratio
        flow.Dstar = D / flow.Dratio
        flow.P0star = P0 / flow.P0ratio
        flow.T0star = T0 / flow.T0ratio

        flow.P = P
        flow.D = D
        flow.T = T
        flow.P0 = P0
        flow.T0 = T0

    def getReal(flow,Pstar,Dstar,Tstar,P0star,T0star):
        flow.P = Pstar * flow.Pratio
        flow.D = Dstar * flow.Dratio
        flow.T = Tstar * flow.Tratio
        flow.P0 = P0star * flow.P0ratio
        flow.T0 = T0star * flow.T0ratio

        flow.Pstar = Pstar
        flow.Tstar = Tstar
        flow.Dstar = Dstar
        flow.P0star = P0star
        flow.T0star = T0star

class fanno:
    def __init__(flow,k):
        flow.k = k

    def M(flow,M):
        M = M ** 2
        k = flow.k
        flow.fL4D = (1 - M) / (k * M) + (k + 1) /(2 * k) * log(((k + 1) * M) / (2 + (k - 1)*M))

        flow.Pratio = (1 / M ** (1/2)) * ((k + 1) / (2 + (k - 1) * M)) ** (1/2)
        flow.Dratio = (1 / M ** (1/2)) * ((2 + (k - 1) * M) / (k + 1)) ** (1/2)
        flow.Tratio = (k + 1) / (2 + (k - 1) * M)
        flow.P0ratio = (1 / M ** (1/2)) * ((2 + (k - 1) * M) / (k + 1)) ** ((k + 1) / (2 * (k-1)))

    def getStar(flow,P,D,T,P0):
        flow.Pstar = P / flow.Pratio
        flow.Tstar = T / flow.Tratio
        flow.Dstar = D / flow.Dratio
        flow.P0star = P0 / flow.P0ratio

        flow.P = P
        flow.D = D
        flow.T = T
        flow.P0 = P0

    def getReal(flow,Pstar,Dstar,Tstar,P0star):
        flow.P = Pstar * flow.Pratio
        flow.D = Dstar * flow.Dratio
        flow.T = Tstar * flow.Tratio
        flow.P0 = P0star * flow.P0ratio

        flow.Pstar = Pstar
        flow.Tstar = Tstar
        flow.Dstar = Dstar
        flow.P0star = P0star

class normal:
    def __init__(shock,k):
        shock.k = k

    def M(shock,M):
        shock.Dratio = ((shock.k+1) * M ** 2) / (2 + (shock.k - 1) * M ** 2)
        shock.Pratio = 1 + (((2 * shock.k) / (shock.k + 1)) * ((M ** 2) - 1))
        shock.Tratio = shock.Pratio * shock.Dratio
        shock.M2 = (((shock.k - 1) * M ** 2 + 2) / (2 * shock.k * M ** 2 - (shock.k - 1))) ** (1/2)

    def get1(shock,P2,D2,T2):
        shock.P1 = P2 / shock.Pratio
        shock.T1 = T2 / shock.Tratio
        shock.D1 = D2 / shock.Dratio

        shock.P2 = P2
        shock.D2 = D2
        shock.T2 = T2

    def get2(shock,P1,D1,T1):
        shock.P2 = P1 * shock.Pratio
        shock.D2 = D1 * shock.Dratio
        shock.T2 = T1 * shock.Tratio

        shock.P1 = P1
        shock.T1 = T1
        shock.D1 = D1

class oblique():
    def __init__(shock,k):
        shock.k = k
        shock.Normal = normal(k)

    def getBeta(shock,theta,M):
        if theta < 0 or theta > 45 or M < 1 or M > 21:
            print("My chart is NOT good enough for that")
            pass

        f = open(r"C:\Users\joels\Documents\Python\ThetaMachBeta.txt",'r')

        dataString = f.read()
        dataString = dataString.splitlines()        
        f.close()

        for i in range(360):
            dataString[i] = dataString[i].split()[0:3]
            for j in range(3):
                dataString[i][j] = float(dataString[i][j])

        for j in range(0,360,40):
            if dataString[j][0] == theta:
                print(dataString[j][0])
                for k in range(0,40):
                    if dataString[j+k][1] > M:
                        print(dataString[j+k][1])
                        if dataString[j+k][2] == '-1':
                            print('detached wave')
                            pass
                        shock.b = linInterp(M,dataString[j+k-1][1],dataString[j+k][1],dataString[j+k-1][2],dataString[j+k][2])
                        return shock.b
            
    def M(shock,M):
        shock.Mn = M * sin(radians(shock.b))
        shock.Normal.M(shock.Mn)

def linInterp(x1,y1,y2,z1,z2):
    # x1: your data point
    # y1, y2: the values bracketing your data point
    # z1, z2: the values bracketing your answer
    return (x1 - y1) / (y2 - y1) * (z2 - z1) + z1

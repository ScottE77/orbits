import math
import numpy as np
import matplotlib.pyplot as plt
#import sympy as sp

G = 6.67408*(10**(-11))
M = 1/G
a = 1
e = 0.5
b = a*math.sqrt(1-e**2)
x = []
y = []
for i in range(801):
    E = e*math.pi*(i/200)
    x.append(math.cos(E)+e)
    y.append(b*math.sin(E))

def FEstep(r0,v0,h): #this is the function I will loop
    v1 = [0,0] #initial definitions
    r1 = [0,0] #initial definitions
    r1[0] = r0[0] +h*v0[0]
    r1[1] = r0[1] +h*v0[1]
    r = math.sqrt((r0[0]**2)+(r0[1]**2))
    v1[0] =v0[0]+h*((-1)*(r0[0])/r**3)
    v1[1] =v0[1]+h*((-1)*(r0[1])/r**3)
    return [r1,v1]

def FE(N,FEPO,r0,v0):#FEPO = force evaluations per orbit
    h = (2*math.pi)/FEPO
    allx0 = []
    ally0 = []
    r = math.sqrt((r0[0]**2)+(r0[1]**2))
    E0 = (1/2)*(v0[0]**2+v0[1]**2)-(1/r)
    allE = [E0]
    allEerror = [0]
    L0 = r0[0]*v0[1]-r0[1]*v0[0]
    allL = [L0]
    allLerror = [0]
    for i in range(N*FEPO):
        allx0.append(r0[0])
        ally0.append(r0[1])
        r = math.sqrt((r0[0]**2)+(r0[1]**2))
        E1 = (1/2)*(v0[0]**2+v0[1]**2)-(1/r)
        allE.append(E1)
        allEerror.append(abs(E1-E0)/abs(E0))
        L1 = r0[0]*v0[1]-r0[1]*v0[0]
        allL.append(L1)
        allLerror.append(abs(L1-L0)/abs(L0))
        q = FEstep(r0,v0,h)
        r0 = q[0]
        v0 = q[1]
        
    t = []
    for i in range(N*FEPO+1):
        t.append(i)  
    
    x = []
    y = []
    for i in range(801):
        E = e*math.pi*(i/200)
        x.append(math.cos(E)+e)
        y.append(b*math.sin(E))
    return [allEerror,allLerror,t,allx0,ally0,r0,v0]

def FEgraph(N,FEPO):
    values = FE(N,FEPO,[a*(1+e),0],[0,math.sqrt((a*(1-(e)))/(a*(1+(e))))])
    plt.figure(1)
    plt.subplot(211)
    plt.plot(values[3],values[4])
    plt.plot(x,y)
    plt.plot(0,0,'r*')
    plt.figure(2)
    plt.loglog(values[2],values[0])
    plt.figure(3)
    plt.loglog(values[2],values[1])
    
def MEstep(r0,v0,h):
    v1 = [0,0] #initial definitions
    r1 = [0,0] #initial definitions
    r1[0] = r0[0] +h*v0[0]
    r1[1] = r0[1] +h*v0[1]
    r = math.sqrt((r1[0]**2)+(r1[1]**2))
    v1[0] =v0[0]+h*((-1)*(r1[0])/r**3)
    v1[1] =v0[1]+h*((-1)*(r1[1])/r**3)
    return [r1,v1]

def ME(N,FEPO,r0,v0):#FEPO = force evaluations per orbit
    h = (2*math.pi)/FEPO
    allx0 = []
    ally0 = []
    r = math.sqrt((r0[0]**2)+(r0[1]**2))
    E0 = (1/2)*(v0[0]**2+v0[1]**2)-(1/r)
    allE = [E0]
    allEerror = [0]
    L0 = r0[0]*v0[1]-r0[1]*v0[0]
    allL = [L0]
    allLerror = [0]
    for i in range(N*FEPO):
        allx0.append(r0[0])
        ally0.append(r0[1])
        r = math.sqrt((r0[0]**2)+(r0[1]**2))
        E1 = (1/2)*(v0[0]**2+v0[1]**2)-(1/r)
        allE.append(E1)
        allEerror.append(abs(E1-E0)/abs(E0))
        L1 = r0[0]*v0[1]-r0[1]*v0[0]
        allL.append(L1)
        allLerror.append(abs(L1-L0)/abs(L0))
        q = MEstep(r0,v0,h)
        r0 = q[0]
        v0 = q[1]
        
    t = []
    for i in range(N*FEPO+1):
        t.append(i)  
    x = []
    y = []
    for i in range(801):
        E = e*math.pi*(i/200)
        x.append(math.cos(E)+e)
        y.append(b*math.sin(E))
    return([allEerror,allLerror,t,allx0,ally0,r0,v0])

def MEgraph(N,FEPO):
    values = ME(N,FEPO,[a*(1+e),0],[0,math.sqrt((a*(1-(e)))/(a*(1+(e))))])
    plt.figure(1)
    plt.subplot(211)
    plt.plot(values[3],values[4])
    plt.plot(x,y)
    plt.plot(0,0,'r*')
    plt.figure(2)
    plt.loglog(values[2],values[0])
    plt.figure(3)
    plt.loglog(values[2],values[1])
    
def Leapfrogstep(r0,v0,h):
    rprime = [0,0]#initial definition
    v1 = [0,0]
    r1 = [0,0]
    rprime[0] = r0[0] + (h/2)*v0[0]
    rprime[1] = r0[1] + (h/2)*v0[1]
    r = math.sqrt((rprime[0]**2)+(rprime[1]**2))
    v1[0]= v0[0]+ h*((-1)*(rprime[0])/r**3)
    v1[1]= v0[1]+ h*((-1)*(rprime[1])/r**3)
    r1[0] = rprime[0] + (h/2)*v1[0]
    r1[1] = rprime[1] + (h/2)*v1[1]
    return [r1,v1]   

def Leapfrog(N,FEPO,r0,v0):#FEPO = force evaluations per orbit
    h = (2*math.pi)/FEPO
    allx0 = []
    ally0 = []
    r = math.sqrt((r0[0]**2)+(r0[1]**2))
    E0 = (1/2)*(v0[0]**2+v0[1]**2)-(1/r)
    allE = [E0]
    allEerror = [0]
    L0 = r0[0]*v0[1]-r0[1]*v0[0]
    allL = [L0]
    allLerror = [0]
    for i in range(N*FEPO):
        allx0.append(r0[0])
        ally0.append(r0[1])
        r = math.sqrt((r0[0]**2)+(r0[1]**2))
        E1 = (1/2)*(v0[0]**2+v0[1]**2)-(1/r)
        allE.append(E1)
        allEerror.append(abs(E1-E0)/abs(E0))
        L1 = r0[0]*v0[1]-r0[1]*v0[0]
        allL.append(L1)
        allLerror.append(abs(L1-L0)/abs(L0))
        q = Leapfrogstep(r0,v0,h)
        r0 = q[0]
        v0 = q[1]
        
    t = []
    for i in range(N*FEPO+1):
        t.append(i)  
    
    x = []
    y = []
    for i in range(801):
        E = e*math.pi*(i/200)
        x.append(math.cos(E)+e)
        y.append(b*math.sin(E))
    return([allEerror,allLerror,t,allx0,ally0,r0,v0])
    
def Leapfroggraph(N,FEPO):
    values = Leapfrog(N,FEPO,[a*(1+e),0],[0,math.sqrt((a*(1-(e)))/(a*(1+(e))))])
    plt.figure(1)
    plt.subplot(211)
    plt.plot(values[3],values[4])
    plt.plot(x,y)
    plt.plot(0,0,'r*')
    plt.figure(2)
    plt.loglog(values[2],values[0])
    plt.figure(3)
    plt.loglog(values[2],values[1])

def RK4step(r0,v0,h):#doesn't work yet, pls help based on runge kutta 4th order method
    r = math.sqrt((r0[0]**2)+(r0[1]**2))
    def f(x):
        return np.array([x[2],x[3],(((-1)/r**3)*(x[0])),(((-1)/r**3)*(x[1]))])
    w0 = [r0[0],r0[1],v0[0],v0[1]]
    k1 = f(w0)
    k1 = [i * h for i in k1]
    k2 = [0,0,0,0]
    for i in range(len(k2)):
        k2[i] = w0[i]+0.5*k1[i]
    k2 = f(k2)
    k2 = [i * h for i in k2]
    k3 = [0,0,0,0]
    for i in range(len(k3)):
        k3[i] = w0[i]+0.5*k2[i]
    k3 = f(k3)
    k3 = [i * h for i in k3]
    k4 = [0,0,0,0]
    for i in range(len(k4)):
        k4[i] = w0[i]+k3[i]
    k4 = f(k4)
    k4 = [i * h for i in k4]
    w1 = [0,0,0,0]
    for i in range(len(w1)):
        w1[i] = w0[i] + (1/6)*(k1[i]+2*k2[i]+2*k3[i]+k4[i])
    return w1

def RK4(N,FEPO,r0,v0):#FEPO = force evaluations per orbit
#    r0 = [a*(1+e),0]
#    v0 = [0,math.sqrt((a*(1-(e)))/(a*(1+(e))))]
    h = (2*math.pi)/(FEPO/4)
    r = math.sqrt((r0[0]**2)+(r0[1]**2))
    allx0 = []
    ally0 = []
    E0 = (1/2)*(v0[0]**2+v0[1]**2)-(1/r)
    allE = [E0]
    allEerror = [0]
    L0 = r0[0]*v0[1]-r0[1]*v0[0]
    allL = [L0]
    allLerror = [0]
    for i in range(N*(FEPO//4)):
        allx0.append(r0[0])
        ally0.append(r0[1])
        r = math.sqrt((r0[0]**2)+(r0[1]**2))
        E1 = (1/2)*(v0[0]**2+v0[1]**2)-(1/r)
        allE.append(E1)
        allEerror.append(abs(E1-E0)/abs(E0))
        L1 = r0[0]*v0[1]-r0[1]*v0[0]
        allL.append(L1)
        allLerror.append(abs(L1-L0)/abs(L0))
        q = RK4step(r0,v0,h)
        r0 = [q[0],q[1]]
        v0 = [q[2],q[3]]
    t = []
    for i in range(N*(FEPO//4)+1):
        t.append(i)
    x = []
    y = []
    for i in range(801):
        E = e*math.pi*(i/200)
        x.append(math.cos(E)+e)
        y.append(b*math.sin(E))
    return([allEerror,allLerror,t,allx0,ally0,r0,v0])
    
def RK4graph(N,FEPO):
    values = RK4(N,FEPO,[a*(1+e),0],[0,math.sqrt((a*(1-(e)))/(a*(1+(e))))])
    plt.figure(1)
    plt.subplot(211)
    plt.plot(values[3],values[4])
    plt.plot(x,y)
    plt.plot(0,0,'r*')
    plt.figure(2)
    plt.loglog(values[2],values[0])
    plt.figure(3)
    plt.loglog(values[2],values[1])
    
    
def errors(N,FEPO,e):
    FEerrors = FE(N,FEPO,[a*(1+e),0],[0,math.sqrt((a*(1-(e)))/(a*(1+(e))))])
    MEerrors = ME(N,FEPO,[a*(1+e),0],[0,math.sqrt((a*(1-(e)))/(a*(1+(e))))])
    Leapfrogerrors = Leapfrog(N,FEPO,[a*(1+e),0],[0,math.sqrt((a*(1-(e)))/(a*(1+(e))))])
    RK4errors = RK4(N,FEPO,[a*(1+e),0],[0,math.sqrt((a*(1-(e)))/(a*(1+(e))))])
    plt.figure(1)
    plt.subplot(211)
    plt.loglog(FEerrors[2],FEerrors[0])
    plt.loglog(MEerrors[2],MEerrors[0])
    plt.loglog(Leapfrogerrors[2],Leapfrogerrors[0])
    plt.loglog(RK4errors[2],RK4errors[0])
    plt.figure(2)
    plt.subplot(212)
    plt.loglog(FEerrors[2],FEerrors[1])
    plt.loglog(MEerrors[2],MEerrors[1])
    plt.loglog(Leapfrogerrors[2],Leapfrogerrors[1])
    plt.loglog(RK4errors[2],RK4errors[1])
    
    
def reversable(functionname):#Leapfrog is the only reversable case
    forward = functionname(10,1000,[a*(1+e),0],[0,math.sqrt((a*(1-(e)))/(a*(1+(e))))])
    backward = functionname(10,1000,forward[5],[-forward[6][0],-forward[6][1]])
    plt.plot(forward[3]+backward[3],forward[4]+backward[4])
    finalpos = backward[5]
    finalvel = [-backward[6][0],-backward[6][1]]
    poserrorx = abs(a*(1+e)-finalpos[0])
    poserrory = abs(0-finalpos[1])
    velerrorx = abs(0-finalvel[0])
    velerrory = abs(math.sqrt((a*(1-(e)))/(a*(1+(e))))-finalvel[1])
    return [poserrorx,poserrory,velerrorx,velerrory]

def allreversable():
    FEreverrors = reversable(FE)
    print ('The postional error for forward Euler is ' +str(math.sqrt(FEreverrors[0]**2+FEreverrors[1]**2)))
    print ('The velocity error for forward Euler is ' +str(math.sqrt(FEreverrors[2]**2+FEreverrors[3]**2)))
    MEreverrors = reversable(ME)
    print ('The postional error for modified Euler is ' +str(math.sqrt(MEreverrors[0]**2+MEreverrors[1]**2)))
    print ('The velocity error for modified Euler is ' +str(math.sqrt(MEreverrors[2]**2+MEreverrors[3]**2)))
    Leapfrogreverrors = reversable(Leapfrog)
    print ('The postional error for Leapfrog is ' +str(math.sqrt(Leapfrogreverrors[0]**2+Leapfrogreverrors[1]**2)))
    print ('The velocity error for Leapfrog is ' +str(math.sqrt(Leapfrogreverrors[2]**2+Leapfrogreverrors[3]**2)))
    RK4reverrors = reversable(RK4)
    print ('The postional error for RK4 is ' +str(math.sqrt(RK4reverrors[0]**2+RK4reverrors[1]**2)))
    print ('The velocity error for RK4 is ' +str(math.sqrt(RK4reverrors[2]**2+RK4reverrors[3]**2)))




def TwoBodyLeapfrogStep(first,second,h):
    rprimeone = [0,0]#initial definitions
    rprimetwo = [0,0]
    first = [first[0],first[1],[0,0],[0,0],[0,0],[0,0],first[2]] #[r0],[v0],[r1],[v1],[relative pos],[relative vel],mass
    second = [second[0],second[1],[0,0],[0,0],[0,0],[0,0],second[2]]
    CentreOfMass = [(first[0][0]*first[6]+second[0][0]*second[6])/(first[6]+second[6]),(first[0][1]*first[6]+second[0][1]*second[6])/(first[6]+second[6]),]#x pos,y pos
    first[4][0] = first[0][0] - CentreOfMass[0]
    first[4][1] = first[0][1] - CentreOfMass[1]
    first[5][0] = first[1][0]
    first[5][1] = first[1][1]
    second[4][0] = second[0][0] - CentreOfMass[0]
    second[4][1] = second[0][1] - CentreOfMass[1]
    second[5][0] = second[1][0]
    second[5][1] = second[1][1]
    rprimeone[0] = first[4][0] + (h/2)*first[5][0]
    rprimeone[1] = first[4][1] + (h/2)*first[5][1]
    r = math.sqrt((rprimeone[0]**2)+(rprimeone[1]**2))
    first[3][0]= first[5][0]+ h*((-1)*(rprimeone[0])/r**3)
    first[3][1]= first[5][1]+ h*((-1)*(rprimeone[1])/r**3)
    first[2][0] = rprimeone[0] + (h/2)*first[3][0]+CentreOfMass[0]
    first[2][1] = rprimeone[1] + (h/2)*first[3][1]+CentreOfMass[1]
    rprimetwo[0] = second[4][0] + (h/2)*second[5][0]
    rprimetwo[1] = second[4][1] + (h/2)*second[5][1]
    r = math.sqrt((rprimetwo[0]**2)+(rprimetwo[1]**2))
    second[3][0]= second[5][0]+ h*((-1)*(rprimetwo[0])/r**3)
    second[3][1]= second[5][1]+ h*((-1)*(rprimetwo[1])/r**3)
    second[2][0] = rprimetwo[0] + (h/2)*second[3][0]+CentreOfMass[0]
    second[2][1] = rprimetwo[1] + (h/2)*second[3][1]+CentreOfMass[1]
    return [[first[2],first[3],first[6]],[second[2],second[3],second[6]]]

def TwoBodyLeapfrog(N,FEPO,first,second):
    h = (2*math.pi)/FEPO
    allfirst = [[],[]]
    allsecond = [[],[]]
    for i in range(int(N*FEPO)):
        allfirst[0].append(first[0][0])
        allfirst[1].append(first[0][1])
        allsecond[0].append(second[0][0])
        allsecond[1].append(second[0][1])
        q = TwoBodyLeapfrogStep(first,second,h)
        first = q[0]
        second = q[1]
    return [allfirst,allsecond,first,second]

def TwoBodyLeapfrogGraph(N,FEPO):
    #values = TwoBodyLeapfrog(N, FEPO, [[1,0],[0,1],1], [[-1,0],[0,-1],1])
    values = TwoBodyLeapfrog(N,FEPO,[[a*(1+e)+3,3],[0,math.sqrt((a*(1-(e)))/(a*(1+(e))))],0.5],[[3,-(a*(1+e))+3],[0,-(math.sqrt((a*(1-(e)))/(a*(1+(e)))))],0.5])
    plt.figure(1)
    plt.subplot(211)
    plt.plot(values[0][0],values[0][1])
    plt.plot(values[1][0],values[1][1])
    

def NBodyLeapfrogStep(n,objects,h):#n = number of objects, objects is the actual objects as a list and it will be like two body leapfrog with [first,second,third,etc]
    rprime = []#inital defs
    TotalMass = 0
    CentreOfMass = [0,0]#x pos, y pos
    for i in range(n):
        rprime.append([0,0])
        objects[i] = [objects[i][0],objects[i][1],[0,0],[0,0],[0,0],[0,0],objects[i][2]]
        TotalMass += objects[i][6]
        CentreOfMass[0] += objects[i][0][0]
        CentreOfMass[1] += objects[i][0][1]
    CentreOfMass[0] = CentreOfMass[0]/TotalMass
    CentreOfMass[1] = CentreOfMass[1]/TotalMass
    for i in range(n):
        objects[i][4][0] = objects[i][0][0] - CentreOfMass[0]
        objects[i][4][1] = objects[i][0][1] - CentreOfMass[1]
        objects[i][5][0] = objects[i][1][0]
        objects[i][5][1] = objects[i][1][1]
        rprime[i][0] = objects[i][4][0] + (h/2)*objects[i][5][0]
        rprime[i][1] = objects[i][4][1] + (h/2)*objects[i][5][1]
        r = math.sqrt((rprime[i][0]**2)+(rprime[i][1]**2))
        objects[i][3][0]= objects[i][5][0]+ h*((-1)*(rprime[i][0])/r**3)
        objects[i][3][1]= objects[i][5][1]+ h*((-1)*(rprime[i][1])/r**3)
        objects[i][2][0] = rprime[i][0] + (h/2)*objects[i][3][0]+CentreOfMass[0]
        objects[i][2][1] = rprime[i][1] + (h/2)*objects[i][3][1]+CentreOfMass[1]
        objects[i] = [objects[i][2],objects[i][3],objects[i][6]]
    #print(CentreOfMass)
    return objects

def NBodyLeapfrog(N,FEPO,n,objects):
    h = (2*math.pi)/FEPO
    allobjects = []
    for i in range(n):
        allobjects.append([[],[]])
    for j in range(int(N*FEPO)):
        for i in range(n):
            allobjects[i][0].append(objects[i][0][0])
            allobjects[i][1].append(objects[i][0][1])
        objects = NBodyLeapfrogStep(n,objects,h)
    return [allobjects,objects]


def NBodyLeapfrogGraph(N,FEPO,n):
    values = NBodyLeapfrog(N,FEPO,n,[[[a*(1+e),0],[0,math.sqrt((a*(1-(e)))/(a*(1+(e))))],1],[[0,-(a*(1+e))],[0,-(math.sqrt((a*(1-(e)))/(a*(1+(e)))))],1],[[a*(1+e),-a*(1+e)],[0,math.sqrt((a*(1-(e)))/(a*(1+(e))))],1],[[0,0],[0,-math.sqrt((a*(1-(e)))/(a*(1+(e))))],1]])
    #values = NBodyLeapfrog(N, FEPO, n, [[[0,1],[1,0],1],[[-2/3,-1/3],[1/3,2/3],1],[[2/3,-1/3],[-1/3,-2/3],1]])
    #values = NBodyLeapfrog(N,FEPO,n,[[[0,0],[0,0],1000000],[[50,0],[0,500],100],[[-40,0],[0,-500],1000]])
    plt.figure(1)
    plt.subplot(211)
    for i in range(n):
        plt.plot(values[0][i][0],values[0][i][1])
        

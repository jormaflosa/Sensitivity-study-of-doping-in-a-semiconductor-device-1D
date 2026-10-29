import numpy as np
import matplotlib.pyplot as plt


def sech(z):
    return 1/np.cosh(z)

def u(x, t):
    z =(2*x-t)/(2*np.sqrt(3))
    u=(-1/(np.sqrt(10)))*sech(z)
    return u


def Q(x, t):
    z=(2*x-t)/(2*np.sqrt(3))
    n=(1/(np.sqrt(5)))*sech(z)
    return n

def D(x):
    D=1e-9
    return D

def n_i(x):
    z =(2*x)/(2*np.sqrt(3))
    n_i=(1/(np.sqrt(5)))*sech(z)
    return n_i

def f(x,t):
    z =(2*x-t)/(2*np.sqrt(3))
    f=-sech(z)**2*(1/(15*np.sqrt(2))-(1/(5*np.sqrt(2)))*np.tanh(z)**2)+sech(z)*(D(x)/(3*np.sqrt(5))-(2*D(x)/(3*np.sqrt(5)))*np.tanh(z)**2-(1/(2*np.sqrt(15)))*np.tanh(z))
    return -f

def g(x,t):
    z = (2*x - t)/(2*np.sqrt(3))
    g=sech(z)*(1/np.sqrt(5)-1/(3*np.sqrt(10))+(np.sqrt(2)/(3*np.sqrt(5)))*np.tanh(z)**2)
    return g

def N(n):
    N=np.diag(np.full(n-1,2))-np.diag(np.ones(n-2),1)-np.diag(np.ones(n-2),-1)
    return N

def r1(t,h,n,u_0,u_1,m):
    x_i=h
    if t==0:
        j=0
        m=np.zeros((n-1,1))
        while j<n-1:
            m[j]=n_i(x_i)
            x_i=x_i+h
            j=j+1
    x_i=h
    r=np.zeros((n-1,1))
    i=0
    while i<n-1:
        if i==0:
            r[i]=(g(x_i,t)-m[i])*(h**2)+u_0
        elif i==n-2:
            r[i]=(g(x_i,t)-m[i])*(h**2)+u_1
        else:
            r[i]=(g(x_i,t)-m[i])*(h**2)
        i=i+1
        x_i=x_i+h
    return r

def qmas(h,m,i):
    q=(m[i+1]-m[i])/h
    return float(q.item())

def qmen(h,m,i):
    q=(m[i]-m[i-1])/h
    return float(q.item())

def B(l,n,s,h):  
    x_i=h  
    i=1
    B=np.empty((0, n-1))
    while i<=n-1:
        a=(qmas(h,s,i)/2)-((abs(qmas(h,s,i)))/2)
        b=(qmas(h,s,i)/2)+((abs(qmas(h,s,i)))/2)
        c=(qmen(h,s,i)/2)-((abs(qmen(h,s,i)))/2)
        d=(qmen(h,s,i)/2)+((abs(qmen(h,s,i)))/2)
        A=[]
        j=1
        while j<=n-1:
            if j==i:
                A.append(((h**2)/l)+(h*b-h*c)+D(x_i+h/2)+D(x_i-h/2))
            elif j==i+1:
                A.append(h*a-D(x_i+h/2))
            elif j>i+1:
                A.append(0)
            elif j==i-1:
                A.append(-h*d-D(x_i-h/2))
            elif j<i-1:
                A.append(0)
            j=j+1
        i=i+1
        x_i=x_i+h
        B = np.vstack([B, A])
    return(B)

def r2(t,h,l,n,m,u_0,u_1,s):
    x_i=h
    if t==tiempos[1]:
        j=0
        m=np.zeros((n-1,1))
        while j<n-1:
            m[j]=n_i(x_i)
            x_i=x_i+h
            j=j+1
    x_i=h
    r=np.zeros((n-1,1))
    i=0
    while i<n-1:
        d=(qmen(h,s,i+1)/2)+((abs(qmen(h,s,i+1)))/2)
        a=(qmas(h,s,i+1)/2)-((abs(qmas(h,s,i+1)))/2)
        if i==0:
            r[i]=(h**2)*(f(x_i,t)+m[i]*((1)/l))+(u_0*((h*d)+D(1)))
        elif i==n-2:
            r[i]=(h**2)*(f(x_i,t)+m[i]*((1)/l))+u_1*(D(1)-h*a)
        else:
            r[i]=(h**2)*(f(x_i,t)+m[i]*((1)/l))
        i=i+1
        x_i=x_i+h
    return r



def DD(p,n,u_0,u_1):
    A=np.zeros((n+1,1))
    A[0,0]=u_0
    A[n,0]=u_1
    i=1
    while i<=n-1:
        A[i,0]=p[i-1,0]
        i=i+1
    return(A)


z=[25,50,75,100,125,150,175,200,225,250,275,300,325,350]
a=0
E=[]
while a<14:
    n=z[a]
    h=(2*np.pi)/n
    x_nodos=np.linspace(0, 2*np.pi, n+1)
    x = np.linspace(0, 2*np.pi, 500)
    Tf=1
    l=(h)/(2)
    k=int(Tf/l)
    j=1
    tiempos = [j*l for j in range(k+1)]
    p=np.linalg.solve(N(n),r1(0,h,n,(-1/np.sqrt(10))*(1/np.cosh((0)/(2*np.sqrt(3)))),(-1/np.sqrt(10))*(1/np.cosh((4*np.pi)/(2*np.sqrt(3)))),0))
    s=DD(p,n,(-1/np.sqrt(10))*(1/np.cosh((0)/(2*np.sqrt(3)))),(-1/np.sqrt(10))*(1/np.cosh((4*np.pi)/(2*np.sqrt(3)))))
    e=np.linalg.solve(B(l,n,s,h),r2(tiempos[1],h,l,n,0,(1/np.sqrt(5))*(1/np.cosh((-tiempos[1])/(2*np.sqrt(3)))),(1/np.sqrt(5))*(1/np.cosh((4*np.pi-tiempos[1])/(2*np.sqrt(3)))),s))
    w=DD(e,n,(1/np.sqrt(5))*(1/np.cosh((-tiempos[1])/(2*np.sqrt(3)))),(1/np.sqrt(5))*(1/np.cosh((4*np.pi-tiempos[1])/(2*np.sqrt(3)))))

    j=1
    while j<np.size(tiempos)-1:
        p=np.linalg.solve(N(n),r1(tiempos[j],h,n,(-1/np.sqrt(10))*(1/np.cosh((-tiempos[j])/(2*np.sqrt(3)))),(-1/np.sqrt(10))*(1/np.cosh((4*np.pi-tiempos[j])/(2*np.sqrt(3)))),e))
        s=DD(p,n,(-1/np.sqrt(10))*(1/np.cosh((-tiempos[j])/(2*np.sqrt(3)))),(-1/np.sqrt(10))*(1/np.cosh((4*np.pi-tiempos[j])/(2*np.sqrt(3)))))
        e=np.linalg.solve(B(l,n,s,h),r2(tiempos[j+1],h,l,n,e,(1/np.sqrt(5))*(1/np.cosh((-tiempos[j+1])/(2*np.sqrt(3)))),(1/np.sqrt(5))*(1/np.cosh((4*np.pi-tiempos[j+1])/(2*np.sqrt(3)))),s))
        w=DD(e,n,(1/np.sqrt(5))*(1/np.cosh((-tiempos[j+1])/(2*np.sqrt(3)))),(1/np.sqrt(5))*(1/np.cosh((4*np.pi-tiempos[j+1])/(2*np.sqrt(3)))))
        j=j+1
    p=np.linalg.solve(N(n),r1(tiempos[j],h,n,(-1/np.sqrt(10))*(1/np.cosh((-tiempos[j])/(2*np.sqrt(3)))),(-1/np.sqrt(10))*(1/np.cosh((4*np.pi-tiempos[j])/(2*np.sqrt(3)))),e))
    s=DD(p,n,(-1/np.sqrt(10))*(1/np.cosh((-tiempos[j])/(2*np.sqrt(3)))),(-1/np.sqrt(10))*(1/np.cosh((4*np.pi-tiempos[j])/(2*np.sqrt(3)))))
    i=0
    sum=0
    dem=0
    while i<=n:
        U=(s[i]-u(x_nodos[i],tiempos[k]))**2
        sum=sum+U
        R=u(x_nodos[i],tiempos[k])**2
        dem=dem+R
        i=i+1  
    E.append(float(np.sqrt(sum/dem))) 
    print(E)
    a=a+1

plt.plot(z, E)

plt.xlabel("n")
plt.ylabel("Error relativo")
plt.title("Error relativo medio")
plt.legend()
plt.grid(True)
plt.show()
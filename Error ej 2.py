import numpy as np
import matplotlib.pyplot as plt


def u(x, t):
    u=np.sin(t)*x*((x-1)**2)
    return u

def Q(x, t):
    n=np.cos(t)*x*(1-x)
    return n


def D(x):
    if x<0.5:
        D=0.025
    elif x==0.5:
        D=(2*0.025*0.1)/(0.025+0.5)
    else:
        D=0.1
    return D

def wp(x_i,h):
    w=(2*D(x_i)*D(x_i+h))/(D(x_i)+D(x_i+h))
    return w

def wn(x_i,h):
    w=(2*D(x_i-h)*D(x_i))/(D(x_i-h)+D(x_i))
    return w

def n_i(x):
    n_i=x*(1-x)
    return n_i

def f(x,t):
    f=np.sin(t)*(x*(x-1))+np.cos(t)*np.sin(t)*(1-x)*(12*(x**2)-9*x+1)+2*D(x)*np.cos(t)
    return f

def g(x,t):
    g=-np.sin(t)*(6*x-4)+np.cos(t)*(x*(1-x))
    return g

def N(n):
    N=np.diag(np.full(n-1,2))-np.diag(np.ones(n-2),1)-np.diag(np.ones(n-2),-1)
    return N

def r1(t,h,n,u_0,u_1,m):
    x_i= np.linspace(h, 1-h, n - 1) 
    if t==0:
        m=np.zeros((n-1,1))
        for j in range(n-1):
            m[j]=n_i(x_i[j])
    r=np.zeros((n-1,1))

    for i in range(n-1):
        if i==0:
            r[i]=(g(x_i[i],t)-m[i])*(h**2)+u_0
        elif i==n-2:
            r[i]=(g(x_i[i],t)-m[i])*(h**2)+u_1
        else:
            r[i]=(g(x_i[i],t)-m[i])*(h**2)
    return r

def qmas(h,m,i):
    q=(m[i+1]-m[i])/h
    return float(q.item())

def qmen(h,m,i):
    q=(m[i]-m[i-1])/h
    return float(q.item())

def B(l,n,s,h):  
    x_i=np.linspace(0, 1, n+1) 
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
                A.append(((h**2)/l)+(h*b-h*c)+wp(x_i[i],h)+wn(x_i[i],h))
            elif j==i+1:
                A.append(h*a-wp(x_i[i],h))
            elif j>i+1:
                A.append(0)
            elif j==i-1:
                A.append(-h*d-wn(x_i[i],h))
            elif j<i-1:
                A.append(0)
            j=j+1
        i=i+1
        B = np.vstack([B, A])
    return(B)

def r2(t,h,l,n,m):
    x_i= np.linspace(h, 1-h, n - 1)
    x2=h 
    if t==tiempos[1]:
        j=0
        m=np.zeros((n-1,1))
        for j in range(n-1):
            m[j]=n_i((x_i[j]+x2)/2)
            x2=x2+h
    x2=h
    r=np.zeros((n-1,1))
    for i in range (n-1):
        if i==0:
            r[i]=(h**2)*(f((x_i[i]+x2)/2,t)+m[i]*((1)/l))
        elif i==n-2:
            r[i]=(h**2)*(f((x_i[i]+x2)/2,t)+m[i]*((1)/l))
        else:
            r[i]=(h**2)*(f((x_i[i]+x2)/2,t)+m[i]*((1)/l))
        x2=x2+h
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

E=[]
for a in z:
    print(a)
    n=a
    h=1/n
    x_nodos=np.linspace(0, 1, n+1)
    x = np.linspace(0, 1, 500)

    Tf=1
    l=(h**2)/(2*0.1)
    k=int(Tf/l)
    j=1
    tiempos = [j*l for j in range(k+1)]
    p=np.linalg.solve(N(n),r1(0,h,n,0,0,0))
    s=DD(p,n,0,0)
    e=np.linalg.solve(B(l,n,s,h),r2(tiempos[1],h,l,n,0))
    w=DD(e,n,0,0)

    j=1
    while j<np.size(tiempos)-1:
        p=np.linalg.solve(N(n),r1(tiempos[j],h,n,0,0,e))
        s=DD(p,n,0,0)
        e=np.linalg.solve(B(l,n,s,h),r2(tiempos[j+1],h,l,n,e))
        w=DD(e,n,0,0)
        j=j+1
    p=np.linalg.solve(N(n),r1(tiempos[j],h,n,0,0,e))
    s=DD(p,n,0,0)
    i=0
    sum=0
    dem=0
    while i<=n:
        U=(s[i]-u(x_nodos[i],tiempos[k]))**2
        sum=sum+U
        R=u(x_nodos[i],tiempos[k])**2
        dem=dem+R
        i=i+1
    E.append(float(np.sqrt(sum.item()/dem.item())))



plt.plot(z, E)

plt.xlabel("phi")
plt.ylabel("Error relativo")
plt.title("Error relativo medio")
plt.legend()
plt.grid(True)
plt.show()

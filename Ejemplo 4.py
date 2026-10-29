import numpy as np
import matplotlib.pyplot as plt

def nd(x):
    if 0 <= x < 0.1:
        n = 5e5
    elif 0.5 < x <= 0.6:
        n = 5e5
    elif 0.1 <= x <= 0.15:
        p = (x - 0.1) / 0.05
        n = 5e5+4.98e5*(2*(p**3)-3*(p**2))
    elif 0.45 <= x <= 0.5:
        p = (x - 0.45) / 0.05
        n = 2e3-4.98e5*(2*(p**3)-3*(p**2))
    elif 0.15 < x < 0.45:
        n = 2e3
    return n

def mu(x):
    mu=0.0088*(1+((14.2273)/(1+((nd(x))/143200))))
    return mu

def D(x):
    k=0.138*(1e-4)
    T=300
    e=0.1602
    D=(mu(x)*k*T)/e
    return D

def n_i(x):
    n_i=nd(x)
    return n_i

def iphi(x):
    k=0.138*(1e-4)
    T=300
    e=0.1602
    
    phi=((k*T)/e)*np.log(5e17/1.4e10)+0.55
    return phi

def dphi(x):
    phi=iphi(x)+(1.5)
    return phi


def f(x,t):
    f=0
    return f

def g(x,t):
    g=nd(x)
    return g

def N(n):
    N=((11.7*8.85418)/0.1602)*(np.diag(np.full(n-1,2))-np.diag(np.ones(n-2),1)-np.diag(np.ones(n-2),-1))
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
            r[i]=(g(x_i,t)-m[i])*(h**2)+((11.7*8.85418)/0.1602)*u_0
        elif i==n-2:
            r[i]=(g(x_i,t)-m[i])*(h**2)+((11.7*8.85418)/0.1602)*u_1
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
        a=mu(x_i+h/2)*((qmas(h,s,i)/2)-((abs(qmas(h,s,i)))/2))
        b=mu(x_i+h/2)*((qmas(h,s,i)/2)+((abs(qmas(h,s,i)))/2))
        c=mu(x_i-h/2)*((qmen(h,s,i)/2)-((abs(qmen(h,s,i)))/2))
        d=mu(x_i-h/2)*((qmen(h,s,i)/2)+((abs(qmen(h,s,i)))/2))
        A=[]
        j=1
        while j<=n-1:
            if j==i:
                A.append(((h**2)/l)+(h*b-h*c)+D(x_i)+D(x_i))
            elif j==i+1:
                A.append(h*a-D(x_i))
            elif j>i+1:
                A.append(0)
            elif j==i-1:
                A.append(-h*d-D(x_i))
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
        d=mu(x_i-h/2)*((qmen(h,s,i+1)/2)+((abs(qmen(h,s,i+1)))/2))
        a=mu(x_i+h/2)*((qmas(h,s,i+1)/2)-((abs(qmas(h,s,i+1)))/2))
        if i==0:
            r[i]=(h**2)*(f(x_i,t)+m[i]*((1)/l))+(u_0*((h*d)+D(x_i)))
        elif i==n-2:
            r[i]=(h**2)*(f(x_i,t)+m[i]*((1)/l))+u_1*(D(x_i)-h*a)
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

#Parametros
n=int(input("Número de divisiones: "))
print(f"El número de nodos es: {n+1}")
x_nodos=np.linspace(0, 0.6, n+1)
h=0.6/n
Tf=1
l=1/n
k=int(Tf/l)
j=1
tiempos = [j*l for j in range(k+1)]
fig, axs = plt.subplots(1, 3, figsize=(18, 5))


#Solver
j=0
while j<np.size(tiempos)-1:
    if j==0:
        p=np.linalg.solve(N(n),r1(0,h,n,iphi(1),dphi(1),0))
        s=DD(p,n,iphi(1),dphi(1))
        e=np.linalg.solve(B(l,n,s,h),r2(tiempos[1],h,l,n,0,5e5,5e5,s))
        w=DD(e,n,5e5,5e5)
    else:
        p=np.linalg.solve(N(n),r1(tiempos[j],h,n,iphi(1),dphi(1),e))
        s=DD(p,n,iphi(1),dphi(1))
        e=np.linalg.solve(B(l,n,s,h),r2(tiempos[j+1],h,l,n,e,5e5,5e5,s))
        w=DD(e,n,5e5,5e5)
    j=j+1
p=np.linalg.solve(N(n),r1(tiempos[j],h,n,iphi(1),dphi(1),e))
s=DD(p,n,iphi(1),dphi(1))

i=0
E=np.zeros((n+1,1))
while i<np.size(s)-1:
    E[i,0]=(-(s[i+1]-s[i])/h).item()
    i=i+1
E[i,0]=(-(s[i]-s[i-1])/h).item()

#Gráficas
axs[0].plot(x_nodos, w, color='red', label=f"Solución númerica de etha en t={tiempos[j-1]}")
axs[1].plot(x_nodos, s, color='red', label=f"Solución númerica de phi en t={tiempos[j-1]}")
axs[2].plot(x_nodos, E, color='red', label=f"Solución númerica de E en t={tiempos[j-1]}")

axs[2].set_title('Solución E')
axs[2].set_xlabel('x')
axs[2].set_ylabel('E')
axs[2].grid(True)

axs[1].set_title('Solución $\Phi$')
axs[1].set_xlabel('x')
axs[1].set_ylabel('$\phi$')
axs[1].grid(True)

axs[0].set_title('Solución $\eta$')
axs[0].set_xlabel('x')
axs[0].set_ylabel('$eta$')
axs[0].grid(True)


plt.tight_layout()
plt.show()
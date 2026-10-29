import numpy as np
import matplotlib.pyplot as plt
import random

def nd(x,delta):
    n=1000*(2*x**2-2*x+1)
    return (1+delta)*n


def mu(x,delta):
    mu=1.45e3
    return mu

def D(x,delta):
    k=1.380649*1e-23
    T=300
    e=1.6022e-19
    D=(mu(x,delta)*k*T)/e
    return D

def n_i(x,delta):
    n_i=nd(x,0)
    return n_i

def iphi(x):
    k=1.380649*1e-23
    T=300
    e=1.6022e-19
    
    phi=((k*T)/e)*np.log(1e16/9.65e9)
    return phi

def dphi(t):
    phi=iphi(1)-0.5*np.sin((np.pi*t)/1e5)
    return phi


def f(x,t):
    f=0
    return f

def g(x,t,delta):
    g=nd(x,delta)
    return g

def N(n):
    N=(1.0536e-12/1.6022e-19)*(np.diag(np.full(n-1,2))-np.diag(np.ones(n-2),1)-np.diag(np.ones(n-2),-1))
    return N

def r1(t,h,n,u_0,u_1,m,delta):
    x_i=h
    if t==0:
        j=0
        m=np.zeros((n-1,1))
        while j<n-1:
            m[j]=n_i(x_i,delta)
            x_i=x_i+h
            j=j+1
    x_i=h
    r=np.zeros((n-1,1))
    i=0
    while i<n-1:
        if i==0:
            r[i]=(g(x_i,t,delta)-m[i])*(h**2)+(1.0536e-12/1.6022e-19)*u_0
        elif i==n-2:
            r[i]=(g(x_i,t,delta)-m[i])*(h**2)+(1.0536e-12/1.6022e-19)*u_1
        else:
            r[i]=(g(x_i,t,delta)-m[i])*(h**2)
        i=i+1
        x_i=x_i+h
    return r

def qmas(h,m,i):
    q=(m[i+1]-m[i])/h
    return float(q.item())

def qmen(h,m,i):
    q=(m[i]-m[i-1])/h
    return float(q.item())

def B(l,n,s,h,delta):  
    x_i=h  
    i=1
    B=np.empty((0, n-1))
    while i<=n-1:
        a=mu(x_i+h/2,delta)*((qmas(h,s,i)/2)-((abs(qmas(h,s,i)))/2))
        b=mu(x_i+h/2,delta)*((qmas(h,s,i)/2)+((abs(qmas(h,s,i)))/2))
        c=mu(x_i-h/2,delta)*((qmen(h,s,i)/2)-((abs(qmen(h,s,i)))/2))
        d=mu(x_i-h/2,delta)*((qmen(h,s,i)/2)+((abs(qmen(h,s,i)))/2))
        A=[]
        j=1
        while j<=n-1:
            if j==i:
                A.append(((h**2)/l)+(h*b-h*c)+D(x_i,delta)+D(x_i,delta))
            elif j==i+1:
                A.append(h*a-D(x_i,delta))
            elif j>i+1:
                A.append(0)
            elif j==i-1:
                A.append(-h*d-D(x_i,delta))
            elif j<i-1:
                A.append(0)
            j=j+1
        i=i+1
        x_i=x_i+h
        B = np.vstack([B, A])
    return(B)

def r2(t,h,l,n,m,u_0,u_1,s,delta):
    x_i=h
    if t==tiempos[1]:
        j=0
        m=np.zeros((n-1,1))
        while j<n-1:
            m[j]=n_i(x_i,delta)
            x_i=x_i+h
            j=j+1
    x_i=h
    r=np.zeros((n-1,1))
    i=0
    while i<n-1:
        d=mu(x_i-h/2,delta)*((qmen(h,s,i+1)/2)+((abs(qmen(h,s,i+1)))/2))
        a=mu(x_i+h/2,delta)*((qmas(h,s,i+1)/2)-((abs(qmas(h,s,i+1)))/2))
        if i==0:
            r[i]=(h**2)*(f(x_i,t)+m[i]*((1)/l))+(u_0*((h*d)+D(x_i,delta)))
        elif i==n-2:
            r[i]=(h**2)*(f(x_i,t)+m[i]*((1)/l))+u_1*(D(x_i,delta)-h*a)
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
#delta=[-0.1,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1]

delta=[]
for a in range(41):
    delta.append(0.1*random.normalvariate(0, 1))
n=int(input("Número de divisiones: "))
print(f"El número de nodos es: {n+1}")
x_nodos=np.linspace(0, 1, n+1)
h=1/n
Tf=0.1
l=1/n
k=int(Tf/l)
j=1
tiempos = [j*l for j in range(k+1)]
print(tiempos)
historial_soluciones1 = []
historial_soluciones2 = []
fig, axs = plt.subplots(1, 3, figsize=(18, 5))
delta.append(0)
delta=sorted(delta)
#Solución 
for v in delta:
    #Solver
    j=0
    while j<np.size(tiempos)-1:
        if j==0:
            p=np.linalg.solve(N(n),r1(0,h,n,iphi(1),dphi(tiempos[j]),0,v))
            s=DD(p,n,iphi(1),dphi(tiempos[j]))
            e=np.linalg.solve(B(l,n,s,h,v),r2(tiempos[1],h,l,n,0,1000,1000,s,v))
            w=DD(e,n,1000,1000)
        else:
            p=np.linalg.solve(N(n),r1(tiempos[j],h,n,iphi(1),dphi(tiempos[j]),e,v))
            s=DD(p,n,iphi(1),dphi(tiempos[j]))
            e=np.linalg.solve(B(l,n,s,h,v),r2(tiempos[j+1],h,l,n,e,1000,1000,s,v))
            w=DD(e,n,1000,1000)
        j=j+1
    p=np.linalg.solve(N(n),r1(tiempos[j],h,n,iphi(1),dphi(tiempos[j]),e,v))
    s=DD(p,n,iphi(1),dphi(tiempos[j]))
    i=0
    E=np.zeros((n+1,1))
    E[i,0]=(-(-3*s[i]+4*s[i+1]-s[i+2])/(2*h)).item()
    i=1
    while i<np.size(s)-1:
        E[i,0]=(-(s[i+1]-s[i])/h).item()
        i=i+1
    E[i,0]=(-(3*s[i]-4*s[i-1]+s[i-2])/(2*h)).item()

    if v==0:
        solucion_original1=w
        solucion_original2=s
        campoe=E
    else:
        axs[0].plot(x_nodos, w, color='blue', label=f"Solución númerica de etha en delta={v}")
        axs[1].plot(x_nodos, s, color='blue', label=f"Solución númerica de phi en delta={v}")
        axs[2].plot(x_nodos, E, color='blue', label=f"Solución númerica de E en delta={v}")
    historial_soluciones1.append(w.copy())
    historial_soluciones2.append(s.copy())

axs[0].plot(x_nodos, solucion_original1, color='red', label=f"Solución númerica de etha en t={tiempos[j-1]}")
axs[1].plot(x_nodos, solucion_original2, color='red', label=f"Solución númerica de phi en t={tiempos[j-1]}")
axs[2].plot(x_nodos, campoe, color='red', label=f"Solución númerica de E en t={tiempos[j-1]}")



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


historial_errores = []
for v in historial_soluciones1:
    numerador = np.linalg.norm(v-solucion_original1)
    denominador = np.linalg.norm(solucion_original1)
    error_actual_eta= numerador / denominador
    historial_errores.append(error_actual_eta)
plt.figure()
plt.plot(delta, historial_errores, '-o', color='blue', label='Error relativo de eta')
plt.xlabel('Magnitud de la perturbación ($\delta$)')
plt.ylabel('Error relativo $\eta$')
plt.title('Sensibilidad del sistema ante perturbaciones')
plt.grid(True)  

historial_errores2=[]
for v in historial_soluciones2:
    numerador = np.linalg.norm(v-solucion_original2)
    denominador = np.linalg.norm(solucion_original2)
    error_actual_eta= numerador / denominador
    historial_errores2.append(error_actual_eta)
plt.figure()
plt.plot(delta, historial_errores2, '-o', color='blue', label='Error relativo de phi')
plt.xlabel('Magnitud de la perturbación ($\delta$)')
plt.ylabel('Error relativo $\phi$')
plt.title('Sensibilidad del sistema ante perturbaciones')
plt.grid(True)  
    
    
    
    
plt.tight_layout()
plt.show()
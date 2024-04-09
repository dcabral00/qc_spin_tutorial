import tt
import tt.ksl

import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.special import jv # Bessel function of the second kind


dpi=300
plt.style.use('plot_style.txt')
plt.rcParams['axes.linewidth'] = 1.5
plt.rcParams['lines.markersize'] = 6 


# Heisenberg Hamiltonian in TT-format
def ham(N):
    # defining the Pauli matrices (in numpy array format)
    x = np.array([[0,1],[1,0]])
    y = np.array([[0,complex(0,-1)],[complex(0,1),0]])
    z = np.array([[1,0],[0,-1]])

    # converting the Pauli matrices to tt.matrix format
    tt_x=tt.matrix(x)
    tt_y=tt.matrix(y)
    tt_z=tt.matrix(z)

    # rounding threshold
    eps  = 1e-14

    # parameters of the model
    eps0 = 1.3/2 # on-site potential coefficent term (for 1st site)
    epsn = 2.0/2 # on-site potential coefficent term
    v0   = 0.75/2 # off-site kinetic coefficent term (for 1st site)
    vn   = 1.0/2 # off-site kinetic coefficent term
    # on-site potential hamiltonian matrix, initialized to zeros
    tt_H=0*tt.eye(2,N)
    tmp=tt_H
    k=0
    while k < N:
        if k == 0:
            tmp=tt.matrix(eps0*z)
            tmp=tt.kron(tmp,tt.eye(2,N-1))
        else:
            tmp=tt.kron(tt.eye(2,k),tt.matrix(epsn*z))
            tmp=tt.kron(tmp,tt.eye(2,N-k-1))
        tt_H=tt_H+tmp
        tt_H=tt_H.round(eps)
        k=k+1
    # off-site potential matrices, initialized to zero
    tt_HNNx=0*tt.eye(2,N) # XX
    tt_HNNy=tt_HNNx # YY
    tt_HNNz=tt_HNNx # ZZ
    k=0
    while k < N-1:
        if k == 0:
            tmpX=tt.kron(tt_x,tt_x)*v0
            tmpY=tt.kron(tt_y,tt_y)*v0
            tmpZ=tt.kron(tt_z,tt_z)*0
        else:
            tmp=tt.eye(2,k);
            tmpX=tt.kron(tmp,tt.kron(tt_x,tt_x)*vn)
            tmpY=tt.kron(tmp,tt.kron(tt_y,tt_y)*vn)
            tmpZ=tt.kron(tmp,tt.kron(tt_z,tt_z)*0)
        if k < N-2:
            tmpX=tt.kron(tmpX,tt.eye(2,N-2-k))
            tmpY=tt.kron(tmpY,tt.eye(2,N-2-k))
            tmpZ=tt.kron(tmpZ,tt.eye(2,N-2-k))
        tt_HNNx=tt_HNNx+tmpX
        tt_HNNy=tt_HNNy+tmpY
        tt_HNNz=tt_HNNz+tmpZ
        k=k+1
    # building the full Heisenberg Hamiltonian matrix and rounding
    tt_H=tt_H+tt_HNNx+tt_HNNy+tt_HNNz
    tt_H=tt_H.round(eps)
    return tt_H


# Taylor's expansion of exp(A): scaling and squaring algorithm
# A TT-matrix
def expA(A,eps):
    global Nspins,rmax
    N=10 # number of terms in Taylor's expansion
    w0=A*(1.0/2**N)
    e=tt.eye(2,Nspins)
    e=tt.matrix(e)
    tm=e
    k=N-1
    while k > 0:
        tm=e+tm*w0*(1.0/k)
        tm=tm.round(eps,rmax)
        k=k-1
    while k < N:
        tm=tm*tm
        tm=tm.round(eps,rmax)
        k=k+1
    return tm


# Chebyshev expansion of exp(A)*v, Clenshaw algorithm
# A TT-matrix, v TT-vector
def exp_clen(A,v):
    eps=1e-14
    a=-1;b=1;n=12
    #  Dp = a+b
    Dm = b-a
    d0=0*v; d1=d0
    for j in range(n):
        d2=d1
        d1=d0
        k = n-1-j
        ck=(-1j)**(k)*jv(k,1j)
    #    d0=ck*v+tt.matvec(A,d1)*(4/Dm)-2*d1*(Dp/Dm)-d2
        d0=ck*v+tt.matvec(A,d1)*(4/Dm)-d2
        d0=d0.round(eps)
    return d0-d2


# Clenshaw algorithm: Generates a tt of an arbitrary function func
# x is in TT format, either TT-vector or TT-matrix
# v could be provided as an argument to compute func*v (as in exp_clen)
def func_clen(x,func,v=None):
    b=1;a=-1; bma = 0.5 * (b - a); bpa = 0.5 * (b + a); n=12; eps=1e-14
    f = [func(math.cos(math.pi * (k + 0.5) / n) * bma + bpa) for k in range(n)]
    fac = 2.0 / n
    #c = [fac * sum([f[k] * math.cos(math.pi * j * (k + 0.5) / n)
    #          for k in range(n)]) for j in range(n)]
    c = [2*(1j)**j * jv(j,-1j) for j in range(n)] # when func = np.exp
    nf=0
    if v == None:
        v=tt.eye(x.n[0],x.d)               # for exp(x), when x is a tt-matrix
        #  v=tt.ones(x.n[0],x.d)              # for exp(x), when x is a tt-vector
    else:
        nf=1

    y=x
    y2 = 2.0 * y
    (d, dd) = (c[-1]*v, 0*v)             # Special case first step for efficiency
    # Clenshaw's recurrence
    for cj in c[-2:0:-1]:
        if nf == 1:
            (d, dd) = (tt.matvec(y2,d) - dd + v * cj, d)    # to compute func*v
        else:
            (d, dd) = (y2 * d - dd + v * cj, d)
        d=d.round(eps)
        dd=dd.round(eps)
    if nf == 1:
        out = tt.matvec(y, d) - dd + 0.5 * c[0]*v        # to compute func*v
    else:
        out = y * d - dd + 0.5 * c[0]*v
    return out.round(eps)


if __name__ == '__main__':
    # Simulation
    N = 3 
    Nspins = N # number of spins in the system
    rmax=4
    su = np.array([1,0]) # spin-up vector
    sd = np.array([0,1]) # spin down vector
    tau=0.1 # time step
    tfin=25 # final time; this does not appear to be used
    nsteps=250 # number of time steps (resulting in total time tau * nsteps = 25)
    eps=1e-14 # rounding threshold
    
    y0=tt.tensor(su, eps) # inital state, spin up on first site
    for j in range(N-1):
        y0=tt.kron(y0,tt.tensor(sd, eps)) # and spin-down on the subsequent sites

    A = ham(N)*1j #complex(0.0,1.0), ie iH
    Atau=A*tau # iHt
    eA= expA(Atau,eps) # Taylor expansion of exp(Atau), exp(iHT)
    eAc=func_clen(Atau,np.exp) # Chebyshev expansion of exp(Atau)
    
    t= np.arange(0,tau*(nsteps+1),tau) # [0,250*0.1, 0.1]
    p = np.empty_like(t)
    pr = np.empty_like(t)
    pc = np.empty_like(t)
    
    radd = 2 #Modify this to increase the rank
    if ( radd > 0 ):
        y0 = y0 + 1e-12*tt.rand(y0.n, y0.d, radd) #Hack, better by initialization
    
    #y0 = y0.round(-1)
    
    y0 = y0.round(1e-14) # initial state
    y = y0.copy()
    yr = y0.copy()
    yc = y0.copy()
    
    for k in range(nsteps+1):
        if k > 0:
            y = tt.ksl.ksl(A, y, tau);
            yr = tt.matvec(eA, yr)
            yr=yr.round(eps,radd)
    #        yc = tt.matvec(eAc, yc) # commented out Chebyshev
            yc=exp_clen(Atau,yc)
    ##        yc=exp_clen(Atau*k,y0) # commented out Chebyshev
            yc=yc.round(eps,radd)
            ###        yc=func_clen(Atau,np.exp,yc)
        p[k] = abs(tt.dot(y,y0))
        pr[k] = abs(tt.dot(yr,y0))
        pc[k] = abs(tt.dot(yc,y0))
    
    np.save(f'TT_{N}_spin_chain_time', t)
    np.save(f'TT_{N}_spin_chain_sa_ksl', p)
    np.save(f'TT_{N}_spin_chain_sa_Taylor', pr)
    np.save(f'TT_{N}_spin_chain_sa_Chebyshev', pc)

    plt.yscale('log')
    plt.plot(t,p,'b',label='ksl')
    plt.plot(t,pr,'-',label='Taylor')
    plt.plot(t,pc,'-',label='Chebyshev')
    plt.title('Heisenberg Model')
    plt.xlabel('Time')
    plt.ylabel('Absolute Value of Survival Amplitude')
    # plt.xscale('log')
    plt.legend()
    # plt.legend(loc='lower center', bbox_to_anchor=(0.5, -0.30),
    #           ncol=3, fancybox=True, shadow=True)
    plt.show()

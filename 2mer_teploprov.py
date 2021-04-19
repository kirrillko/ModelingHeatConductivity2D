import numpy as np

Nx=20 #число разбиений по координате х
Ny=20 #число разбиений по координате у
xn=0 #начальная точка х
xk=0.1 #конечная точка х
yn=0 #начальная точка у
yk=0.1 #конечная точка у
c_t=700 #теплоёмкость
ro_t=7500 #плотность
lm_t=30 #теплопроводность
time_k=1000 #конечное время расчёта
Tb=300 #начальная температура
Txl=300 #температура на левой границе
Txp=1000 #температура на правой границе
Tyn=500 #температура на нижней границе
Tyv=600 #температура на верхней границе
dt=0.001 #шаг по времени
hx=(xk-xn)/Nx #шаг по пространству в направлении х
hy=(yk-yn)/Ny #шаг по пространству в направлении у
kappa = lm_t/c_t/ro_t
kur_x=dt*kappa/hx/hx
kur_y=dt*kappa/hy/hy
tk=time_k/dt #количество шагов по времени
time=[]
kp=10000

#начальные данных
T = np.zeros(Nx,Ny,tk)
T12 = np.zeros(Nx,Ny,tk)
for i in range(0,Nx):
    for j in range(0,Ny):
        T[i,j,0]=Tb

#граничные условия
for k in range(0,tk):
    for j in range(0,Ny):
        T[0,j,k]=Txl
        T[Nx-1,j,k]=Txp
        T12[0, j, k] = Txl
        T12[Nx - 1, j, k] = Txp
    for i in range(0,Nx):
        T[i,0,k]=Tyn
        T[i,Ny-1,k]=Tyv
        T12[i, 0, k] = Tyn
        T12[i, Ny - 1, k] = Tyv

#коэффициенты прогонки по Х
ax = 1/(hx*hx)
bx = (2/(hx*hx)+1/(0.5*dt))
cx = 1/(hx*hx)
def dx(T,i,j,k):
    return -1/(hy*hy)*T[i,j-1,k]+(2/(hy*hy)-1/(0.5*dt))*T[i,j,k]-1/(hy*hy)*T[i,j+1,k]
def ksiSledX(ksi):
    return cx/(bx-ax*ksi)
def etaSledX(ksi,d,eta):
    return (ax*eta-d)/(bx-ax*ksi)
ksiX = np.zeros(Nx)
ksiX[0]=0
for i in range(1,Nx-1):
    ksiX[i] = ksiSledX(ksiX[i-1])
etaX = np.zeros(Nx)
etaX[0]=0
dX=np.zeros(Nx)

#коэффициенты прогонки по Y
ay = 1/(hy*hy)
by = (1/(hy*hy)+1/(0.5*dt))
cy = 1/(hy*hy)
def dy(T,i,j,k):
    return -1/(hx*hx)*T[i-1,j,k]+(2/(hx*hx)+1/(0.5*dt))*T[i,j,k]-1/(hx*hx)*T[i+1,j,k]
def ksiSledY(ksi):
    return cy/(by-ay*ksi)
def etaSledY(ksi,d,eta):
    return (ay*eta-d)/(by-ay*ksi)
ksiY = np.zeros(Ny)
ksiY[0]=0
for i in range(1,Ny-1):
    ksiY[i] = ksiSledY(ksiY[i-1])
etaY = np.zeros(Ny)
etaY[0]=0
dY=np.zeros(Ny)

#основной цикл
for t in range(1,tk):
    for j in range(1,Ny-1):
        for i in range(0,Nx-2):
            dX[i]=dx(T,i,j,k)
        for i in range(1,Nx-1):
            ksiX[i]=ksiSledX(ksi[i-1])
            etaX[i]=etaSledX(ksiX[i-1],dX[i-1],etaX[i-1])
        for i in range(1,Nx-1):
            T12[i,j,k]=T12[i-1,j,k]*((bx-ax*ksiX[i-1])/cx)+(dX[i-1]-ax*eta[i-1])/cx





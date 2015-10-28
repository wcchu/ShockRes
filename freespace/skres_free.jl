# *** FT convention: g(w) = int{exp(i*w*t)f(t)dt} / sqrt(2*pi)

time0 = time()

afs = 0.024188843
am = 5.2917721E-11
acm = 5.2917721E-9
aum = 5.2917721E-5
lam = 45.5633
aPHz = 6.57969
aTHz = 6.57969E3
aGHz = 6.57969E6
aWcm2 = 7.02E16
c = 137.036

# material parameters
O1 = 80.0/aTHz # resonance strength, Xe: 80THz
w1 = lam/136.0 # resonance position, Xe: 136nm
T1 = 40.0/afs # decay lifetime
K1 = 1.0/T1/2.0 # half width
u1 = sqrt(w1*w1-K1*K1)
O2 = 0.0000259687 #6.0/aTHz # resonance strength
w2 = lam/1925.0 # resonance position
T2 = 40.0/afs # decay lifetime
K2 = 1.0/T2/2.0 # half width
u2 = sqrt(w2*w2-K2*K2)
n2 = 5.8E-19*aWcm2
@printf(STDOUT,"n2 = %14.6F \n",n2)
@printf(STDOUT,"u1 = %14.6F, u2 = %14.6F \n",u1,u2)

# input field
amp = sqrt(1.0E13/aWcm2)
w0 = lam/2000.0
dur = 40.0/afs
phi = 0.0

# dimension parameters
L = 12.0/acm
dzi = 30.0/aum
dzlo = 1.0/aum
dzhi = 0.5/acm
tollo = 7.0E-5
tol0 = 1.0E-4
tolhi = 2.0E-4
nz = 1000
ti = -300.0/afs
tf = 350.0/afs
dt = 0.15/afs
nt = int((tf-ti)/dt+1)
mt = 1
tconv = 10.0/afs
nt2 = int(tconv/dt+1)
tmarg = 20.0/afs # time-window margin
wi = 0.01/aPHz
wf = 0.90/aPHz
dw = 0.001/aPHz
nw = int((wf-wi)/dw+1)
mw = 1
wmarg = 0.05/aPHz # frequency-window margin
@printf(STDOUT,"nt, nw = %8i %8i \n",nt,nw)

# t, w
t = zeros(Float64,nt)
E = zeros(Float64,nt) # E(z,t)
w = zeros(Float64,nw)
Y = zeros(Complex{Float64},nw) # E(z,w)
for it = 1:nt
  t[it] = ti+dt*(it-1)
end
for iw = 1:nw
  w[iw] = wi+dw*(iw-1)
end

# dispersion curve
X1b = zeros(Complex{Float64},nw) # background
X1r = zeros(Complex{Float64},nw) # resonance
nb = zeros(Float64,nw)
n = zeros(Float64,nw)
b0b = zeros(Float64,nw)
b0 = zeros(Float64,nw)
b1b = zeros(Float64,nw)
b1 = zeros(Float64,nw)
b2b = zeros(Float64,nw)
b2 = zeros(Float64,nw)
outn = open("n.dat","w")
for iw = 1:nw
  X1b[iw] = O1*O1/(w1*w1-w[iw]*w[iw]-im*w[iw]*2.0*K1) ### defines dispersion curve
  X1r[iw] = O2*O2/(w2*w2-w[iw]*w[iw]-im*w[iw]*2.0*K2)
  nb[iw] = sqrt(1.0+real(X1b[iw]))
  n[iw] = sqrt(1.0+real(X1b[iw]+X1r[iw]))
  b0b[iw] = w[iw]/c*nb[iw]
  b0[iw] = w[iw]/c*n[iw]
  @printf(outn,"%14.6E %14.6E %14.6E \n",w[iw]*aPHz,nb[iw]-1,n[iw]-1)
end
flush(outn)
close(outn)
np = 1.0
ng = 1.0
outb1 = open("b1.dat","w")
outb2 = open("b2.dat","w")
for iw = 2:nw-1
  b1b[iw] = (b0b[iw+1]-b0b[iw-1])/2.0/dw
  b1[iw] = (b0[iw+1]-b0[iw-1])/2.0/dw
  b2b[iw] = (b0b[iw+1]-2.0*b0b[iw]+b0b[iw-1])/dw/dw
  b2[iw] = (b0[iw+1]-2.0*b0[iw]+b0[iw-1])/dw/dw
  if abs(w[iw]-w0) < 0.5*dw
    np = nb[iw] # phase velocity at pump
    ng = b1b[iw]*c # group velocity at pump
    @printf(STDOUT,"np - 1 = %14.6E \nng - 1 = %14.6E \n",np-1.0,ng-1.0)
  end
  @printf(outb1,"%14.6E %14.6E %14.6E \n",w[iw]*aPHz,b1b[iw]*afs/am,b1[iw]*afs/am)
  @printf(outb2,"%14.6E %14.6E %14.6E \n",w[iw]*aPHz,b2b[iw]*afs*afs/am,b2[iw]*afs*afs/am)
end
@printf(outb1,"# np = %14.6E, ng = %14.6E at pump",np,ng)
flush(outb1)
close(outb1)
flush(outb2)
close(outb2)

# define supergaussian function
sg(x,y) = exp(-0.5*(x/y)^4)
# time window array
twin = ones(Float64,nt)
tt1 = t[1]+tmarg
tt2 = t[nt]-tmarg
for it = 1:nt
  if t[it] < tt1
    twin[it] = sg(t[it]-tt1,tmarg/2.0)
  elseif t[it] > tt2
    twin[it] = sg(t[it]-tt2,tmarg/2.0)
  end
end
# frequency window array
fwin = ones(Float64,nw)
ww1 = w[1]+wmarg
ww2 = w[nw]-wmarg
for iw = 1:nw
  if w[iw] < ww1
    fwin[iw] = sg(w[iw]-ww1,wmarg/2.0)
  elseif w[iw] > ww2
    fwin[iw] = sg(w[iw]-ww2,wmarg/2.0)
  end
end

# input field
tt = dur/sqrt(2.0*log(2.0))
field(x) = amp*exp(-0.5*abs(x/tt)^2)*cos(w0*x+phi)
for it = 1:nt
  E[it] = field(t[it])
end

# preparation for propagation
# 1. define frame velocity, reference refractive index
nf = ng # use group vel. as frame vel.
n0 = np
# 2. coefficients in prop. equation
coeff1 = 0.5*(nf-1.0/nf)/c
coeff2 = -n0/nf*2.0*n2/c
coeff3 = -0.5/nf/c
# 3. convolution function for background dispersion
bb = zeros(Complex{Float64},nw)
B = zeros(Float64,nt2*2+1)
for iw = 1:nw
  bb[iw] = 0.5/nf*w[iw]/c*(1.0+X1b[iw]+X1r[iw]-nf*nf)*fwin[iw]
end
for it = -nt2:nt2
  tempc = 0.0+0.0im
  for iw = 1:nw
    tempc = tempc+exp(-im*w[iw]*dt*it)*bb[iw]*dw
  end
  B[nt2+1+it] = real(im/pi*tempc)
end

# propagation equation
function dEdz(inE)
  dd1 = zeros(Float64,nt)
  dd2 = zeros(Float64,nt)
  dd1[1] = coeff2*inE[1]*inE[1]*(inE[2]-inE[1])/dt
  for it = 2:nt-1
    dd1[it] = coeff2*inE[it]*inE[it]*(inE[it+1]-inE[it-1])/2.0/dt
  end
  dd1[nt] = coeff2*inE[nt]*inE[nt]*(inE[nt]-inE[nt-1])/dt
  for it = 1:nt
    temp = 0.0
    for jt = it-nt2:it+nt2
      if jt >= 1 && jt <= nt
        temp = temp+dt*B[it-jt+nt2+1]*inE[jt]
      end
    end
    dd2[it] = temp
  end
  return dd1+dd2
end

# Runge-Kutta step
function ERK(indz,inE)
  k1 = dEdz(inE)
  k2 = dEdz(inE+0.5*indz.*k1)
  k3 = dEdz(inE+0.5*indz.*k2)
  k4 = dEdz(inE+indz.*k3)
  kk = (k1+2*k2+2*k3+k4)/6
  return inE+indz.*kk
end

# propagation in z with adaptive step
outE = open("E.dat","w")
outY = open("Y.dat","w")
#outdE = open("dE.dat","w")
z = 0.0
dz = dzi
tol = 1.0
for iz = 1:nz
  E = E.*twin
# report E(t) and I(w) at z
  for iw = 1:nw # Ew
    aa = 0+0im
    for it = 1:nt
      aa = aa+dt*exp(im*w[iw]*t[it])*E[it]
    end
    Y[iw] = aa/sqrt(2*pi)
  end
  for it = 1:mt:nt
    @printf(outE,"%14.6F %14.6E %14.6E \n",z*acm,t[it]*afs,E[it])
  end
  @printf(outE," \n")
  for iw = 1:mw:nw
    @printf(outY,"%14.6F %14.6E %14.6E \n",z*acm,w[iw]*aPHz,abs(Y[iw])^2)
  end
  @printf(outY," \n")
  if z > L
    break
  end
# adaptive step
  E0 = E # field before adaptive step
  E1 = ERK(dz,E0) # single step with dz
  Ex = ERK(dz/2.0,E0)
  E2 = ERK(dz/2.0,Ex) # two steps with dz/2
  dE = E2-E1
  tol = norm(dE)/norm(E0)
  if tol > tolhi || tol < tollo # dz has to change; recalculate dE
    dz0 = dz*abs(tol0/tol)^0.2
    if dz0 <= dzlo
      dz = dzlo
    elseif dz0 >= dzhi
      dz = dzhi
    else
      dz = dz0
    end # now the "correct" dz for current E(z,t) has been determined to achieve tol0
    E1 = ERK(dz,E0)
    Ex = ERK(dz/2.0,E0)
    E2 = ERK(dz/2.0,Ex)
    dE = E2-E1
    tol = norm(dE)/norm(E0)
  end
  @printf(STDOUT,"%10.4F %10.2E  tol = %10.2E  cpu_t (min) = %8.2F \n",z*acm,dz*aum,tol,(time()-time0)/60.0)
  #for it = 1:mt:nt
    #@printf(outdE,"%14.6F %14.6E %14.6E %14.6E \n",z*acm,t[it]*afs,E0[it],dE[it])
  #end
  #@printf(outdE," \n")
  E = E2+dE/15.0 # E(z+dz,t) calculated based on corrected dz
  z = z+dz # go to next z
end
flush(outE)
close(outE)
flush(outY)
close(outY)
#flush(outdE)
#close(outdE)

@printf(STDOUT,"total cpu time (min) = %8.2F",(time()-time0)/60.0)

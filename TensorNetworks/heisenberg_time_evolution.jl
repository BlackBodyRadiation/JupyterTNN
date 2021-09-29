using QuantumOptics
using PyPlot

b = SpinBasis(1//2)

sx = sigmax(b)
sy = sigmay(b)
sz = sigmaz(b)

N = 6

spinchain = b ⊗ b

for ispin in 3:N
    global spinchain
    spinchain = spinchain ⊗ b
end

XX = sx ⊗ sx 
YY = sy ⊗ sy 
ZZ = sz ⊗ sz 

J=1

HHeis = embed(spinchain, [1,2] , J*(XX+YY+ZZ))
for ispin in 2:N-1
    global HHeis
    HHeis += embed(spinchain, [ispin,ispin+1] , J*(XX+YY+ZZ))
end

# make sure H is numerically Hermitian
HHeis = ( HHeis + dagger(HHeis) )/2 

s5z = embed(spinchain, [5] , sz)
s6z = embed(spinchain, [6] , sz)

sup = spinup(b)
sdo = spindown(b)

psi0 = sup
for ispin in 2:N
    global psi0
    if ispin == 5
       psi0 = psi0 ⊗ sdo
    else 
       psi0 = psi0 ⊗ sup
    end 
end

times = 0:.1:1.5

tout, psi_t = timeevolution.schroedinger(times, psi0, HHeis)

m5z = expect(s5z, psi_t)
m6z = expect(s6z, psi_t)

plot(tout, m5z)
plot(tout, m6z)

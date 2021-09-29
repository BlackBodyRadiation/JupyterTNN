using QuantumOptics
using PyPlot

b = SpinBasis(1//2)
sx = 1/2*sigmax(b)
sy = 1/2*sigmay(b)
sz = 1/2*sigmaz(b)
XX = sx ⊗ sx 
YY = sy ⊗ sy 
ZZ = sz ⊗ sz

Jxy=1
Jz=1

function basis(N)    
    spinchain=b
    for ispin in 2:N
        spinchain = spinchain ⊗ b
        
    end
    return spinchain # same local
end

function CrtModelH(N)
    Htwoint=Jxy*(XX+YY)+Jz*ZZ
    HXXZ = embed(basis(N), [1,2] , Htwoint)
    for s in 2:N-1
        HXXZ += embed(basis(N), [s,s+1] , Htwoint)
    end
    return HXXZ
end

s5z = embed(basis(10) , [5] , sz)
s6z = embed(basis(10) , [6] , sz)
sup = spinup(b)
sdo = spindown(b)
psi0 = sup
for ispin in 2:10
    global psi0
    if ispin == 5
       psi0 = psi0 ⊗ sdo
    else 
       psi0 = psi0 ⊗ sup
    end 
end

times = 0:.1:12
tout, psi_t = timeevolution.schroedinger(times, psi0, CrtModelH(10))

m5z = expect(s5z, psi_t)
m6z = expect(s6z, psi_t)



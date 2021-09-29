using QuantumOptics

b = SpinBasis(1//2)

sx = sigmax(b)
sy = sigmay(b)
sz = sigmaz(b)

N = 5

spinchain = b ⊗ b

for ispin in 3:N
    global spinchain
    spinchain = spinchain ⊗ b
end

XX = sx ⊗ sx 
YY = sy ⊗ sy 
ZZ = sz ⊗ sz 

Jz=1

HIsing = embed(spinchain, [1,2] , Jz*ZZ)
for ispin in 2:N-1
    global HIsing
    HIsing += embed(spinchain, [ispin,ispin+1] , Jz*ZZ)
end

# make sure H is numerically Hermitian
HIsing = ( HIsing + dagger(HIsing) )/2 

# find first 5 eigenvalues using sparse matrix
E5 = eigenenergies(HIsing, 5)
println(" ")
println("Lowest five eigenvalues")
println(E5)
println(" ")

# full spectrum
EE, UU = eigenstates(dense(HIsing))
println(" ")
println("All Eigenvalues")
println(EE)
println(" ")

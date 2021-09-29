using QuantumOptics

b = SpinBasis(1//2);
sx = 1/2*sigmax(b);
sy = 1/2*sigmay(b);
sz = 1/2*sigmaz(b);
XX = sx ⊗ sx ;
YY = sy ⊗ sy ;
ZZ = sz ⊗ sz;
Jxy=1;
Jz=1;

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

function solvemodl(N)
    E5 = eigenenergies(CrtModelH(N), 5)
    println(" ")
    #println("Lowest five eigenvalues");
    #println(E5);
    println(" ")
    println(" ")
    # full spectrum
    #EE, UU = eigenstates(dense(CrtModelH(N)))
    println(" ")
    #println("All Eigenvalues");
    #println(EE);
    println(" ")
    return E5
end



CrtModelH(2)

dense(CrtModelH(2))

solvemodl(3)

#Solves the for N=7,8,9,10
Nlist = 7:10;
E0s=[]; 
for N in Nlist
    #HM = CrtModelH(N)
    E0 = solvemodl(N)
    push!(E0s,real(E0)/N)
end 
E0pQO=[E0s[i][1] for i=1:length(Nlist)];

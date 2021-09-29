using NBInclude

@nbinclude("GroundState(QOptics).ipynb")

using ITensors

using PyPlot

function HHeisdmrg(N, Jxy, Jz, Nsweep, Bonddim, truncerr)
  sites = siteinds("S=1/2",N)

  ampo = OpSum()
  for j=1:N-1 
    ampo += Jz , "Sz",j,"Sz",j+1
    ampo += Jxy, "Sx",j,"Sx",j+1
    ampo += Jxy, "Sy",j,"Sy",j+1
  end
H = MPO(ampo,sites)
psi0 = randomMPS(sites,10)

  sweeps = Sweeps(Nsweep)
  setmaxdim!(sweeps, Bonddim)
  setcutoff!(sweeps, truncerr)

  energy, psi = dmrg(H,psi0, sweeps)

  return energy, psi
end

Jxy=1
Jz=1
truncerr=1E-12
nsweep = 5
Bonddim= 5
bondlist = 1:2:10
e0s = []

Nlist = 7:10;
E0sIT=[]
for l in Nlist
    e0, psi = HHeisdmrg(l, Jxy, Jz, nsweep, Bonddim, truncerr)
    push!(E0sIT,real(e0)/l)
end 
E0p=[E0sIT[i][1] for i=1:length(Nlist)];

plot(Nlist,E0pQO,color="blue",linestyle = ":",linewidth = 3)
plot(Nlist,E0p,color="red",linewidth = 3)
scatter(Nlist,E0pQO,facecolor = "blue", edgecolors="blue", s=35) 
scatter(Nlist,E0p,facecolor = "red", edgecolors="red", s=35) 
legend(["ED Data(QOptics)", "ED Data(ITensors)","E0(QOptics)", "E0(ITensors)"])
xlabel("N")
ylabel("E0/N")
title("QOptics/ITensors-GS Energies(N)")
savefig( "GS(QOIT).pdf", bbox_inches = "tight", pad_inches = 0.1 )

Nlist = 7:10;
E0s=[]; 
for N in Nlist
    for bonddim in bondlist
        e0, psi = HHeisdmrg(N, Jxy, Jz, nsweep, bonddim, truncerr)
        push!(E0s,real(e0)/N)
    end
    
end 
E0pIT=[E0s[i][1] for i=1:length(Nlist)]

e07=E0s[1:5]
e08=E0s[6:10]
e09=E0s[11:15]
e01=E0s[16:20]

plot(bondlist,e07,color = "black",linewidth = 0.6)
plot(bondlist,e08,color = "green",linewidth = 0.6)
plot(bondlist,e09,color = "blue",linewidth = 0.6)
plot(bondlist,e01,color = "red",linewidth = 0.6)
legend(["E07","E08","E09","E010"])
xlabel("Bond Dim")
ylabel("E0/N")
savefig( "ITensorE0(Bondim).pdf", bbox_inches = "tight", pad_inches = 0.1 )



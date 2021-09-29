using NBInclude


@nbinclude("Time Evolution(QOptics).ipynb")

using ITensors
using PyPlot

  N = 10
  cutoff = 1E-8
  tau = 0.05
  ttotal = 1.5*4

function TEvolHeis(N, cutoff,tau, ttotal)

  # Compute the number of steps to do
  Nsteps = Int(ttotal/tau)

  # Make an array of 'site' indices
  s = siteinds("S=1/2",N)

  # Make gates (1,2),(2,3),(3,4),...
  gates = ITensor[]
  for j=1:N-1
    s1 = s[j]
    s2 = s[j+1]
    hj =       op("Sz",s1) * op("Sz",s2) +
               op("Sx",s1) * op("Sx",s2) +
               op("Sy",s1) * op("Sy",s2)
    Gj = exp(-1.0im * tau/2 * hj)
    push!(gates,Gj)
  end
  # Include gates in reverse order too
  # (N,N-1),(N-1,N-2),...
  append!(gates,reverse(gates))

  c = div(N,2) # center site

  # Initialize psi to be a product state (alternating up and down)
  psi = productMPS(s, n -> n!=c ? "Up" : "Dn")

  Szc=[]
  Szc2=[]
  # Compute and print initial <Sz> value on site c
  t = 0.0
  Sz  = ITensors.expect(psi,"Sz";site_range=c:c)
  Sz2 = ITensors.expect(psi,"Sz";site_range=c+1:c+1)
  println("$t $Sz $Sz2")
  append!(Szc,Sz)
  append!(Szc2,Sz2)

  # Do the time evolution by applying the gates
  # for Nsteps steps and printing <Sz> on site c
  for step=1:Nsteps
    psi = apply(gates, psi; cutoff=cutoff)
    t += tau
    Sz  = ITensors.expect(psi,"Sz";site_range=c:c)
    Sz2 = ITensors.expect(psi,"Sz";site_range=c+1:c+1)
    println("$t $Sz $Sz2")
    append!(Szc,Sz)
    append!(Szc2,Sz2)
  end
  times = 0:tau:ttotal
  return Szc, Szc2, times
end


#Clear demostration of time-evolution
Szc1=TEvolHeis(10, 1E-8, 0.05, 1.5*4)[1];
Szc2=TEvolHeis(10, 1E-8, 0.05, 1.5*4)[2];
timess=TEvolHeis(10, 1E-8, 0.05, 1.5*4)[3];
#Data of longer time-step
Szc1t=TEvolHeis(10, 1E-8, 0.05, 1.5*8)[1];
Szc2t=TEvolHeis(10, 1E-8, 0.05, 1.5*8)[2];
timesst=TEvolHeis(10, 1E-8, 0.05, 1.5*8)[3];
#Data of bigger cut-off
Szc1f=TEvolHeis(10, 1E-4, 0.05, 1.5*4)[1];
Szc2f=TEvolHeis(10, 1E-4, 0.05, 1.5*4)[2];
timessf=TEvolHeis(10, 1E-4, 0.05, 1.5*4)[3];

#Plotting both QOptics and ITensors data of time-evolution 

PyPlot.plot(timess, Szc1)
PyPlot.plot(timess, Szc2)
PyPlot.plot(tout, m5z,linestyle = ":",linewidth = 5)
PyPlot.plot(tout, m6z,linestyle = ":",linewidth = 5)
legend(["m5z(ITensors)", "m6z(ITensors)","m5z(QOptics)","m6z(QOptics)"])
xlabel("Time")
ylabel(L"$\langle S^{Z} \rangle$") 

savefig( "time-evolution.pdf", bbox_inches = "tight", pad_inches = 0.1 )

#Plot of longer Time-step  
PyPlot.plot(timesst, Szc1t)
PyPlot.plot(timesst, Szc2t)
title("Longer Time Scale")
xlabel("Time")
ylabel(L"$\langle  S^{Z}  \rangle$") 
legend(["m5z(ITensors)","m6z(ITensors)"])
savefig( "Longer Time Scale.pdf", bbox_inches = "tight", pad_inches = 0.1 )

#Plot of bigger cut-off
PyPlot.plot(timessf, Szc1f)
PyPlot.plot(timessf, Szc2f)
title(L"Cut-Off($10^{-4}$)")
xlabel("Time")
ylabel(L"$\langle S^{Z} \rangle$") 
legend(["m5z(ITensors)","m6z(ITensors)"])
savefig( "smallercut-off.pdf", bbox_inches = "tight", pad_inches = 0.1 )



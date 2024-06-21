# The main workhorse
using Plots
using DifferentialEquations
using ProgressLogging
using Sundials
include("ecme-dox.jl")

u0 = collect(u0_tup)

param = Params(
  pDOX = DOXParams(DOX=0.0),
  pCK=CKParams(V_ATPASE_CYTO= 1E-5),
  pC1=C1Params(SCALE=1e3 * 50),
  pC3=C3Params(SCALE=60e3 * 15),
  pC4=C4Params(SCALE=60e3 * 15),
  pC5 = C5Params(ρF1=1.5),
  pANT = ANTParams(VMAX=5E-3),
  pSODi = SODParams(ET=1.43E-3 * 0.65),
  BCL=1000.0)

@unpack A_CAP_V_MYO_F = param

tspan = (0.0, 10000.0)
prob = ODEProblem(rhs, u0, tspan, param)

solPacing = solve(prob, CVODE_BDF(); dt=0.01, progress=true, tstops=tspan[1]:1000.0:tspan[end])
pVolt = plot(solPacing, vars=(0, nameLUT[:vm]), label="Membrane Potential", lw=1)
plot!(pVolt, solPacing, vars=(0, nameLUT[:dpsi]), label="Mitochondrial Potential", lw=1, xlabel="time(ms)", ylabel = "Voltage (mV)")
savefig(pVolt, "voltage.png")

plot(solPacing, vars=(0, [nameLUT[:adp_i], nameLUT[:adp_m]]), label=["ADPi" "ADPm"], lw=1)
plot(solPacing, vars=(0, nameLUT[:nadh]), label="[NADH]", lw=1)
plot(solPacing, vars=(0, nameLUT[:ca_m]), label="[Ca]m", lw=1)
plot(solPacing, vars=(0, nameLUT[:isoc]), label="[ISOC]", lw=1)
plot(solPacing, vars=(0, nameLUT[:akg]), label="[AKG]", lw=1)
plot(solPacing, vars=(0, nameLUT[:scoa]), label="[SCoA]", lw=1)
plot(solPacing, vars=(0, nameLUT[:suc]), label="[SUC]", lw=1)
plot(solPacing, vars=(0, nameLUT[:fum]), label="[FUM]", lw=1)
plot(solPacing, vars=(0, nameLUT[:mal]), label="[MAL]", lw=1)
plot(solPacing, vars=(0, nameLUT[:oaa]), label="[OAA]", lw=1)
plot(solPacing, vars=(0, [nameLUT[:Q_n], nameLUT[:QH2_n]]), lw=1, label=["Q_n" "QH2_n"])
plot(solPacing, vars=(0, [nameLUT[:b1], nameLUT[:b2], nameLUT[:b3], nameLUT[:b4]]), lw=1)
plot(solPacing, vars=(0, [nameLUT[:fes_ox], nameLUT[:cytc1_ox], nameLUT[:cytc_ox]]), lw=1, label=["fes_ox" "cytc1_ox" "cytc_ox"])
plot(solPacing, vars=(0, [nameLUT[:sox_i], nameLUT[:sox_m]]), lw=1, label=["sox_i" "sox_m"])
plot(solPacing, vars=(0, nameLUT[:sox_i]), lw=1, label=["sox_i"])
plot(solPacing, vars=(0, nameLUT[:gsh_i]), lw=1)

plot()
plot(t-> ecme_dox(solPacing(t), param, t)[rateMap[:iKatp]], solPacing.t[1], solPacing.t[end], label="iKatp")
plot!(solPacing.t, 8 .- solPacing[nameLUT[:adp_i],:], label="ATPi")
plot!(solPacing, vars=(0, nameLUT[:adp_i]), lw=1, label="ADPi", xlabel="Time (ms)", ylabel=" uA/uF or mM")
savefig("katp.png")

plot!(t-> ecme_dox(solPacing(t), param, t)[rateMap[:vANT]], solPacing.t[1], solPacing.t[end], label="vANT")
plot!(t-> ecme_dox(solPacing(t), param, t)[rateMap[:vATPase]], solPacing.t[1], solPacing.t[end], label="vATPase")
plot!(t-> ecme_dox(solPacing(t), param, t)[rateMap[:vSL]], solPacing.t[1], solPacing.t[end], label="vSL")
plot!(t-> ecme_dox(solPacing(t), param, t)[rateMap[:jUp]], solPacing.t[1], solPacing.t[end], label="jUp")
plot!(t-> A_CAP_V_MYO_F * ecme_dox(solPacing(t), param, t)[rateMap[:iPCa]], solPacing.t[1], solPacing.t[end], label="iPCa")
plot!(t-> A_CAP_V_MYO_F * ecme_dox(solPacing(t), param, t)[rateMap[:iNaK]], solPacing.t[1], solPacing.t[end], label="iNaK")
plot!(t-> ecme_dox(solPacing(t), param, t)[rateMap[:vAM]], solPacing.t[1], solPacing.t[end], label="vAM")
plot!(t-> ecme_dox(solPacing(t), param, t)[rateMap[:d_adp_i]], solPacing.t[1], solPacing.t[end], label="d_adp_i")

plot(sol, vars=(0, nameLUT[:dpsi]), label="Mitochondrial Potential", lw=1)
plot(sol, vars=(0, nameLUT[:nadh]), label="[NADH]", lw=1)
plot(sol, vars=(0, nameLUT[:adp_m]), label="[ADP]m", lw=1)
plot(sol, vars=(0, nameLUT[:isoc]), label="[ISOC]", lw=1)
plot(sol, vars=(0, nameLUT[:akg]), label="[AKG]", lw=1)
plot(sol, vars=(0, nameLUT[:scoa]), label="[SCoA]", lw=1)
plot(sol, vars=(0, nameLUT[:suc]), label="[SUC]", lw=1)
plot(sol, vars=(0, nameLUT[:fum]), label="[FUM]", lw=1)
plot(sol, vars=(0, nameLUT[:mal]), label="[MAL]", lw=1)
plot(sol, vars=(0, nameLUT[:oaa]), label="[OAA]", lw=1)
plot(sol, vars=(0, [nameLUT[:Q_n], nameLUT[:QH2_n]]), lw=1, label=["Q_n" "QH2_n"])
plot(sol, vars=(0, [nameLUT[:b1], nameLUT[:b2], nameLUT[:b3], nameLUT[:b4]]), lw=1)
plot(sol, vars=(0, [nameLUT[:fes_ox], nameLUT[:cytc1_ox], nameLUT[:cytc_ox]]), lw=1, label=["fes_ox" "cytc1_ox" "cytc_ox"])
plot(sol, vars=(0, [nameLUT[:sox_i], nameLUT[:sox_m]]), lw=1, label=["sox_i" "sox_m"])
plot(sol, vars=(0, nameLUT[:sox_i]), lw=1, label=["sox_i"])
plot(sol, vars=(0, nameLUT[:gsh_i]), lw=1)


# Get the rates
plot(t-> ecme_dox(sol(t), param, t)[end-15], sol.t[1], sol.t[end])
plot!(t-> ecme_dox(sol(t), param, t)[end-22], sol.t[1], sol.t[end])
plot(t-> ecme_dox(sol(t), param, t)[end], tspan[1], tspan[end])

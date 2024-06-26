using ProgressLogging
using OrdinaryDiffEq
using ModelingToolkit
using Plots
using DataFrames
using CSV
using ECMEDox
using ECMEDox: second, nernst, mM, Hz
Plots.default(lw=2, size=(800, 800), fmt=:png)

tend = 1000.0second
bcl = 1second
@named sys = build_model(; bcl, tend)
@unpack DOX = sys
sts = unknowns(sys)
obs = [i.lhs for i in observed(sys)]
u0 = build_u0(sys)
# The model is very sensitive to aconitase activity
# The phase transition (of the Q cycle) is between 250uM to 260uM of DOX
prob0 = ODEProblem(sys, u0, tend, [DOX => 0.26mM])
alg = FBDF()
@time sol = solve(prob0, alg; reltol=1e-7, abstol=1e-7, progress=true, maxiters=1e8)

@unpack atp_i, adp_i, vm, na_o, na_i, atp_i, adp_i, atp_m, adp_m, sox_i, sox_m = sys
@unpack ca_i, ca_nsr, ca_jsr, ca_ss = sys
plot(sol, idxs=vm)
p1 = plot(sol, idxs=[ca_i, ca_ss], tspan=(990, tend))
p2 = plot(sol, idxs=[ca_nsr, ca_jsr], tspan=(990, tend))
plot(p1, p2, layout=(2, 1))
@unpack iK1, iK, iKp, iKatp, iNa, iNaCa, iCaL = sys
plot(sol, idxs=[iK1, iK, iKp, iKatp, iCaL])
plot(sol, idxs=[atp_i / adp_i, atp_m / adp_m])

@unpack vHres, vHu, vANT, vHleak, vNaCa, vUni, vIMAC, vROS, vTrROS = sys
@unpack Q_n, Qdot_n, QH2_n, QH2_p, Qdot_p, Q_p, fes_ox, fes_rd, cytc_ox, cytc_rd = sys
plot(sol, idxs=[Q_n, Qdot_n, QH2_n, QH2_p, Qdot_p, Q_p])
plot(sol, idxs=[sys.dpsi * 1000])
plot(sol, idxs=[vHres, vHu, vANT, vHleak, vNaCa, vUni, vTrROS], tspan=(740, 800))
plot(sol, idxs=vANT, tspan=(740, 800))
plot(sol, idxs=[fes_ox / fes_rd], tspan=(400, 500))
plot(sol, idxs=[cytc_ox / cytc_rd], tspan=(400, 500))
plot(sol, idxs=[sys.sox_i, sys.h2o2_i], tspan=(740, 800))
plot(sol, idxs=[sys.vSOD_i, sys.vGPX_i, sys.vGR_i, sys.vCAT], tspan=(400, 500))
plot(sol, idxs=[sys.vROSC1, sys.vROSC3], tspan=(740, 800))

plot(sol, idxs=[sys.vROS / (sys.vO2 + sys.vROS)])
@unpack nad_m, nadh_m, vC1, vIDH, vKGDH, vMDH = sys
plot(sol, idxs=[nad_m / nadh_m])

@unpack cit, isoc, oaa, akg, scoa, suc, fum, mal = sys
@unpack vCS, vACO, vIDH, vKGDH, vSL, vFH, vMDH, vAAT, vSDH = sys
plot(sol, idxs=[cit, isoc, oaa, akg, scoa, suc, fum, mal])

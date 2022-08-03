using OdsIO
using Plots
plotlyjs()
##Loading data files
Bp_s = ods_read("BP.ods", sheetName = "Sheet1")
Bp_b = ods_read("BP.ods", sheetName = "Sheet2")
#Bp_d = ods_read("BP.ods", sheetName = "Sheet3")
Bp_f = ods_read("BP.ods", sheetName = "Sheet4")
##Stretching constants
const D4 = convert(Vector{Float64},Bp_s[3:19,5])
const B4 = convert(Vector{Any},Bp_s[3:19,6])
const R04 = convert(Vector{Any},Bp_s[3:19,7])
const Rf_r2 = convert(Float64,Bp_s[19,8])
const bsw = 1.55
##Bending constants
const KO4 = convert(Vector{Float64},Bp_b[3:40,6])#eV/rad^2
const O04 = Bp_b[3:40,7]#in degrees
const Of4 = Bp_b[3:40,8]
const Of4_abInitio = Bp_b[3:40,9]
const Bp_lab = Bp_b[3:40,4]
##Dihedral bending
const K4d = 4.4 #dihedral constnts doe bicyclopentane eV/rad^2
const T04 = 0.0
##Flaping
const a = convert(Vector{Float64},Bp_f[2:8,2])
const alp_eq = 67.26
##
F_sw(alpha,alpha_eq=alp_eq) = (tanh(bsw*(alpha^2 - alpha_eq^2)*(3.141/180)^2))^2#converting to radin necessary
g(alpha, r_4, r0_4=R04[17], rf_4=Rf_r2, alpha_eq=alp_eq) = r_4 .- r0_4 .- (rf_4 .- r0_4)*F_sw.(alpha, alpha_eq)
##
Vs(r_i,r0_i,D_i,B_i) = D_i*( exp.(-2B_i*(r_i .- r0_i)) -2exp.(-B_i*(r_i .- r0_i)))# r_i is the variable
function V3_stretch(R,Alpha,D=D4,B=B4,R0=R04)
    V3s = 0.0
    for i in 1:length(D)-1
        if D4[i] != 0
                V3s += Vs(R[i], R0[i], D[i], B[i])
        end
    end
    V3s += Vs(g(Alpha,R[17]),0,D[17],B[17])
    return V3s
end
##
function V3d_dum(Alpha,R=R04)
    return V3_stretch(R,Alpha)
end
##
Vb(oi,ki,o0i) = 0.5ki*((oi .- o0i)*pi/180).^2 #pi/180 is used to convert degree to radians
function V3_bend(O,Alpha,O0=O04,K=KO4,Of=Of4,alpha_eq=alp_eq,Bp_label=Bp_lab)
    V3b = 0.0
    for i in 1:length(K)
        if K[i] != 0
            label = parse(Int8,Bp_lab[i][2:3])
            if label < 9 || label > 22
                V3b += Vb(O[i], K[i], O0[i])
            else
                V3b += 0.5*K[i]*( ( O[i] - O0[i] - ( Of[i] - O0[i] )*cos((Alpha/alpha_eq)*pi/2) )*pi/180)^2
            end
        end
    end
    return V3b
end
##
function V3d_dumB(Alpha,O=O04)
    return V3_bend(O,Alpha) + V3_flap(Alpha)
end
##Phi is the C1-C2-C3-C4 dihedral angle
V3_dih(Phi,kd=K4d) = 0.5*kd*Phi^2
##
function V3_flap(Alpha,A = a)
        V3f = 0.0
        i = 0
        for Ai in A
            V3f += Ai*cosd(i*Alpha)
            i += 1
        end
        return V3f
end
## Alpha is the dihedral angles C2-C1-C4-C5
function V3(Alpha, R=R04, O=O04, T=T04)
    V3s = V3_stretch(R, Alpha, D4, B4, R04)
    V3b = V3_bend(O, Alpha, O04, KO4, Of4, alp_eq, Bp_lab)
    Phi = T[1]
    V3d = V3_dih(Phi,K4d)
    V3f = V3_flap(Alpha,a)
    V3 = V3s + V3d + V3f + V3b
    return V3
end

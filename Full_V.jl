#all distances in Angstrom
# in defining the constants everything after '_' is a subscript, numbers not written after '_' is super script
using OdsIO
using Plots
plotlyjs()
const w_1 = 3.8
const w_2 = 3.901
const w_3 = 8.0
const w_4 = 3.0
const w_5 = 2.0
const w_6 = 2.0
const w_7 = 6.0
##
#constants D(eV),B(A^-1) for all the streches  for molecule (1) and (2)
#reading data from spreadsheet
DBH_stretch = ods_read("Stretch.ods", sheetName = "Sheet1")
const D1  = convert(Vector{Float64},DBH_stretch[3:18,4])
const B1  = convert(Vector{Float64},DBH_stretch[3:18,5])
const R01 = convert(Vector{Float64},DBH_stretch[3:18,6])
DB_stretch  = ods_read("Stretch.ods", sheetName = "Sheet2")
const D2  = convert(Vector{Float64},DB_stretch[3:18,4])
const B2  = convert(Vector{Any},DB_stretch[3:18,5])
const R02 = convert(Vector{Float64},DB_stretch[3:18,6])
DB3_stretch = ods_read("Stretch.ods", sheetName = "Sheet3")
const D3  = convert(Vector{Float64},DB3_stretch[3:18,4])
const B3  = convert(Vector{Any},DB3_stretch[3:18,5])
const R03 = convert(Vector{Float64},DB3_stretch[3:18,6])
## Abstracting constants K,O(theta) from file
Bend = ods_read("bending.ods")
const K1  = convert(Vector{Float64}, Bend[3:34,5])
const O01 = convert(Vector{Float64}, Bend[3:34,6])
const K2  = convert(Vector{Float64}, Bend[3:34,7])
const O02 = convert(Vector{Float64}, Bend[3:34,8])
const K3  = convert(Vector{Float64}, Bend[3:34,9])
const O03 = convert(Vector{Float64}, Bend[3:34,10])
##
Dih = ods_read("Dih.ods")
const K1d = convert(Vector{Float64}, Dih[3:14,6])#dihedral force constants
const T01 = convert(Vector{Float64}, Dih[3:14,7])#idicates equilibrium dihedral angles for 1
const K2d = convert(Vector{Float64}, Dih[3:14,8])
const T02 = convert(Vector{Float64}, Dih[3:14,9])#idicates equilibrium dihedral angles for 2
const K3d = convert(Vector{Float64}, Dih[3:14,10])
const T03 = convert(Vector{Float64}, Dih[3:14,11])
const O_1 = convert(Vector{UInt8}, Dih[3:14,12])
const O_2 = convert(Vector{UInt8}, Dih[3:14,13])
##
#pageNo:7916
const B_NN = 2.686035 #in Angstrom inverse
const D_NN = 9.94388 #in ev
const r0_NN = 1.097
##for Bicyclopentane
Bp_s = ods_read("BP_corrected.ods", sheetName = "Sheet1")
Bp_b = ods_read("BP_corrected.ods", sheetName = "Sheet2")
#Bp_d = ods_read("BP.ods", sheetName = "Sheet3")
Bp_f = ods_read("BP_corrected.ods", sheetName = "Sheet4")
##Stretching constants for Bicyclopentane
const D4 = convert(Vector{Float64},Bp_s[3:19,5])
const B4 = convert(Vector{Any},Bp_s[3:19,6])
const R04 = convert(Vector{Any},Bp_s[3:19,7])
const Rf_r2 = convert(Float64,Bp_s[19,8])
const bsw = 1.55
##Bending constants for Bicyclopentane
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
#since S functions(switching functions) have only one susbscript they are given simply as numbers without '_'
S3(r_5, r01_5=R01[5]) = tanh.(w_3*(r_5 .- r01_5).^2)
S4(r_5, r01_5=R01[5]) = exp.(-w_4*(r_5 .- r01_5).^2)
S5(r_6, r02_6=R02[6]) = 1 .- tanh.(w_5*(r_6 .- r02_6))
S6(r_5, r_6, r01_5=R01[5], r01_6=R01[6]) = exp.( -w_6*( (r_5 .- r01_5).^2 .+ (r_6 .- r01_6).^2 ) )
S7(r_5, r02_5=R02[5]) = exp.(-w_7*(r_5 .- r02_5).^2)
##
# The strecthing function for V1
C_i(r_5  ,C1 ,C2 ,r01_5 = R01[5]) = C1 .+ (C2 .- C1).*S3(r_5 ,r01_5)# a single function which can all the jobs of the D_i,B_i,r0_i ...
D_i(r_5 ,D1_i ,D2_i ,r01_5 = R01[5]) = D1_i .+ (D2_i .- D1_i).*S3(r_5 ,r01_5)
B_i(r_5 ,B1_i ,B2_i ,r01_5 = R01[5]) = B1_i .+ (B2_i .- B1_i).*S3(r_5 ,r01_5)
r0_i(r_5 ,r01_i ,r02_i ,r01_5 = R01[5]) = r01_i .+ (r02_i .- r01_i).*S3(r_5 ,r01_5) #0(in r0_i) is supercript indicating the equilibrium bond length
##Stretching
Vs(r_i,r0_i,D_i,B_i) = D_i*( exp.(-2B_i*(r_i .- r0_i)) -2exp.(-B_i*(r_i .- r0_i)))# r_i is the variable
function V1_stretch(R,D=D1,B=B1,R0=R01,n=5)#the full stetching function for V1
    #R is the full bond length vector for all bonds
    #n idicates which bond controls the transition from 1 to 2/3
    # n=5 indicates we are evalutaing the V1 potential between [1] and [2]
    # n=6 idicates transition between [1] and [3]
    D1 = D[1]
    D2 = D[2]# D2/3
    B1 = B[1]
    B2 = B[2]# B2/3
    R01 = R0[1]
    R02 = R0[2]#R02/3
    r5  = R[n] # controlls which bond we monitor
    r01_5 = R01[n]
    V1s = 0.0
    for i in 1:length(D1)
        if D2[i] == 0
            Di = D1[i]
            Bi = B1[i]
            r0i = R01[i]
        else
            Di  = C_i(r5 ,D1[i] ,D2[i] ,r01_5)
            Bi  = C_i(r5 ,B1[i] ,B2[i] ,r01_5)
            r0i = C_i(r5 ,R01[i] ,R02[i] ,r01_5)
        end
        V1s += Vs(R[i] ,r0i ,Di ,Bi)
    end
    return V1s
end
F_sw(alpha,alpha_eq=alp_eq) = (tanh(bsw*(alpha^2 - alpha_eq^2)*(pi/180)^2))^2#converting to radin necessary
g(alpha, r_4, r0_4=R04[17], rf_4=Rf_r2, alpha_eq=alp_eq) = r_4 - r0_4 - (rf_4 - r0_4)*F_sw(alpha, alpha_eq)
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
##Bending
Vb(oi,ki,o0i) = 0.5ki*((oi - o0i)*pi/180)^2 #pi/180 is used to convert degree to radians
function V_bend(O,K,O0,r,r0=[R01[5],R02[6]])
    #r contains both the 5th and 6th bond length
    #r_5 is the parametr that consts the transition between [1] and [2]
    V1b = 0.0
    K1, K2 = K[1], K[2]
    O01, O02 = O0[1], O0[2]
    r_5, r_6 = r[1], r[2]
    r01_5, r02_6 = r0[1], r0[2]
    for i in 1:length(K1)
        if K2[i] == 0
            ki = K1[i]*S4(r_5,r01_5)
            O0i = O01[i]
        else
            if i in [3,4,5,6,9,10,31,32]
                k2i = K2[i]*S5(r_6,r02_6)
            else
                k2i = K2[i]
            end
            ki  = C_i(r_5 ,K1[i] ,k2i ,r01_5)
            O0i = C_i(r_5 ,O01[i] ,O02[i] ,r01_5)
        end
        V1b += Vb(O[i], ki, O0i)
    end
    return V1b
end
function V3_bend(O,Alpha,O0=O04,K=KO4,Of=Of4,alpha_eq=alp_eq,Bp_label=Bp_lab)
    V3b = 0.0
    for i in 1:length(K)
        if K[i] != 0
            label = parse(Int8,Bp_lab[i][2:3])
            if label < 9 || label > 22
                V3b += Vb(O[i], K[i], O0[i])
            else
                V3b += 0.5*K[i]*( pi/180*( O[i] - O0[i] - ( Of[i] - O0[i] )*cos((Alpha/alpha_eq)*pi/2)) )^2
            end
        end
    end
    return V3b
end
##Dihedral Bending
function V_dih(T,O,O_1,O_2,Kd,T0,r,r0=[R01[5],R02[5],R01[6]])
    #T, O, R is the dihedral angle, Angle, and bond length vectors
    #O_1, O_2 are the label vectors for the angles within the dihedral angle
    K1d, K2d = Kd[1], Kd[2]
    T_01, T_02 = T0[1], T0[2]
    r01_5, r02_5, r01_6 = r0[1], r0[2], r0[3]
    r_5, r_6 = r[1], r[2]
    vd1 = 0.0 #Different parts of the dihedral potential
    vd2 = 0.0
    vd3 = 0.0
    for i in 1:8
        if (i < 3)
            vd1 += K1d[i]*(sind(O[O_1[i]])^3)*(sind(O[O_2[i]])^3)*sind(T[i])^2
            #print("vd1=",vd1)
        else
            vd2 += K1d[i]*(sind(O[O_1[i]])^3)*(sind(O[O_2[i]])^3)*(cosd(T[i]) - cosd(T_01[i]))^2
            #print("\tvd2=",vd2)
        end
        vd3 += K2d[i]*(sind(O[O_1[i]])^3)*(sind(O[O_2[i]])^3)*(cosd(T[i]) - cosd(T_02[i]))^2
        #print("\tvd3=",vd3,"\n")
    end
    Vdih = (vd1 + vd2)*S6(r_5,r_6,r01_5,r01_6) + vd3*S7(r_5,r02_5)
    return Vdih
end
V3_dih(Phi,kd=K4d) = 0.5*kd*(Phi*pi/180)^2
##Flapping
function V3_flap(Alpha,A = a)
        V3f = 0.0
        i = 0
        for Ai in A
            V3f += Ai*cosd(i*Alpha)
            i += 1
        end
        return V3f
end
##V1
function V1(Var=[R01, O01, T01])
    R = Var[1]
    O = Var[2]
    T = Var[3]
    D = [D1, D2]
    B = [B1, B2]
    R0 = [R01, R02]
    n = 5# the bond which is treated as breaking
    V1s = V1_stretch(R, D, B, R0, n)
    K = [K1, K2]
    O0 = [O01, O02]
    r = [ R[5], R[6] ]
    r0 = [ R01[5], R02[6] ]
    V1b = V_bend( O, K, O0, r, r0 )
    Kd = [ K1d, K2d ]
    r0 = [ R01[5], R02[5], R01[6] ]
    T0 = [ T01, T02]
    V1d = V_dih( T, O, O_1, O_2, Kd, T0, r, r0 )
    V1 = V1s + V1b + V1d
end
##V2
function V2(Var=[R01, O01, T01])
    R = Var[1]
    O = Var[2]
    T = Var[3]
    D = [D1, D3]
    B = [B1, B3]
    R0 = [R01, R03]
    n = 6
    V2s = V1_stretch(R, D, B, R0, n)
    K = [K1, K3]
    O0 = [O01, O03]
    r = [R[6], R[5]]
    r0 = [R01[6], R03[5]]
    V2b = V_bend( O, K, O0, r, r0 )
    Kd = [ K1d, K3d ]
    r0 = [ R01[6], R03[6], R01[5] ]
    T0 = [ T01, T03]
    V2d = V_dih( T, O, O_1, O_2, Kd, T0, r, r0 )
    V2 = V2s + V2b + V2d
end
##V3
function V3(Var=[alp_eq, R04, O04, T04])
    Alpha = Var[1]
    R = Var[2]
    O = Var[3]
    T = Var[4]
    V3s = V3_stretch(R, Alpha, D4, B4, R04)
    V3b = V3_bend(O, Alpha, O04, KO4, Of4, alp_eq, Bp_lab)
    Phi = T[1]
    V3d = V3_dih(Phi,K4d)
    V3f = V3_flap(Alpha,a)
    V3 = V3s + V3d + V3f + V3b
    return V3
end
##V4
V4(r) = D_NN*(exp.(-2B_NN*(r-r0_NN)) .- 2exp.(-B_NN*(r-r0_NN)))
##
A, B, C = 28, 0.98, 1.52
S1(r_5, r_6, r01_5=R01[5], r01_6 = R01[6]) = ( 1.0 - tanh( w_1*(r_5 -r01_5)^2 )*tanh( w_1*(r_6 - r01_6)^2 ) )*( 1 - B*(r_5 - r01_5)^2*(r_6 - r01_6)^2 )*exp( -C*( (r_5 - r01_5)^2 + (r_6 - r01_6)^2 ) )
S2(r_5, r_6, r01_5=R01[5], r01_6 = R01[6]) = tanh( w_2*(r_5 - r01_5)^2 )*tanh( w_2*(r_6 - r01_6)^2 )
##
function V(Var=[R01, O01, T01, alp_eq])
    R = Var[1]
    O = Var[2]
    T = Var[3]
    Alpha = Var[4]
    v1 = V1([R, O, T])
    v2 = V2([R, O, T])
    v3 = V3([Alpha, R, O, T])
    v4 = V4(R[16])
    s1 = S1(R[5],R[6],R01[5],R01[6])
    s2 = S2(R[5],R[6],R01[5],R01[6])
    #print("\n$s1\t$s2\n")
    V = ( (v1 + v2*exp(-A*(v2 - v1)))/(1 + exp(-A*(v2 - v1))) )*s1 + (v3 + v4)*s2
    return V
end
##
#dummy Varaibles
R = zeros(17)
R[1:16] = R02
O = zeros(38)
O[1:32] = O02
O[33:38] = O04[33:38]

#all distances in Angstrom
# in defining the constants everything after '_' is a subscript, numbers not written after '_' is super script
using OdsIO
using Plots
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
const B2  = convert(Vector{Float64},DB_stretch[3:18,5])
const R02 = convert(Vector{Float64},DB_stretch[3:18,6])
DB3_stretch = ods_read("Stretch.ods", sheetName = "Sheet3")
const D3  = convert(Vector{Float64},DB3_stretch[3:18,4])
const B3  = convert(Vector{Float64},DB3_stretch[3:18,5])
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
Vs(r_i,r0_i,D_i,B_i) = D_i*( exp.(-2B_i*(r_i .- r0_i)) -2exp.(-B_i*(r_i .- r0_i)))# r_i is the variable
##
function V1_stretch(R,D,B,R0,n=5)#the full stetching function for V1
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
Vb(oi,ki,o0i) = 0.5ki*((oi .- o0i)*pi/180).^2 #pi/180 is used to convert degree to radians

function V_bend(O,K,O0,r,r0=[R01[5],R02[6]])
    # r contains both the 5th and 6th bond length
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
##
function V_dih(T,O,O_1,O_2,Kd,T0,r,r0=[R01[5],R02[5],R01[6]])
    #T, O, R is the dihedral angle, Angle, and bond length vectors
    #O_1, O_2 are the label vectors for the angles within the dihedral angle
    K1d, K2d = Kd[1], Kd[2]
    T01, T02 = T0[1], T0[2]
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
            vd2 += K1d[i]*(sind(O[O_1[i]])^3)*(sind(O[O_2[i]])^3)*(cosd(T[i]) - cosd(T01[i]))^2
            #print("\tvd2=",vd2)
        end
        vd3 += K2d[i]*(sind(O[O_1[i]])^3)*(sind(O[O_2[i]])^3)*(cosd(T[i]) - cosd(T02[i]))^2
        #print("\tvd3=",vd3,"\n")
    end
    Vdih = (vd1 + vd2)*S6(r_5,r_6,r01_5,r01_6) + vd3*S7(r_5,r02_5)
    return Vdih
end
##
V4(r) = D_NN*(exp.(-2B_NN*(r-r0_NN)) .- 2exp.(-B_NN*(r-r0_NN)))
##
function V1(R, O, T)
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
##
function V2(R, O, T)
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

##
V1(R03,O03,T01)

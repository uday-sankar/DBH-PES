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
## Abstracting constants K,O(theta) from file
Bend = ods_read("bending.ods")
const K1  = convert(Vector{Float64}, Bend[3:34,5])
const O01 = convert(Vector{Float64}, Bend[3:34,6])
const K2  = convert(Vector{Float64}, Bend[3:34,7])
const O02 = convert(Vector{Float64}, Bend[3:34,8])
##
Dih = ods_read("Dih.ods")
const K1d = convert(Vector{Float64}, Dih[3:14,6])#dihedral force constants
const T01 = convert(Vector{Float64}, Dih[3:14,7])#idicates equilibrium dihedral angles for 1
const K2d = convert(Vector{Float64}, Dih[3:14,8])
const T02 = convert(Vector{Float64}, Dih[3:14,9])#idicates equilibrium dihedral angles for 2
const O_1 = convert(Vector{UInt8}, Dih[3:14,10])
const O_2 = convert(Vector{UInt8}, Dih[3:14,11])
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
function V1_stretch(R,D,B,R0)#the full stetching function for V1
    #R is the full bond length vector for all bonds
    D1 = D[1]
    D2 = D[2]
    B1 = B[1]
    B2 = B[2]
    R01 = R0[1]
    R02 = R0[2]
    r5  = R[5]
    r01_5 = R01[5]
    V1s = 0.0
    for i in 1:length(R)
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
##
function V1_stretchM(Rn,D1,D2,B1,B2,R01,R02)
    Vr = zeros(length(Rn))
    i = 1
    for r in Rn
        Vr[i] = V_stretch(r,D1,D2,B1,B2,R01,R02)
        i += 1
    end
    return Vr
end
##

Vb(o ,o0 ,ki::Float64) = 0.5ki*(o .- o0).^2 #bending modes; ki is a float or array

function V_bend(O,r_5,K,O0,r01_5=R01[5])
    V1b = 0.0
    K1 = K[1]
    K2 = K[2]
    O01 = O0[1]
    O02 = O0[2]
    for i in 1:length(O)
        if K2[i] == 0
            ki = K1[i]*S4(r_5,r01_5)
            O0i = O10[i]
        else
            ki  = C_i(r_5 ,K1[i] ,K2[i] ,r01_5)
            O0i = C_i(r_5 ,O10[i] ,O10[i] ,r01_5)
        end
        V1b += Vb(O[i] ,O0i ,ki)
    end
    return V1b
end
##
function V_dih(T,R,O,O_1,O_2,K1d,K2d,T01,T02)
    #T, O, R is the dihedral angle, Angle, and bond length vectors
    #O_1, O_2 are the label vectors for the angles within the dihedral angle
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
    Vdih = (vd1 + vd2)*S6(R[5],R[6]) + vd3*S7(R[5])
    return Vdih
end

##
V4(r) = D_NN*(exp.(-2B_NN*(r-r0_NN)) .- 2exp.(-B_NN*(r-r0_NN)))

##
x=1.4:0.05:34
#y=Vs(x,1.5,1,1)

plotlyjs()
##
#plot(x,Vb(x,10,4.0))
r1 = 0.9:0.001:2
r2 = 1.15:0.001:2
x = []
y = []
V = []
R = copy(R01)
for i in 1:length(r1)
    for j in 1:length(r2)
        R[1] = copy(r1[i])
        R[2] = copy(r2[j])
        append!(x ,r1[i])
        append!(y ,r2[j])
        append!(V ,V1_stretch(R,D1,D2,B1,B2,R01,R02))
    end
end
##
function Vp_bend(r1,r2=102,p1=pos1,p2=pos2)
    r_5 = 1.503
    R = copy(O01)
    R[p1] = r1
    R[p2] = r2
    V_bend(R,r_5,K1,K2,O01,O02,R01[5])
end

function Vp_str(r1,r2=R01[2],p1=pos1,p2=pos2)
    R = copy(R01)
    R[p1] = r1
    #R[p2] = r2
    V1_stretch(R,[D1,D2],[B1,B2],[R01,R02])
end
o1 = 90:1:110
o2 = 50:1:150
pos1 = 1
pos2 = 3
r1 = 1:0.05:3
r2 = 1:0.05:4
#plotlyjs()
#p = plot(r2 ,Vp_str)
plot(o1 ,Vp_bend, xaxis = (label = "x"))
##
function Vp_dih(r1,r2,p1=pos1,p2=pos2)
    T = copy(T01)
    O = copy(O01)
    T[p1] = r1
    T[p2] = r2
    V_dih(T,R01,O,O_1,O_2,K1d,K2d,T01,T02)
end
d1 = 0:1:360
pos1 = 1
pos2 = 9
#surface(d1,d1 ,Vp_dih, xaxis = (label = "x"))

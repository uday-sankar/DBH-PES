#optimize(V3, [30,ones(17), ones(38), 40])
using Statistics
function Derivative(F,Variables=[alp_eq,R04,O04,T04],dx=0.01,pos=2)
    VarAr = copy(Variables)
    val_init = F(VarAr)
    #print("Initial Value = $val_init")
    Var = copy(VarAr[pos])
    if Var isa Array
        dF = zeros(length(Var))
        for i in 1:length(Var)
            if isa(Var[i], Number)
                #print("\nPOS:$i")
                var_1 = Var[i]
                var_0 = var_1 - dx #The previous point
                var_2 = var_1 + dx #The next point
                Var[i] = var_2
                VarAr[pos] = Var
                #print("\nVariable $(Var[i])")
                dFf = (F(VarAr) -val_init)
                #print("\ndF = $dFf")
                Var[i] = var_0
                VarAr[pos] = Var
                #print("\nVariable $(Var[i])\n")
                dFb = (val_init - F(VarAr) )
                #print("\ndF = $dFb")
                dF[i] = (dFf+dFb)/(2dx)
                Var[i] = var_1
                VarAr[pos] = Var
            end
        end
    elseif Var isa Number
        dF = 0.0
        var_1 = Var
        var_0 = var_1 - dx #The previous point
        var_2 = var_1 + dx #The next point
        Var = var_2
        VarAr[pos] = Var
        #print("\nVariable $(Var)")
        dFf = (F(VarAr) -val_init)
        #print("\ndF = $dFf")
        Var = var_0
        VarAr[pos] = Var
        #print("\nVariable $(Var)\n")
        dFb = (val_init - F(VarAr) )
        #print("\ndF = $dFb")
        dF = (dFf+dFb)/(2dx)
        Var = var_1
        VarAr[pos] = Var
    else
        dF = 0
    end
    return dF
end
##
function Optimize(F,Variables=[alp_eq,R04,O04,T04];Threshold=0.01,counter=3,dx=0.1,pos=1)
    OptVar = copy(Variables[pos])
    VarAr = copy(Variables)
    V_init = F(VarAr)
    V_last = V_init
    dF = Derivative(F, Variables,dx,pos)
    c = 0
    if OptVar isa Number
        if dF == 0
            print("dF=0\n")
            OptVar = OptVar .+ dx
            VarAr[pos] = OptVar
            dF = Derivative(F,VarAr,dx,pos)
            if dF > 0
                print("Optimization reached")
                return Variables[pos]
            end
        end
        #print("Initial:$V_init\t$OptVar\t$dx\t$dF\n")
        while abs(dF)>Threshold
            if c > counter
                break
            end
            #print("\n$OptVar")
            OptVar = OptVar .- dx*sign(dF)
            VarAr[pos] = OptVar
            V_n = F(VarAr)
            if V_n > V_last
                #print("\nTrue")
                OptVar = OptVar .+ dx*sign(dF)
                VarAr[pos] = OptVar
                dx = dx/10
                c += 1
            end
            V_last = V_n
            dF = Derivative(F,VarAr,dx,pos)
            print("\ndF \t$dF")
        end
    elseif OptVar isa Array
        while maximum(abs.(dF)) > Threshold
            if c > counter
                break
            end
            #print("\n$OptVar\n$dF")
            OptVar = OptVar .- dx*dF/maximum(abs.(dF))
            VarAr[pos] = OptVar
            V_n = F(VarAr)
            if V_n > V_last
                #print("\nTrue")
                OptVar = OptVar .+ dx*dF./maximum(abs.(dF))
                VarAr[pos] = OptVar
                dx = dx/10
                c += 1
            end
            V_last = V_n
            dF = Derivative(F,VarAr,dx,pos)
            print("\ndF \t$(mean(dF))")
        end
    end
    print("\nOptimization done\n")#\nOptimized Value:$OptVar\n")
    return OptVar
end
##
function Flap_p(flap_a)
    V = zeros(length(flap_a))
    i = 1
    R0 = replace(R04, "..."=>0)
    O0 = replace(O04, "..."=>0)
    T0 = T04
    print("\nStarted")
    for f in flap_a
        print("\n Flap_angle = $f")
        R0_i = Optimize(V3,[f,R0,O0,T0]; Threshold = 1e-5,counter =4 , dx = 0.1, pos = 2)
        print("\nR0 Optimization complete")
        O0_i = Optimize(V3,[f,R0,O0,T0]; Threshold = 1e-5,counter =2 , dx = 0.1, pos = 3)
        print("\nO0 Optimization complete")
        T0_i = Optimize(V3,[f,R0,O0,T0]; Threshold = 1e-5,counter =2 , dx = 0.1, pos = 4)
        print("\nT0 Optimization complete\n")
        V[i] = V3(f,R0_i,O0_i,T0_i)
        i += 1
    end
    return V
end

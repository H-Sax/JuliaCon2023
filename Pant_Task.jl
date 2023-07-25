using DifferentialEquations, ModelingToolkit, Plots,  CirculatorySystemModels

begin
    τ = 1.0
    ### Shi params for ventricle, atrium 
    v0_rv = 10.0
    p0_rv = 1.0
    Emin_rv = 0.1
    Emax_rv = 1.15
    τes_rv = 0.3
    τed_rv = 0.45
    Eshift_rv = 0.0
    v0_ra = 4.0
    p0_ra = 1.0
    Emin_ra = 0.15
    Emax_ra = 0.25
    τpwb_ra = 0.92
    τpww_ra = 0.09
    τes_ra = τpww_ra/2
    τed_ra = τpww_ra
    Eshift_ra = τpwb_ra

    ## PAITENT A PARAMETERS WITH NO REGURITATION 
    # Paitent A is tricuspid and Pulmonary valve 
    # Valve parameters
    ρ = 1.06
    #AOV # Pulmonary
    Amax_Aov = 1.57
    Amin_Aov = 0.0 + eps()
    Kvo_P = 0.02
    Kvc_P = 0.02
    Leff_P = 1
    #AVV # tricuspid 
    Amax_Avv = 1.1
    Amin_Avv = 0.0 + eps()
    Kvo_T = 0.03
    Kvc_T = 0.04
    Leff_T = 0


    #Sarcome chamber 
    l0 = 1.9
    la0 = 1.65
    lam = 2.2
    lae = 2.4
    laf = 2.95
    ν0 = 10
    cp = 12
    Ta0 = 55
    Tp0 = 0.9
    cv = 0
    t1 = 0.05
    Tc = 0.51
    ϕ_ = 0.275
    μ_ = 0.3
    # SV
    Vw_SV = 33.9
    V0_SV = 12.4
    Ea_SV = 0.5
    ta_SV = 0.31
    # SA
    Vw_SA = 1.3
    V0_SA = 6.37
    Ea_SA = 1.42
    ta_SA = 0.15
    # Heart - To go in here will be valve and chamber params 
    C_AO = 4.6e-02

    # Lower Body 
    R1_LB = 1.08
    K_LB = 5.0e-03
    L1_LB = 0.027
    C1_LB = 0.14
    R2_LB = 6.62
    C2_LB = 3.09
    R3_LB = 0.18

    # Lungs 
    R_SH = 0.79
    K_SH = 0.25 
    C_SH = 1.0e-04
    R1_LU = 7.9e-03
    L1_LU = 0.4
    C1_LU = 0.11
    R2_LU = 0.54
    C2_LU = 0.06
    R3_LU = 0.23

    # Upper Body 
    R1_UB = 0.6
    K_UB = 3.3e-03
    L1_UB = 4.0e-04
    C1_UB = 0.11
    R2_UB = 3.26
    C2_UB = 0.12
    R3_UB = 0.76

    ## Paitent B parameters Left ventricle 
    # Shi parameters for ventricle and Atrium 
    v0_lv = 5.0
    p0_lv = 1.0
    Emin_lv = 0.1
    Emax_lv = 2.5
    τes_lv = 0.3
    τed_lv = 0.45
    Eshift_lv = 0.0
    v0_la = 4.0
    p0_la = 1.0
    Emin_la = 0.15
    Emax_la = 0.25
    τpwb_la = 0.92
    τpww_la = 0.09
    τes_la = τpww_la/2
    τed_la = τpww_la
    Eshift_la = τpwb_la
    #AOV # Aortic
    Amax_AovB = 0.92
    Amin_AovB = 0.0 + eps()
    Kvo_A =  0.012
    Kvc_A =  0.012
    Leff_A = 1
    #AVV # Mitral
    Amax_AvvB = 2.64
    Amin_AvvB = 0.0 + eps()
    Armax_AvvB = 0.036
    Kvo_M = 0.03
    Kvc_M = 0.04
    Leff_M = 0

    #Sarcome chamber 
    l0 = 1.9
    la0 = 1.65
    lam = 2.2
    lae = 2.4
    laf = 2.95
    ν0 = 10
    cp = 12
    Ta0 = 55
    Tp0 = 0.9
    cv = 0
    t1 = 0.05
    Tc = 0.51
    ϕ_ = 0.275
    μ_ = 0.3
    # SV
    Vw_SV_B = 40.0
    V0_SV_B = 10.84
    Ea_SV_B = 0.48
    ta_SV_B = 0.32
    # SA
    Vw_SA_B = 2.55
    V0_SA_B = 11.41
    Ea_SA_B = 1.72
    ta_SA_B = 0.26
    # Heart - To go in here will be valve and chamber params 
    C_AO_B = 2.2e-02

    # Lower Body 
    R1_LB_B = 3.42
    K_LB_B= 0.03
    L1_LB_B = 5.5e-03
    C1_LB_B = 0.045
    R2_LB_B = 15.24
    C2_LB_B = 3.35
    R3_LB_B = 0.044

    # Lungs 
    R_SH_B = 0.79 # no valve use paper 3
    K_SH_B = 0.25 # no value use paper 3
    C_SH_B = 0.02
    R1_LU_B = 0.11
    L1_LU_B = 0.045
    C1_LU_B = 0.11
    R2_LU_B = 0.21
    C2_LU_B = 0.64
    R3_LU_B = 0.08

    # Upper Body 
    R1_UB_B = 0.4
    K_UB_B = 1.97e-02
    L1_UB_B = 1.06e-02
    C1_UB_B = 0.1
    R2_UB_B = 6.5
    C2_UB_B = 0.3
    R3_UB_B = 0.16
end

function PantValve_SemiLunar(; name, ρ, Leff, Amax, Amin, Kvc, Kvo, Δpr, Kvc_r, Kvo_r, Armax)
    @named oneport = OnePort()
    @unpack Δp, q = oneport
    ps = @parameters ρ = ρ Leff = Leff Amax = Amax Amin = Amin Kvc = Kvc Kvo = Kvo Δpr=Δpr Kvc_r=Kvc_r Kvo_r=Kvc_r Armax=Armax
    sts = @variables Aeff(t) = 0.0 ζ(t) = 0.0 B(t) = 0.0 L(t) = 0.0 
    D = Differential(t) 
    Δp = -1333.22*Δp
    # make θmax the real opening angle and define a θmaxopen for a healthy valve
    # that means we can use θmax as a stenosis parameter
    eqs = [                      
        # Opening ratio
        D(ζ) ~   (Δp ≥ 0)*((1-ζ)*Kvo*Δp) + (Δp < 0)*(ζ ≥ 0)*((ζ*Kvc*Δp)) + (Δp ≤ Δpr)*((1+ζ)*Kvo_r*(Δp-Δpr)) + (Δp > Δpr)*(ζ < 0)*(-ζ*Kvc_r*(Δp-Δpr))
        Aeff ~ ifelse(ζ ≥ 0, ((Amax - Amin)*ζ + Amin) , ((-Armax - Amin)*ζ + Amin))
        # Flow equation
        B ~ ρ/(2*Aeff^2)  
        L ~ ρ*Leff/Aeff
        D(q) ~ (Δp - B*q*abs(q))*1/L
        
    ]
    extend(ODESystem(eqs, t, sts, ps; name=name), oneport)
end

function PantValve_Atrioventricular(; name, ρ, Amax, Amin, Kvc, Kvo, Δpr, Kvc_r, Kvo_r, Armax)
    @named oneport = OnePort()
    @unpack Δp, q = oneport
    ps = @parameters ρ = ρ Amax=Amax  Amin = Amin Kvc = Kvc Kvo = Kvo Δpr = Δpr Kvc_r=Kvc_r Kvo_r=Kvo_r Armax=Armax
    sts = @variables Aeff(t) = 0.0 ζ(t) = 0.0 B(t) = 0.0 Aeff_min(t) = 0.0 Aeff_max(t) = 0.0 L(t) = 0.0 
    D = Differential(t) 
    Δp = -1333.22*Δp
    # make θmax the real opening angle and define a θmaxopen for a healthy valve
    # that means we can use θmax as a stenosis parameter
    eqs = [                      
        # Opening ratio
        D(ζ) ~   (Δp ≥ 0)*((1-ζ)*Kvo*Δp) + (Δp < 0)*(ζ ≥ 0)*((ζ*Kvc*Δp)) + (Δp ≤ Δpr)*((1+ζ)*Kvo_r*(Δp-Δpr)) + (Δp > Δpr)*(ζ < 0)*(-ζ*Kvc_r*(Δp-Δpr))
        Aeff ~ ifelse(ζ ≥ 0, ((Amax - Amin)*ζ + Amin) , ((-Armax - Amin)*ζ + Amin))
        # Flow equation
        B ~ ρ/(2*Aeff^2)  
        q ~ sqrt(1/B*abs(Δp))*sign(Δp)
        
    ]
    extend(ODESystem(eqs, t, sts, ps; name=name), oneport)
end

@parameters t
@named SV = ShiChamber(V₀=v0_rv, p₀ = p0_rv, Eₘᵢₙ=Emin_rv, Eₘₐₓ=Emax_rv, τ=τ, τₑₛ=τes_rv, τₑₚ=τed_rv, Eshift=0.0)
## The atrium can be defined either as a ShiChamber with changed timing parameters, or as defined in the paper
@named SA = ShiChamber(V₀=v0_ra, p₀ = p0_ra, Eₘᵢₙ=Emin_ra, Eₘₐₓ=Emax_ra, τ=τ, τₑₛ=τpww_ra/2, τₑₚ =τpww_ra, Eshift=τpwb_ra)

@named AOV = PantValve_SemiLunar(ρ = ρ, Leff = Leff_P, Amax = Amax_Aov, Amin = Amin_Aov, Kvc = Kvc_P, Kvo = Kvo_P, Δpr = 0.0, Kvc_r = 0.0, Kvo_r = 0.0, Armax = 0.0) #Pulmonary
@named AVV = PantValve_Atrioventricular(ρ = ρ, Amax = Amax_Avv, Amin = Amin_Avv, Kvc = Kvc_T, Kvo = Kvo_T, Δpr = 0.0, Kvc_r = 0.0, Kvo_r = 0.0, Armax = 0.0) #Tricuspid

@named Cao = Compliance(C=C_AO)

# Lower body 
@named R1_lb = Resistor(R=R1_LB)
@named K_lb = QResistor(K=K_LB)
@named L1_lb = Inductance(L=L1_LB)
@named C1_lb = Compliance(C=C1_LB)

@named R2_lb = Resistor(R=R2_LB)
@named C2_lb = Compliance(C=C2_LB)

@named R3_lb = Resistor(R=R3_LB)

#Lungs 
@named R_sh = Resistor(R=R_SH) 
@named K_sh = QResistor(K=K_SH)
@named C_sh = Compliance(C=C_SH)

@named R1_lu = Resistor(R=R1_LU)
@named L1_lu = Inductance(L=L1_LU)
@named C1_lu = Compliance(C=C1_LU)

@named R2_lu = Resistor(R = R2_LU)
@named C2_lu = Compliance(C = C2_LU)

@named R3_lu = Resistor(R = R3_LU)

# Uppper Body 
@named R1_ub = Resistor(R=R1_UB)
@named K_ub = QResistor(K=K_UB)
@named L1_ub = Inductance(L=L1_UB)
@named C1_ub = Compliance(C=C1_UB)

@named R2_ub = Resistor(R=R2_UB)
@named C2_ub = Compliance(C=C2_UB)

@named R3_ub = Resistor(R=R3_UB)


circ_eqs = [
    #Vnetricle out 
    connect(SV.out, AOV.in)
    connect(AOV.out, Cao.in)
    connect(Cao.out,R1_lb.in, R_sh.in, R1_ub.in)
    # lower body loop
    connect(R1_lb.out, K_lb.in)
    connect(K_lb.out, L1_lb.in)
    connect(L1_lb.out,C1_lb.in)
    connect(C1_lb.out, R2_lb.in)
    connect(R2_lb.out, C2_lb.in)
    connect(C2_lb.out, R3_lb.in)
    connect(R3_lb.out, SA.in)
    # Lungs Loops
    connect(R_sh.out, K_sh.in)
    connect(K_sh.out, C_sh.in)
    connect(C_sh.out, R1_lu.in)
    connect(R1_lu.out, L1_lu.in)
    connect(L1_lu.out, C1_lu.in)
    connect(C1_lu.out, R2_lu.in)
    connect(R2_lu.out, C2_lu.in)
    connect(C2_lu.out, R3_lu.in)
    connect(R3_lu.out, SA.in)
    ## Upper body loop
    connect(R1_ub.out, K_ub.in)
    connect(K_ub.out, L1_ub.in)
    connect(L1_ub.out,C1_ub.in)
    connect(C1_ub.out, R2_ub.in)
    connect(R2_ub.out, C2_ub.in)
    connect(C2_ub.out, R3_ub.in)
    connect(R3_ub.out, SA.in)
    # Heart Close loop
    connect(SA.out, AVV.in)
    connect(AVV.out, SV.in)
]

## Compose the whole ODAE system
@named _circ_model = ODESystem(circ_eqs, t)

##
@named circ_model = compose(_circ_model,
    [SA, AVV, SV, AOV, Cao, R1_lb, K_lb, L1_lb, C1_lb, R2_lb, C2_lb, R3_lb, R_sh, K_sh, C_sh, R1_lu, L1_lu, C1_lu, R2_lu, C2_lu, R3_lu, R1_ub, K_ub, L1_ub, C1_ub, R2_ub, C2_ub, R3_ub])
    
## And simplify it
@time circ_sys = structural_simplify(circ_model)

u0 = [0.0, 0.0, 0.0, 50.0, 0.0, 40.0, 30.0, 40.0, 0.0, 40.0, 30.0, 0.0, 40.0, 30.0, 8.0, 20.0] # IC for Shi chamber and Mynard valves


@time prob = ODAEProblem(circ_sys,u0, (0.0, 20.0))
##
@time sol = solve(prob, AutoVern7(Rodas4()), reltol = 1e-8, abstol = 1e-8)


plot(sol, idxs = [SA.V, SV.V], tspan = (10*τ, 15*τ))
plot(sol, idxs = [AVV.q, AOV.q], tspan = (10*τ, 15*τ)) 
plot(sol, idxs = [AVV.ζ, AOV.ζ], tspan = (10*τ, 15*τ)) 
plot(sol, idxs = [Cao.p, Cao.in.q], tspan = (17*τ, 18*τ)) 
plot(sol, idxs = [C1_ub.in.q,C1_lb.in.q, C1_lu.in.q], tspan = (17*τ, 18*τ)) 
plot(sol, idxs = [SV.p], tspan = (17*τ, 18*τ)) 
plot(sol, idxs = [SA.p], tspan = (17*τ, 18*τ)) 
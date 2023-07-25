using CirculatorySystemModels, ModelingToolkit, DifferentialEquations, Plots

## Go through basic equations of CirculatorySystemModels ## 

## Two element Windkessel ##

@parameters t

begin
	function heaviside(t)
	   0.5 * (sign(t) + 1)
	end
	
	function I_sin(t)
		sin(2pi*(t))^2 * (heaviside(mod(t,1)) - heaviside(mod(t,1)-0.5))
	end
end

plot(I_sin, xlim=[0,2], label="I", xlabel="time", ylabel="I [ml/s]")

@named source = DrivenFlow(Q=-100.0, fun=I_sin)
@named R = Resistor(R = 1.0)
@named C = Capacitor(C = 1.3)
@named ground = Ground()

circ_eqs = [
    connect(source.out, C.in, R.in)
    connect(R.out, source.in, C.out, ground.g) 
]

## Compose the whole ODAE system
@named _circ_model = ODESystem(circ_eqs, t)

##
@named circ_model = compose(_circ_model,[source, ground, C, R])

## And simplify it
@time circ_sys = structural_simplify(circ_model)
# Examine equation 

# give mean pressures as IC
u0 = [10.0]
@time prob = ODAEProblem(circ_sys, u0, (0.0, 10.0))

@time sol = solve(prob)

plot(sol)


# Two element Windkessel compliance #
@parameters t

begin
	function heaviside(t)
	   0.5 * (sign(t) + 1)
	end
	
	function I_sin(t)
		sin(2pi*(t))^2 * (heaviside(mod(t,1)) - heaviside(mod(t,1)-0.5))
	end
end

plot(I_sin, xlim=[0,2], label="I", xlabel="time", ylabel="I [ml/s]")

@named source = DrivenFlow(Q=100.0, fun=I_sin)
@named R = Resistor(R = 1.0)
@named C = Compliance(C = 1.3)
@named ground = Ground()

circ_eqs = [
    connect(source.out, C.in) 
	connect(C.out, R.in)
    connect(R.out, source.in, ground.g) 
]

## Compose the whole ODAE system
@named _circ_model = ODESystem(circ_eqs, t)

##
@named circ_model = compose(_circ_model,[source, ground, C, R])

## And simplify it
@time circ_sys = structural_simplify(circ_model)
# Examine Equation

# give mean pressures as IC
u0 = [9.0]
@time prob = ODAEProblem(circ_sys, u0, (0.0, 10.0))

@time sol = solve(prob)

plot(sol)

# Three element windkessel #
@parameters t

begin
	function heaviside(t)
	   0.5 * (sign(t) + 1)
	end
	
	function I_sin(t)
		sin(2pi*(t))^2 * (heaviside(mod(t,1)) - heaviside(mod(t,1)-0.5))
	end
end

plot(I_sin, xlim=[0,2], label="I", xlabel="time", ylabel="I [ml/s]")

@named source = DrivenFlow(Q=100.0, fun=I_sin)
@named Rp = Resistor(R = 0.1)
@named Rc = Resistor(R = 1.0)
@named C = Compliance(C = 1.3)
@named ground = Ground()

circ_eqs = [
    connect(source.out, Rp.in) 
	connect(Rp.out, C.in)
	connect(C.out, Rc.in)
    connect(Rc.out, source.in, ground.g) 
]

## Compose the whole ODAE system
@named _circ_model = ODESystem(circ_eqs, t)

##
@named circ_model = compose(_circ_model,[source, ground, C, Rp, Rc])

## And simplify it
@time circ_sys = structural_simplify(circ_model)
# Examine Equation

# give mean pressures as IC
u0 = [9.0]
@time prob = ODAEProblem(circ_sys, u0, (0.0, 10.0))

@time sol = solve(prob)

plot(sol)
plot(sol, idxs = Rp.q)

## 4 Element Windkessel series ##
@parameters t

begin
	function heaviside(t)
	   0.5 * (sign(t) + 1)
	end
	
	function I_sin(t)
		sin(2pi*(t))^2 * (heaviside(mod(t,1)) - heaviside(mod(t,1)-0.5))
	end
end

plot(I_sin, xlim=[0,2], label="I", xlabel="time", ylabel="I [ml/s]")

@named source = DrivenFlow(Q=100.0, fun=I_sin)
@named L = Inductance(L = 1e-2)
@named Rp = Resistor(R = 0.1)
@named Rc = Resistor(R = 1.0)
@named C = Compliance(C = 1.3)
@named ground = Ground()

circ_eqs = [
    connect(source.out, L.in)
	connect(L.out, Rp.in) 
	connect(Rp.out, C.in)
	connect(C.out, Rc.in)
    connect(Rc.out, source.in, ground.g) 
]

## Compose the whole ODAE system
@named _circ_model = ODESystem(circ_eqs, t)

##
@named circ_model = compose(_circ_model,[source, ground, C, Rp, Rc, L])

## And simplify it
@time circ_sys = structural_simplify(circ_model)
# Examine Equation

# give mean pressures as IC
u0 = [9.0]
@time prob = ODAEProblem(circ_sys, u0, (0.0, 10.0))

@time sol = solve(prob)

plot(sol)
plot(sol, idxs = L.q)
plot(sol, idxs = C.V)
# - Ply for 20 minutes examining different elements etc 

# - Jump to presentation about valves and cardiac chambers 
## Nikolai ODE model 

# Function for the valves
function Valve(R, deltaP)
    q = 0.0
    if (-deltaP) < 0.0 
        q =  deltaP/R
    else
        q = 0.0
    end
    return q

end

# Functions for the Ventricle
function ShiElastance(t, Eₘᵢₙ, Eₘₐₓ, τ, τₑₛ, τₑₚ, Eshift)
    τₑₛ = τₑₛ*τ
    τₑₚ = τₑₚ*τ
    #τ = 4/3(τₑₛ+τₑₚ)
    tᵢ = rem(t + (1 - Eshift) * τ, τ)

    Eₚ = (tᵢ <= τₑₛ) * (1 - cos(tᵢ / τₑₛ * pi)) / 2 +
         (tᵢ > τₑₛ) * (tᵢ <= τₑₚ) * (1 + cos((tᵢ - τₑₛ) / (τₑₚ - τₑₛ) * pi)) / 2 +
         (tᵢ <= τₑₚ) * 0

    E = Eₘᵢₙ + (Eₘₐₓ - Eₘᵢₙ) * Eₚ

    return E
end

function DShiElastance(t, Eₘᵢₙ, Eₘₐₓ, τ, τₑₛ, τₑₚ, Eshift)

    τₑₛ = τₑₛ*τ
    τₑₚ = τₑₚ*τ
    #τ = 4/3(τₑₛ+τₑₚ)
    tᵢ = rem(t + (1 - Eshift) * τ, τ)

    DEₚ = (tᵢ <= τₑₛ) * pi / τₑₛ * sin(tᵢ / τₑₛ * pi) / 2 +
          (tᵢ > τₑₛ) * (tᵢ <= τₑₚ) * pi / (τₑₚ - τₑₛ) * sin((τₑₛ - tᵢ) / (τₑₚ - τₑₛ) * pi) / 2
    (tᵢ <= τₑₚ) * 0
    DE = (Eₘₐₓ - Eₘᵢₙ) * DEₚ

    return DE
end


#Parameters
Eshift = 0.0
Eₘᵢₙ = 0.03
τₑₛ = 0.3
τₑₚ = 0.45 
Eₘₐₓ = 1.5
Rmv = 0.06
τ = 1.0

function NIK!(du, u, p, t)
    pLV, psa, psv, Vlv, Qav, Qmv, Qs = u # Model variables
    τₑₛ, τₑₚ, Rmv, Zao, Rs, Csa, Csv, Eₘₐₓ, Eₘᵢₙ = p # Model parameters
    # the differential equations
    du[1] = (Qmv - Qav) * ShiElastance(t, Eₘᵢₙ, Eₘₐₓ, τ, τₑₛ, τₑₚ, Eshift) + pLV / ShiElastance(t, Eₘᵢₙ, Eₘₐₓ, τ, τₑₛ, τₑₚ, Eshift) * DShiElastance(t, Eₘᵢₙ, Eₘₐₓ, τ, τₑₛ, τₑₚ, Eshift)
    # 1 Left Ventricle
    du[2] = (Qav - Qs ) / Csa #Systemic arteries     
    du[3] = (Qs - Qmv) / Csv # Venous
    du[4] = Qmv - Qav # volume
    du[5]    = Valve(Zao, (pLV - psa)) - Qav # AV 
    du[6]   = Valve(Rmv, (psv - pLV)) - Qmv # MV
    du[7]     = (du[2] - du[3]) / Rs # Systemic flow
    nothing 
end
# Mass Matrix as ODAE problem
M = [1.  0  0  0  0  0  0
     0  1.  0  0  0  0  0
     0  0  1.  0  0  0  0
     0  0  0  1.  0  0  0
     0  0  0  0  0  0  0
     0  0  0  0  0  0  0 
     0  0  0  0  0  0  1. ]
     
Nik_ODE = ODEFunction(NIK!,mass_matrix=M)

# Inintal conditions
u0 = [8.0, 8.0, 8.0, 265.0, 0.0, 0.0, 0.0]
# Model Parameters
p = [0.3, 0.45, 0.06, 0.033, 1.11, 1.13, 11.0, 1.5, 0.03]
# Time span we solve over 
tspan = (0, 20)
# Define ODE problem
prob = ODEProblem(Nik_ODE, u0, tspan, p)

# Solve the problem 
# ODAE problem use Rodas4 solver
# autidiff = false due to undifferentiable valve statements 
@time sol = solve(prob, Rodas4(autodiff = false),  reltol = 1e-10, abstol = 1e-10)
# Plot solution 
plot(sol)

# Nikolai Symbolic model 

Emin = 0.03
τ = 1
τes = 0.3*τ
τep = 0.45*τ
Emax = 1.5
Rmv = 0.006
Zao = 0.033
Rs = 1.11
Rab = 1.2
Cao = 1.13
Csv = 11.0 

@parameters t

# Define elements of the model 
# Ventricle
@named LV = ShiChamber(V₀ = 20.0, p₀ = 1.0, Eₘᵢₙ = Emin, Eₘₐₓ = Emax, τ = τ, τₑₛ = τes, τₑₚ = τep, Eshift=0.0)

#Valves 
@named AVa = ResistorDiode(R = Zao)
@named MV = ResistorDiode(R = Rmv)

#Chambers 
@named AO = Compliance(C = Cao) # Aorta
@named VS = Compliance(C = Csv) # venous System
@named SR = Resistor(R = Rs) # Systemic Resistance


# Connect the equations given the direction of the flow 
circ_eqs = [
    connect(LV.out, AVa.in)
    connect(AVa.out, AO.in)
    connect(AO.out, SR.in) 
    connect(SR.out, VS.in)
    connect(VS.out, MV.in)
    connect(MV.out, LV.in)
]

## Compose the whole ODAE system
@named _circ_model = ODESystem(circ_eqs, t)

# Assign variable names 
@named circ_model = compose(_circ_model,[LV, AVa, MV, AO, SR, VS])
# equations of the model 
equations(circ_model) #unusable in this form 

## Simplfy equations - eliminating variables 
@time circ_sys = structural_simplify(circ_model)
# Print equations 
equations(circ_sys) # Same as what would have given we formed by hand - note LV.V 
observed(circ_sys) # Gives equations for whole system 
# notice now all algebraic variables have been eliminated as we left them in previosly 

# Inital conditions
u0 = [100.0, 90.0, 8.0]
@time prob = ODAEProblem(circ_sys, u0, (0.0, 20.0))

# Vern9 is a high order ode solv er not compatible with algebraic equations 
@time sol = solve(prob, Vern9(), reltol=1e-10, abstol=1e-10)
plot(sol)

# We now plot specifc variables symbolically 

# Aortic Valve flow 
plot(sol, idxs = AVa.q)
# Venous volume 
plot(sol, idxs = VS.V)
# Aortic Pressure 
plot(sol, idxs = AO.p)

# Nik DH elstance 
Emin = 0.03
τ = 1
τ1 = 0.303*τ
τ2 = τ*0.508
n1 = 1.32
n2 = 21.9
nstep = 1000
t = LinRange(0, τ, nstep)
k = 1  / maximum((t ./ τ1).^n1 ./ (1 .+ (t ./ τ1).^n1) .* 1 ./ (1 .+ (t ./ τ2).^n2))
Emax = 1.5
Rmv = 0.006
Zao = 0.033
Rs = 1.11
Rab = 1.2
Cao = 1.13
Csv = 11.0 

@parameters t

# Define elements of the model 
# Ventricle
@named LV = DHChamber(V₀ = 0.0,  Eₘᵢₙ = Emin, Eₘₐₓ = Emax, n₁ = n1, n₂ = n2, τ = τ, τ₁ = τ1, τ₂ = τ2, k = k)
#Valves 
@named AVa = ResistorDiode(R = Zao)
@named MV = ResistorDiode(R = Rmv)

#Chambers 
@named AO = Compliance(C = Cao) # Aorta
@named VS = Compliance(C = Csv) # venous System
@named SR = Resistor(R = Rs) # Systemic Resistance


# Connect the equations given the direction of the flow 
circ_eqs = [
    connect(LV.out, AVa.in)
    connect(AVa.out, AO.in)
    connect(AO.out, SR.in) 
    connect(SR.out, VS.in)
    connect(VS.out, MV.in)
    connect(MV.out, LV.in)
]

## Compose the whole ODAE system
@named _circ_model = ODESystem(circ_eqs, t)

# Assign variable names 
@named circ_model = compose(_circ_model,[LV, AVa, MV, AO, SR, VS])

## Simplfy equations - eliminating variables 
@time circ_sys = structural_simplify(circ_model)
# Print equations 
equations(circ_sys) # Same as what would have given we formed by hand - note LV.p now in pressure 
observed(circ_sys) # Gives equations for whole system 
# notice now all algebraic variables have been eliminated as we left them in previosly 

# Inital conditions
u0 = [5.0, 70.0, 6.0]
@time prob = ODAEProblem(circ_sys, u0, (0.0, 20.0))

# Vern9 is a high order ode solv er not compatible with algebraic equations 
@time sol = solve(prob, Vern9(), reltol=1e-10, abstol=1e-10)
plot(sol)

# We now plot specifc variables symbolically 

# Aortic Valve flow 
plot(sol, idxs = AVa.q)
# Venous volume 
plot(sol, idxs = VS.V)
# Aortic Pressure 
plot(sol, idxs = AO.p)

## More physiological adaptation constant pressure on venous base plate 

Emin = 0.03
τ = 1
τes = 0.3*τ
τep = 0.45*τ
Emax = 1.5
Rmv = 0.006
Zao = 0.033
Rs = 1.11
Rab = 1.2
Cao = 1.13
Csv = 11.0 

@parameters t

# Define elements of the model 
# Ventricle
@named LV = ShiChamber(V₀ = 20.0, p₀ = 1.0, Eₘᵢₙ = Emin, Eₘₐₓ = Emax, τ = τ, τₑₛ = τes, τₑₚ = τep, Eshift=0.0)

#Valves 
@named AVa = ResistorDiode(R = Zao)
@named MV = ResistorDiode(R = Rmv)

#Chambers 
@named AO = Compliance(C = Cao) # Aorta
@named VS = Compliance(C = Csv) # venous System
@named SR = Resistor(R = Rs) # Systemic Resistance
@named BP = ConstantPressure(P = 6.0)
@named ground = Ground()
# Connect the equations given the direction of the flow 
circ_eqs = [
    connect(LV.out, AVa.in)
    connect(AVa.out, AO.in)
    connect(AO.out, SR.in) 
    connect(SR.out, VS.in)
	connect(BP.out, VS.in)
    connect(VS.out, MV.in)
    connect(MV.out, LV.in)
	connect(BP.in, ground.g)
]

## Compose the whole ODAE system
@named _circ_model = ODESystem(circ_eqs, t)

# Assign variable names 
@named circ_model = compose(_circ_model,[LV, AVa, MV, AO, SR, VS, BP, ground])
# equations of the model 
equations(circ_model) #unusable in this form 

## Simplfy equations - eliminating variables 
@time circ_sys = structural_simplify(circ_model)
# Print equations 
equations(circ_sys) # Same as what would have given we formed by hand - note LV.V 
observed(circ_sys) # Gives equations for whole system 
# notice now all algebraic variables have been eliminated as we left them in previosly 

# Inital conditions
u0 = [100.0, 90.0]
@time prob = ODAEProblem(circ_sys, u0, (0.0, 20.0))

# Vern9 is a high order ode solv er not compatible with algebraic equations 
@time sol = solve(prob, Vern9(), reltol=1e-10, abstol=1e-10)
plot(sol)

# We now plot specifc variables symbolically 

# Aortic Valve flow 
plot(sol, idxs = AVa.q)
# Venous volume 
plot(sol, idxs = VS.V)
# Aortic Pressure 
plot(sol, idxs = AO.p)
# Left Ventricular Pressure 
plot(sol, idxs = LV.p)

### Have a play with different elements trying different things then we move onto a more complex model ###



## Combine circulation NIK

function CRC(; name, C=1.0, R=1.0, C1=1.0)
    @named in = Pin()
    @named out = Pin()

    sts = @variables Δp(t) = 0.0 q(t) = 0.0
    ps = []

    # These are the components the subsystem is made of:
    @named C = Compliance(C=C)
    @named R = Resistor(R=R)
    @named C1 = Compliance(C = C1)

    # The equations for the subsystem are created by
    # 'connect'-ing the components

    eqs = [
            Δp ~ out.p - in.p
            q ~ in.q
            connect(in, C.in)
            connect(C.out, R.in)
            connect(R.out, C1.in)
            connect(C1.out, out)
    ]

    # and finaly compose the system
    compose(ODESystem(eqs, t, sts, ps; name=name), in, out, C, R, C1)
end


Emin = 0.03
τ = 1
τes = 0.3*τ
τep = 0.45*τ
Emax = 1.5
Rmv = 0.006
Zao = 0.033
Rs = 1.11
Rab = 1.2
Cao = 1.13
Csv = 11.0 

@parameters t

# Define elements of the model 
# Ventricle
@named LV = ShiChamber(V₀ = 20.0, p₀ = 1.0, Eₘᵢₙ = Emin, Eₘₐₓ = Emax, τ = τ, τₑₛ = τes, τₑₚ = τep, Eshift=0.0)

#Valves 
@named AVa = ResistorDiode(R = Zao)
@named MV = ResistorDiode(R = Rmv)

#Systemic circulation
@named sys_circ = CRC(C = Cao, R = Rs, C1 = Csv)

# Connect the equations given the direction of the flow 
circ_eqs = [
    connect(LV.out, AVa.in)
    connect(AVa.out, sys_circ.in)
    connect(sys_circ.out, MV.in)
    connect(MV.out, LV.in)
]

## Compose the whole ODAE system
@named _circ_model = ODESystem(circ_eqs, t)

# Assign variable names 
@named circ_model = compose(_circ_model,[LV, AVa, MV, sys_circ])
# equations of the model 
equations(circ_model) #unusable in this form 

## Simplfy equations - eliminating variables 
@time circ_sys = structural_simplify(circ_model)
# Print equations 
equations(circ_sys) # Same as what would have given we formed by hand - note LV.V 
observed(circ_sys) # Gives equations for whole system 
# notice now all algebraic variables have been eliminated as we left them in previosly 

# Inital conditions
u0 = [100.0, 90.0, 8.0]
@time prob = ODAEProblem(circ_sys, u0, (0.0, 20.0))

# Vern9 is a high order ode solv er not compatible with algebraic equations 
@time sol = solve(prob, Vern9(), reltol=1e-10, abstol=1e-10)
plot(sol)

# We now plot specifc variables symbolically 

# Aortic Valve flow 
plot(sol, idxs = AVa.q)
# Venous volume 
plot(sol, idxs = sys_circ.C1.V)
# Aortic Pressure 
plot(sol, idxs = sys_circ.C.p)

## Shi systemic model Diode Valves ##

begin 
    τ = 1.0
    v0_lv = 5.0
    p0_lv = 1.0
    Emin_lv = 0.1
    Emax_lv = 2.5
    τes_lv = 0.3
    τed_lv = 0.45
    Eshift_lv = 0.0
    LV_Vt0 = 500
    ### LA Atrium Parameters #### Checked
    v0_la = 4.0
    p0_la = 1.0
    Emin_la = 0.15
    Emax_la = 0.25
    τpwb_la = 0.92
    τpww_la = 0.09
    τes_la = τpww_la/2
    τed_la = τpww_la
    Eshift_la = τpwb_la
    LA_Vt0 = 20
    #####
    Csas = 0.08
    Rsas = 0.003
    Lsas = 6.2e-5
    pt0sas = 100.0
    qt0sas = 0.0
    ## Systemic Artery #### Checked
    Csat = 1.6
    Rsat = 0.05
    Lsat = 0.0017
    pt0sat = 100.0
    qt0sat = 0.0
    ## Systemic Arteriole #### Checked
    Rsar = 0.5
    ## Systemic Capillary #### Checked 
    Rscp = 0.52
    ## Systemic Vein #### Checked
    Csvn = 20.5
    Rsvn = 0.075
    pt0svn = 0.0
    qt0svn = 0.0
    ## Valve Parameters 
    Zao = 0.033
    Rmv = 0.06
end 
  
  
@parameters t
@named LV = ShiChamber(V₀=v0_lv, p₀ = p0_lv, Eₘᵢₙ=Emin_lv, Eₘₐₓ=Emax_lv, τ=τ, τₑₛ=τes_lv, τₑₚ=τed_lv, Eshift=0.0)
@named LA = ShiChamber(V₀=v0_la, p₀ = p0_la, Eₘᵢₙ=Emin_la, Eₘₐₓ=Emax_la, τ=τ, τₑₛ=τpww_la/2, τₑₚ=τpww_la, Eshift=τpwb_la)
# Examine how atrium is different - Show there is an explicit element but can use same one 

@named AV = ResistorDiode(R = Zao)
@named MV = ResistorDiode(R = Rmv)

####### Systemic Loop #######

## Systemic Aortic Sinus ##
@named SAS = CRL(C=Csas, R=Rsas, L=Lsas)


## Systemic Artery ##
@named SAT = CRL(C=Csat, R=Rsat, L=Lsat)

## Systemic Arteriole ##
@named SAR = Resistor(R=Rsar)

## Systemic Capillary ##
@named SCP = Resistor(R=Rscp)

## Systemic Vein ##
@named SVN = CR(R=Rsvn, C=Csvn)


circ_eqs = [
    connect(LV.out, AV.in)
    connect(AV.out, SAS.in)
    connect(SAS.out, SAT.in)
    connect(SAT.out, SAR.in)
    connect(SAR.out, SCP.in)
    connect(SCP.out, SVN.in)
    connect(SVN.out, LA.in)
    connect(LA.out, MV.in)
    connect(MV.out, LV.in)
]

## Compose the whole ODAE system
@named _circ_model = ODESystem(circ_eqs, t)

##
@named circ_model = compose(_circ_model,
    [LV, LA, AV, MV, SAS, SAT, SAR, SCP, SVN])

## And simplify it
@time circ_sys = structural_simplify(circ_model)

u0 =  [LV_Vt0, LA_Vt0, pt0sas, qt0sas , pt0sat, qt0sat, pt0svn]
@time prob = ODAEProblem(circ_sys, u0, (0.0, 20.0))

@time sol = solve(prob, Vern7(), reltol = 1e-9, abstol = 1e-9)
plot(sol, tspan = (15, 16))

plot(sol, idxs = [LV.p, SAS.C.p], tspan = (15, 16))

# Shi with systemic simplified 

begin 
    τ = 1.0
    v0_lv = 5.0
    p0_lv = 1.0
    Emin_lv = 0.1
    Emax_lv = 2.5
    τes_lv = 0.3
    τed_lv = 0.45
    Eshift_lv = 0.0
    LV_Vt0 = 500
    ### LA Atrium Parameters #### Checked
    v0_la = 4.0
    p0_la = 1.0
    Emin_la = 0.15
    Emax_la = 0.25
    τpwb_la = 0.92
    τpww_la = 0.09
    τes_la = τpww_la/2
    τed_la = τpww_la
    Eshift_la = τpwb_la
    LA_Vt0 = 20
    #####
    Csas = 0.08
    Rsas = 0.003
    Lsas = 6.2e-5
    pt0sas = 100.0
    qt0sas = 0.0
    ## Systemic Artery #### Checked
    Csat = 1.6
    Rsat = 0.05
    Lsat = 0.0017
    pt0sat = 100.0
    qt0sat = 0.0
    ## Systemic Arteriole #### Checked
    Rsar = 0.5
    ## Systemic Capillary #### Checked 
    Rscp = 0.52
    ## Systemic Vein #### Checked
    Csvn = 20.5
    Rsvn = 0.075
    pt0svn = 0.0
    qt0svn = 0.0
    ## Valve Parameters 
    Zao = 0.033
    Rmv = 0.06
end 
    
@parameters t
@named LV = ShiChamber(V₀=v0_lv, p₀ = p0_lv, Eₘᵢₙ=Emin_lv, Eₘₐₓ=Emax_lv, τ=τ, τₑₛ=τes_lv, τₑₚ=τed_lv, Eshift=0.0)
@named LA = ShiChamber(V₀=v0_la, p₀ = p0_la, Eₘᵢₙ=Emin_la, Eₘₐₓ=Emax_la, τ=τ, τₑₛ=τpww_la/2, τₑₚ=τpww_la, Eshift=τpwb_la)
# Examine how atrium is different 

@named AV = ResistorDiode(R = Zao)
@named MV = ResistorDiode(R = Rmv)

@named Sys_loop = ShiSystemicLoop(SAS_C=Csas, SAS_R=Rsas, SAS_L=Lsas, SAT_C=Csat, SAT_R=Rsat, SAT_L=Lsat, SAR_R=Rsar, SCP_R=Rscp, SVN_C=Csvn, SVN_R=Rsvn)


function ShiSystemicLoop(; name, SAS_C, SAS_R, SAS_L, SAT_C, SAT_R, SAT_L, SAR_R, SCP_R, SVN_C, SVN_R)
    @named in = Pin()
    @named out = Pin()

    sts = @variables Δp(t) = 0.0 q(t) = 0.0
    ps = []
    # No parameters in this function

    # These are the components the subsystem is made of:
    ## Systemic Aortic Sinus ##
    @named SAS = CRL(C=SAS_C, R=SAS_R, L=SAS_L)
    ## Systemic Artery ##
    @named SAT = CRL(C=SAT_C, R=SAT_R, L=SAT_L)
    ## Systemic Arteriole ##
    @named SAR = Resistor(R=SAR_R)
    ## Systemic Capillary ##
    @named SCP = Resistor(R=SCP_R)
    ## Systemic Vein ##
    @named SVN = CR(C=SVN_C, R=SVN_R)

    # The equations for the subsystem are created by
    # 'connect'-ing the components
    eqs = [
            Δp ~ out.p - in.p
            q ~ in.q
            connect(in, SAS.in)
            connect(SAS.out, SAT.in)
            connect(SAT.out, SAR.in)
            connect(SAR.out, SCP.in)
            connect(SCP.out, SVN.in)
            connect(SVN.out, out)
    ]

    # and finaly compose the system
    compose(ODESystem(eqs, t, sts, ps; name=name), in, out, SAS, SAT, SAR, SCP, SVN)
end

circ_eqs = [
    connect(LV.out, AV.in)
    connect(AV.out, Sys_loop.in)
    connect(Sys_loop.out, LA.in)
    connect(LA.out, MV.in)
    connect(MV.out, LV.in)
]

## Compose the whole ODAE system
@named _circ_model = ODESystem(circ_eqs, t)

##
@named circ_model = compose(_circ_model,[LV, LA, AV, MV, Sys_loop ])

## And simplify it
@time circ_sys = structural_simplify(circ_model)

u0 =  [LV_Vt0, LA_Vt0, pt0sas, qt0sas , pt0sat, qt0sat, pt0svn]

@time prob = ODAEProblem(circ_sys, u0, (0.0, 20.0))

@time sol = solve(prob, Vern7(), reltol = 1e-9, abstol = 1e-9)
plot(sol, tspan = (15, 16))

plot(sol, idxs = [LV.p, Sys_loop.SAS.C.p], tspan = (15, 16))


### Shi Valve with CallBacks  - Systemic Loop 
function ShiValve(; name, CQ, Kp, Kf, Kb, Kv, θmax, θmin)
    @named oneport = OnePort()
    @unpack Δp, q = oneport
    ps = @parameters CQ = CQ Kp = Kp Kf = Kf Kb = Kb Kv = Kv θmax = θmax θmin = θmin
    sts = @variables θ(t) = 0.0 ω(t) = 0.0 AR(t) = 0.0 Fp(t) = 0.0 Ff(t) = 0.0 Fb(t) = 0.0 Fv(t) = 0.0 F(t) = 0.0
    D = Differential(t)
    limits = [
            [(θ ~ θmax)] => [ω ~ 0]
            [(θ ~ θmin)] => [ω ~ 0]
    ]

    # make θmax the real opening angle and define a θmaxopen for a healthy valve
    # that means we can use θmax as a stenosis parameter
    θmaxopen = 75 * pi / 180

    eqs = [
            # Forces/Moments
            Fp ~ Kp * -Δp * cos(θ)                 # pressure
            Ff ~ -Kf * ω                           # friction
            Fb ~ Kb * q * cos(θ)                       # Fluid Velocity
            Fv ~ -Kv * q * (q > 0) * sin(2θ)       # vortex behind leaflets
            F ~ Fp + Ff + Fb + Fv                  # total force/moment on leaflets
            #ODEs
            D(θ) ~ ω
            D(ω) ~ F * ((θ < θmax) * (F > 0) + (θ > θmin) * (F < 0))
            # Opening ratio
            #AR ~ ((1 - cos(θ))^2) / ((1 - cos(θmax))^2)
            AR ~ ((1 - cos(θ))^2) / ((1 - cos(θmaxopen))^2)
            # Flow equation
            q ~ -sign(Δp) * CQ * AR * sqrt(abs(Δp))
    ]

    # include the `continuous_events` definition `limits` in the ODE system
    # this is the MTK equivalent to callbacks
    extend(ODESystem(eqs, t, sts, ps; name=name, continuous_events=limits), oneport)
end


begin 
    τ = 1.0
    v0_lv = 5.0
    p0_lv = 1.0
    Emin_lv = 0.1
    Emax_lv = 2.5
    τes_lv = 0.3
    τed_lv = 0.45
    Eshift_lv = 0.0
    LV_Vt0 = 500
    ### LA Atrium Parameters #### Checked
    v0_la = 4.0
    p0_la = 1.0
    Emin_la = 0.15
    Emax_la = 0.25
    τpwb_la = 0.92
    τpww_la = 0.09
    τes_la = τpww_la/2
    τed_la = τpww_la
    Eshift_la = τpwb_la
    LA_Vt0 = 20
    #####
    Csas = 0.08
    Rsas = 0.003
    Lsas = 6.2e-5
    pt0sas = 100.0
    qt0sas = 0.0
    ## Systemic Artery #### Checked
    Csat = 1.6
    Rsat = 0.05
    Lsat = 0.0017
    pt0sat = 100.0
    qt0sat = 0.0
    ## Systemic Arteriole #### Checked
    Rsar = 0.5
    ## Systemic Capillary #### Checked 
    Rscp = 0.52
    ## Systemic Vein #### Checked
    Csvn = 20.5
    Rsvn = 0.075
    pt0svn = 0.0
    qt0svn = 0.0
    ## Valve Parameters 
    Zao = 0.033
    Rmv = 0.06
end 

θmax = 75.0*pi/180
θmin = 5.0*pi/180
Kp_av = 5500.0 
Kf_av = 50.0
Kb_av = 2.0
Kv_av = 7.0
Kf_mv = 50.0
Kp_mv = 5500.0
Kb_mv = 2.0
Kv_mv = 3.5
CQ_AV = 350.0
CQ_MV = 400.0

@parameters t
@named LV = ShiChamber(V₀=v0_lv, p₀ = p0_lv, Eₘᵢₙ=Emin_lv, Eₘₐₓ=Emax_lv, τ=τ, τₑₛ=τes_lv, τₑₚ=τed_lv, Eshift=0.0)
@named LA = ShiChamber(V₀=v0_la, p₀ = p0_la, Eₘᵢₙ=Emin_la, Eₘₐₓ=Emax_la, τ=τ, τₑₛ=τpww_la/2, τₑₚ=τpww_la, Eshift=τpwb_la)

@named AV = ShiValve(CQ = CQ_AV, Kp = Kp_av, Kf = Kf_av, Kb = Kb_av, Kv = Kv_av, θmax = θmax, θmin = θmin)
@named MV = ShiValve(CQ = CQ_MV, Kp = Kp_mv, Kf = Kf_mv, Kb = Kb_mv, Kv = Kv_mv, θmax = θmax, θmin = θmin)

@named Sys_loop = ShiSystemicLoop(SAS_C=Csas, SAS_R=Rsas, SAS_L=Lsas, SAT_C=Csat, SAT_R=Rsat, SAT_L=Lsat, SAR_R=Rsar, SCP_R=Rscp, SVN_C=Csvn, SVN_R=Rsvn)

circ_eqs = [
    connect(LV.out, AV.in)
    connect(AV.out, Sys_loop.in)
    connect(Sys_loop.out, LA.in)
    connect(LA.out, MV.in)
    connect(MV.out, LV.in)
]

## Compose the whole ODAE system
@named _circ_model = ODESystem(circ_eqs, t)

##
@named circ_model = compose(_circ_model,[LV, LA, AV, MV, Sys_loop ])

## And simplify it
@time circ_sys = structural_simplify(circ_model)

u0 =  [LV_Vt0, LA_Vt0, 0.0, 0.0, 0.0, 0.0, pt0sas, qt0sas , pt0sat, qt0sat, pt0svn]

@time prob = ODAEProblem(circ_sys, u0, (0.0, 20.0))

@time sol = solve(prob, Vern7(), reltol = 1e-9, abstol = 1e-9)

plot(sol, idxs = [LV.p, Sys_loop.SAS.C.p], tspan = (15, 16))

# Don't have to write new ode systems can focus on a single element 


## Add pulmonary Loop
begin 
    τ = 1.0
    Eshift=0.0
    Ev=Inf
    #### LV chamber parameters #### Checked
    v0_lv = 5.0
    p0_lv = 1.0
    Emin_lv = 0.1
    Emax_lv = 2.5
    τes_lv = 0.3
    τed_lv = 0.45
    Eshift_lv = 0.0
    #### RV Chamber parameters #### Checked
    v0_rv = 10.0
    p0_rv = 1.0
    Emin_rv = 0.1
    Emax_rv = 1.15
    τes_rv = 0.3
    τed_rv = 0.45
    Eshift_rv = 0.0
    ### LA Atrium Parameters #### Checked
    v0_la = 4.0
    p0_la = 1.0
    Emin_la = 0.15
    Emax_la = 0.25
    τpwb_la = 0.92
    τpww_la = 0.09
    τes_la = τpww_la/2
    τed_la = τpww_la
    Eshift_la = τpwb_la
    ### RA Atrium parameters #### Checked
    v0_ra = 4.0
    p0_ra = 1.0
    Emin_ra = 0.15
    Emax_ra = 0.25
    τpwb_ra = 0.92
    τpww_ra = 0.09
    τes_ra = τpww_ra/2
    τed_ra = τpww_ra
    Eshift_ra = τpwb_ra
    #### Valve parameters #### Checked
    CQ_AV = 350.0
    CQ_MV = 400.0
    CQ_TV = 400.0
    CQ_PV = 350.0
    ## Systemic Aortic Sinus #### Checked
    Csas = 0.08
    Rsas = 0.003
    Lsas = 6.2e-5
    pt0sas = 100.0
    qt0sas = 0.0
    ## Systemic Artery #### Checked
    Csat = 1.6
    Rsat = 0.05
    Lsat = 0.0017
    pt0sat = 100.0
    qt0sat = 0.0
    ## Systemic Arteriole #### Checked
    Rsar = 0.5
    ## Systemic Capillary #### Checked 
    Rscp = 0.52
    ## Systemic Vein #### Checked
    Csvn = 20.5
    Rsvn = 0.075
    pt0svn = 0.0
    qt0svn = 0.0
    ## Pulmonary Aortic Sinus #### Checked
    Cpas = 0.18
    Rpas = 0.002
    Lpas = 5.2e-5
    pt0pas = 30.0
    qt0pas = 0.0
    ## Pulmonary Artery #### Checked
    Cpat = 3.8
    Rpat = 0.01
    Lpat = 0.0017
    pt0pat = 30.0
    qt0pat = 0.0
    ## Pulmonary Arteriole #### Checked
    Rpar = 0.05
    ## Pulmonary Capillary #### Checked
    Rpcp = 0.25
    ## Pulmonary Vein #### Checked
    Cpvn = 20.5
    Rpvn =  0.0006 #   0.0006        # this was 0.006 originally and in the paper, seems to be wrong in the paper!
    # CHANGED THIS IN THE CELLML MODEL AS WELL TO MATCH THE PAPER!!!!!
    pt0pvn = 0.0
    qt0pvn = 0.0
    ## KG diaphragm ## Not in cellML model
    # left heart #
    Kst_la = 2.5
    Kst_lv = 20.0
    Kf_sav = 0.0004
    Ke_sav = 9000.0
    M_sav = 0.0004
    A_sav = 0.00047
    # right heart # 
    Kst_ra = 2.5
    Kst_rv = 20.0
    Kf_pav = 0.0004
    Ke_pav = 9000.0
    M_pav = 0.0004
    A_pav = 0.00047
    #
    #### Diff valve params #### not in cellML model
    Kp_av = 5500.0 # *  57.29578 # Shi Paper has values in radians!
    Kf_av = 50.0
    Kf_mv = 50.0
    Kp_mv = 5500.0 # *  57.29578 
    Kf_tv = 50.0
    Kp_tv = 5500.0 # *  57.29578
    Kf_pv = 50.0
    Kp_pv = 5500.0 #*  57.29578
    Kb_av = 2.0
    Kv_av = 7.0
    Kb_mv = 2.0
    Kv_mv = 3.5
    Kb_tv = 2.0
    Kv_tv = 3.5
    Kb_pv = 2.0
    Kv_pv = 3.5
    θmax_av = 75.0 * pi / 180
    θmax_mv = 75.0 * pi / 180
    θmin_av = 5.0 * pi / 180
    θmin_mv = 5.0 * pi / 180
    θmax_pv = 75.0 * pi / 180
    θmax_tv = 75.0 * pi / 180
    θmin_pv = 5.0 * pi / 180
    θmin_tv = 5.0 * pi / 180

    ## pressure force and frictional force is the same for all 4 valves 

    # Initial conditions #### Checked against cellML model

    LV_Vt0 = 500
    RV_Vt0 = 500
    LA_Vt0 = 20
    RA_Vt0 = 20
end 

@named AV = ShiValve(CQ=CQ_AV, Kp=Kp_av, Kf=Kf_pv, Kb=Kb_av, Kv=Kv_av, θmax=θmax, θmin=θmin)
@named MV = ShiValve(CQ=CQ_MV, Kp=Kp_mv, Kf=Kf_mv, Kb=Kb_mv, Kv=Kv_mv, θmax=θmax, θmin=θmin)
@named TV = ShiValve(CQ=CQ_TV, Kp=Kp_tv, Kf=Kf_tv, Kb=Kb_tv, Kv=Kv_tv, θmax=θmax, θmin=θmin)
@named PV = ShiValve(CQ=CQ_PV, Kp=Kp_pv, Kf=Kf_pv, Kb=Kb_pv, Kv=Kv_pv, θmax=θmax, θmin=θmin)

@named LV = ShiChamber(V₀=v0_lv, p₀ = p0_lv, Eₘᵢₙ=Emin_lv, Eₘₐₓ=Emax_lv, τ=τ, τₑₛ=τes_lv, τₑₚ=τed_lv, Eshift=0.0)
@named LA = ShiChamber(V₀=v0_la, p₀ = p0_la, Eₘᵢₙ=Emin_la, Eₘₐₓ=Emax_la, τ=τ, τₑₛ=τpww_la/2, τₑₚ=τpww_la, Eshift=τpwb_la)
@named RV = ShiChamber(V₀=v0_rv, p₀ = p0_rv, Eₘᵢₙ=Emin_rv, Eₘₐₓ=Emax_rv, τ=τ, τₑₛ=τes_rv, τₑₚ=τed_rv, Eshift=0.0)
@named RA = ShiChamber(V₀=v0_ra, p₀ = p0_ra, Eₘᵢₙ=Emin_ra, Eₘₐₓ=Emax_ra, τ=τ, τₑₛ=τpww_ra/2, τₑₚ =τpww_ra, Eshift=τpwb_ra)

@named Sys_loop = ShiSystemicLoop(SAS_C=Csas, SAS_R=Rsas, SAS_L=Lsas, SAT_C=Csat, SAT_R=Rsat, SAT_L=Lsat, SAR_R=Rsar, SCP_R=Rscp, SVN_C=Csvn, SVN_R=Rsvn)
@named Pul_loop = ShiPulmonaryLoop(PAS_C=Cpas, PAS_R=Rpas, PAS_L=Lpas, PAT_C=Cpat, PAT_R=Rpat, PAT_L=Lpat, PAR_R=Rpar, PCP_R=Rpcp, PVN_C=Cpvn, PVN_R=Rpvn)

circ_eqs = [
    connect(LV.out, AV.in)
    connect(AV.out, Sys_loop.in)
    connect(Sys_loop.out, RA.in)
    connect(RA.out, TV.in)
    connect(TV.out, RV.in)
    connect(RV.out, PV.in)
    connect(PV.out, Pul_loop.in)
    connect(Pul_loop.out, LA.in)
    connect(LA.out, MV.in)
    connect(MV.out, LV.in)
]

## Compose the whole ODAE system
@named _circ_model = ODESystem(circ_eqs, t)

##
@named circ_model = compose(_circ_model,
    [LV, RV, LA, RA, AV, MV, PV, TV, Sys_loop, Pul_loop])

## And simplify it
@time circ_sys = structural_simplify(circ_model)


u0 =  [LV_Vt0, RV_Vt0, LA_Vt0, RA_Vt0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, pt0sas, qt0sas , pt0sat, qt0sat, pt0svn, pt0pas, qt0pas, pt0pat, qt0pat, pt0pvn]
@time prob = ODAEProblem(circ_sys, u0, (0.0, 20.0))
##
@time sol = solve(prob)
plot(sol, idxs = [LV.p, Sys_loop.SAS.C.p], tspan = (15, 16))

## Show how shi heart is created 
begin 
    τ = 1.0
    Eshift=0.0
    Ev=Inf
    #### LV chamber parameters #### Checked
    v0_lv = 5.0
    p0_lv = 1.0
    Emin_lv = 0.1
    Emax_lv = 2.5
    τes_lv = 0.3
    τed_lv = 0.45
    Eshift_lv = 0.0
    #### RV Chamber parameters #### Checked
    v0_rv = 10.0
    p0_rv = 1.0
    Emin_rv = 0.1
    Emax_rv = 1.15
    τes_rv = 0.3
    τed_rv = 0.45
    Eshift_rv = 0.0
    ### LA Atrium Parameters #### Checked
    v0_la = 4.0
    p0_la = 1.0
    Emin_la = 0.15
    Emax_la = 0.25
    τpwb_la = 0.92
    τpww_la = 0.09
    τes_la = τpww_la/2
    τed_la = τpww_la
    Eshift_la = τpwb_la
    ### RA Atrium parameters #### Checked
    v0_ra = 4.0
    p0_ra = 1.0
    Emin_ra = 0.15
    Emax_ra = 0.25
    τpwb_ra = 0.92
    τpww_ra = 0.09
    τes_ra = τpww_ra/2
    τed_ra = τpww_ra
    Eshift_ra = τpwb_ra
    #### Valve parameters #### Checked
    CQ_AV = 350.0
    CQ_MV = 400.0
    CQ_TV = 400.0
    CQ_PV = 350.0
    ## Systemic Aortic Sinus #### Checked
    Csas = 0.08
    Rsas = 0.003
    Lsas = 6.2e-5
    pt0sas = 100.0
    qt0sas = 0.0
    ## Systemic Artery #### Checked
    Csat = 1.6
    Rsat = 0.05
    Lsat = 0.0017
    pt0sat = 100.0
    qt0sat = 0.0
    ## Systemic Arteriole #### Checked
    Rsar = 0.5
    ## Systemic Capillary #### Checked 
    Rscp = 0.52
    ## Systemic Vein #### Checked
    Csvn = 20.5
    Rsvn = 0.075
    pt0svn = 0.0
    qt0svn = 0.0
    ## Pulmonary Aortic Sinus #### Checked
    Cpas = 0.18
    Rpas = 0.002
    Lpas = 5.2e-5
    pt0pas = 30.0
    qt0pas = 0.0
    ## Pulmonary Artery #### Checked
    Cpat = 3.8
    Rpat = 0.01
    Lpat = 0.0017
    pt0pat = 30.0
    qt0pat = 0.0
    ## Pulmonary Arteriole #### Checked
    Rpar = 0.05
    ## Pulmonary Capillary #### Checked
    Rpcp = 0.25
    ## Pulmonary Vein #### Checked
    Cpvn = 20.5
    Rpvn =  0.0006 #   0.0006        # this was 0.006 originally and in the paper, seems to be wrong in the paper!
    # CHANGED THIS IN THE CELLML MODEL AS WELL TO MATCH THE PAPER!!!!!
    pt0pvn = 0.0
    qt0pvn = 0.0
    ## KG diaphragm ## Not in cellML model
    # left heart #
    Kst_la = 2.5
    Kst_lv = 20.0
    Kf_sav = 0.0004
    Ke_sav = 9000.0
    M_sav = 0.0004
    A_sav = 0.00047
    # right heart # 
    Kst_ra = 2.5
    Kst_rv = 20.0
    Kf_pav = 0.0004
    Ke_pav = 9000.0
    M_pav = 0.0004
    A_pav = 0.00047
    #
    #### Diff valve params #### not in cellML model
    Kp_av = 5500.0 # *  57.29578 # Shi Paper has values in radians!
    Kf_av = 50.0
    Kf_mv = 50.0
    Kp_mv = 5500.0 # *  57.29578 
    Kf_tv = 50.0
    Kp_tv = 5500.0 # *  57.29578
    Kf_pv = 50.0
    Kp_pv = 5500.0 #*  57.29578
    Kb_av = 2.0
    Kv_av = 7.0
    Kb_mv = 2.0
    Kv_mv = 3.5
    Kb_tv = 2.0
    Kv_tv = 3.5
    Kb_pv = 2.0
    Kv_pv = 3.5
    θmax_av = 75.0 * pi / 180
    θmax_mv = 75.0 * pi / 180
    θmin_av = 5.0 * pi / 180
    θmin_mv = 5.0 * pi / 180
    θmax_pv = 75.0 * pi / 180
    θmax_tv = 75.0 * pi / 180
    θmin_pv = 5.0 * pi / 180
    θmin_tv = 5.0 * pi / 180

    ## pressure force and frictional force is the same for all 4 valves 

    # Initial conditions #### Checked against cellML model

    LV_Vt0 = 500
    RV_Vt0 = 500
    LA_Vt0 = 20
    RA_Vt0 = 20
end 

function ShiHeart(; name, τ, LV_V₀, LV_p0, LV_Emin, LV_Emax, LV_τes, LV_τed, LV_Eshift, RV_V₀, RV_p0, RV_Emin, RV_Emax, RV_τes, RV_τed, RV_Eshift, LA_V₀, LA_p0, LA_Emin, LA_Emax, LA_τes, LA_τed, LA_Eshift, RA_V₀, RA_p0, RA_Emin, RA_Emax, RA_τes, RA_τed, RA_Eshift, AV_CQ, AV_Kp, AV_Kf, AV_Kb, AV_Kv, AV_θmax, AV_θmin, PV_CQ, PV_Kp, PV_Kf, PV_Kb, PV_Kv, PV_θmax, PV_θmin, MV_CQ, MV_Kp, MV_Kf, MV_Kb, MV_Kv, MV_θmax, MV_θmin, TV_CQ, TV_Kp, TV_Kf, TV_Kb, TV_Kv, TV_θmax, TV_θmin)
    @named LHin = Pin()
    @named LHout = Pin()
    @named RHin = Pin()
    @named RHout = Pin()
    @named in = Pin()
    @named out = Pin()

    sts = @variables Δp(t) = 0.0 q(t) = 0.0
    # sts = []
    # ps = @parameters Rc=Rc Rp=Rp C=C
    # No parameters in this function
    # Parameters are inherited from subcomponents
    ps = []

    # These are the components the subsystem is made of:
    # Ventricles and atria
    @named LV = ShiChamber(V₀=LV_V₀, p₀=LV_p0, Eₘᵢₙ=LV_Emin, Eₘₐₓ=LV_Emax, τ=τ, τₑₛ=LV_τes, τₑₚ=LV_τed, Eshift=LV_Eshift)
    @named LA = ShiChamber(V₀=LA_V₀, p₀=LA_p0, Eₘᵢₙ=LA_Emin, Eₘₐₓ=LA_Emax, τ=τ, τₑₛ=LA_τes, τₑₚ=LA_τed, Eshift=LA_Eshift)
    @named RV = ShiChamber(V₀=RV_V₀, p₀=RV_p0, Eₘᵢₙ=RV_Emin, Eₘₐₓ=RV_Emax, τ=τ, τₑₛ=RV_τes, τₑₚ=RV_τed, Eshift=RV_Eshift)
    @named RA = ShiChamber(V₀=RA_V₀, p₀=RA_p0, Eₘᵢₙ=RA_Emin, Eₘₐₓ=RA_Emax, τ=τ, τₑₛ=RA_τes, τₑₚ=RA_τed, Eshift=RA_Eshift)
    # Valves
    @named AV = ShiValve(CQ=AV_CQ, Kp=AV_Kp, Kf=AV_Kf, Kb=AV_Kb, Kv=AV_Kv, θmax=AV_θmax, θmin=AV_θmin)
    @named MV = ShiValve(CQ=MV_CQ, Kp=MV_Kp, Kf=MV_Kf, Kb=MV_Kb, Kv=MV_Kv, θmax=MV_θmax, θmin=MV_θmin)
    @named TV = ShiValve(CQ=TV_CQ, Kp=TV_Kp, Kf=TV_Kf, Kb=TV_Kb, Kv=TV_Kv, θmax=TV_θmax, θmin=TV_θmin)
    @named PV = ShiValve(CQ=PV_CQ, Kp=PV_Kp, Kf=PV_Kf, Kb=PV_Kb, Kv=PV_Kv, θmax=PV_θmax, θmin=PV_θmin)
    # The equations for the subsystem are created by
    # 'connect'-ing the components

    eqs = [
        Δp ~ out.p - in.p
        q ~ in.q
        connect(LHin, LA.in) 
        connect(LA.out, MV.in)
        connect(MV.out, LV.in)
        connect(LV.out, AV.in)
        connect(AV.out, LHout)
        connect(RHin, RA.in) 
        connect(RA.out, TV.in)
        connect(TV.out, RV.in)
        connect(RV.out, PV.in)
        connect(PV.out, RHout)
        ]
    # and finaly compose the system
    compose(ODESystem(eqs, t, sts, ps; name=name), in, out, LHin, LHout, RHin, RHout, LV, RV, LA, RA, AV, MV, TV, PV)
end


@named Heart = ShiHeart(τ = τ, LV_V₀ = v0_lv, LV_p0 = p0_lv, LV_Emin = Emin_lv, LV_Emax = Emax_lv, LV_τes = τes_lv, LV_τed = τed_lv, LV_Eshift = Eshift_lv, RV_V₀ = v0_rv, RV_p0 = p0_rv, RV_Emin = Emin_rv, RV_Emax = Emax_rv, RV_τes = τes_rv, RV_τed = τed_rv, RV_Eshift = Eshift_rv, LA_V₀ = v0_la, LA_p0 = p0_la, LA_Emin = Emin_la, LA_Emax = Emax_la, LA_τes = τes_la, LA_τed = τed_la, LA_Eshift = Eshift_la, RA_V₀ = v0_ra, RA_p0 = p0_ra, RA_Emin = Emin_ra, RA_Emax = Emax_ra, RA_τes = τes_ra, RA_τed = τed_ra, RA_Eshift = Eshift_ra, AV_CQ = CQ_AV, AV_Kp = Kp_av, AV_Kf = Kf_av, AV_Kb = Kb_av, AV_Kv = Kv_av, AV_θmax = θmax_av, AV_θmin = θmin_av, PV_CQ = CQ_PV, PV_Kp = Kp_pv, PV_Kf = Kf_pv, PV_Kb = Kb_pv, PV_Kv = Kv_pv, PV_θmax = θmax_pv, PV_θmin = θmin_pv, MV_CQ = CQ_MV, MV_Kp = Kp_mv, MV_Kf = Kf_mv, MV_Kb = Kb_mv, MV_Kv = Kv_mv, MV_θmax = θmax_mv, MV_θmin = θmin_mv, TV_CQ = CQ_TV, TV_Kp = Kp_tv, TV_Kf = Kf_tv, TV_Kb = Kb_tv, TV_Kv = Kv_tv, TV_θmax = θmax_tv, TV_θmin = θmin_tv)

@named Sys_loop = ShiSystemicLoop(SAS_C=Csas, SAS_R=Rsas, SAS_L=Lsas, SAT_C=Csat, SAT_R=Rsat, SAT_L=Lsat, SAR_R=Rsar, SCP_R=Rscp, SVN_C=Csvn, SVN_R=Rsvn)
@named Pul_loop = ShiPulmonaryLoop(PAS_C=Cpas, PAS_R=Rpas, PAS_L=Lpas, PAT_C=Cpat, PAT_R=Rpat, PAT_L=Lpat, PAR_R=Rpar, PCP_R=Rpcp, PVN_C=Cpvn, PVN_R=Rpvn)

circ_eqs = [
    connect(Heart.LHout, Sys_loop.in)
    connect(Sys_loop.out, Heart.RHin)
    connect(Heart.RHout, Pul_loop.in)
    connect(Pul_loop.out, Heart.LHin)
]


## Compose the whole ODAE system
@named _circ_model = ODESystem(circ_eqs, t)

##
@named circ_model = compose(_circ_model, [Heart, Sys_loop, Pul_loop])

## And simplify it
@time circ_sys = structural_simplify(circ_model)

## Setup ODE
# Initial Conditions for Shi Valve
u0 = [LV_Vt0, RV_Vt0, LA_Vt0, RA_Vt0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, pt0sas, qt0sas , pt0sat, qt0sat, pt0svn, pt0pas, qt0pas, pt0pat, qt0pat, pt0pvn]

@time prob = ODAEProblem(circ_sys, u0, (0.0, 20.0))

@time sol = solve(prob, Vern7(), reltol=1e-10, abstol=1e-10)

plot(sol, idxs=[Heart.AV.q, Heart.MV.q, Heart.PV.q, Heart.TV.q], tspan=(16 * τ, 18 * τ))

plot(sol, idxs=[Heart.LV.p, Heart.RV.p], tspan=(16 * τ, 18 * τ))
plot(sol, idxs=[Heart.LA.p, Heart.RA.p], tspan=(16 * τ, 18 * τ))

plot(sol, idxs=[Sys_loop.SAS.C.p], tspan=(16 * τ, 18 * τ))

plot(sol, idxs=[Heart.LV.V, Heart.RV.V], tspan=(16 * τ, 18 * τ))
plot(sol, idxs=[Heart.LA.V, Heart.RA.V], tspan=(16 * τ, 18 * τ))

## Remake a single parameter 
parameters(circ_sys)
prob.p[2]

p_ = [
    Heart.LV.p₀ => 5.0
]

prob_ = remake(prob, p = p_)

prob_.p[2]

@time sol_ = solve(prob_, Vern7(), reltol=1e-10, abstol=1e-10)


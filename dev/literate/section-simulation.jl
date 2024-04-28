
using Aeroelasticity, DifferentialEquations, LinearAlgebra

## define non-dimensional parameters
V = range(1e-6, 3.1, length=1000) # = U/(b*ωθ) (reduced velocity)
a = -0.8725 # reference point normalized location
e = 0.1275 # center of mass normalized location
μ = 229.18 # = m/(ρ*pi*b^2) (mass ratio)
r2 = 2.91548 # = Iθ/(m*b^2) (radius of gyration about P)
σ = 0.26116 # = ωh/ωθ (natural frequency ratio)
xθ = e - a # distance from center of mass to reference point
a0 = 2*pi # lift curve slope
α0 = 0 # zero lift angle
cd0 = 0 # drag coefficient
cm0 = 0 # moment coefficient

## choose dimensional parameters
b = 0.115 # semichord
ρ = 1.225 # air density
ωθ = 407.193 # pitch natural frequency
c = 343 # air speed of sound

## calculate dimensionalized parameters
U = V*b*ωθ # freestrean velocity
m = μ*ρ*pi*b^2 # mass
Sθ = m*xθ*b # mass imbalance
Iθ = r2*m*b^2 # inertia
ωh = σ*ωθ # plunge natural frequency
kh = m*ωh^2 # plunge spring constant
kθ = Iθ*ωθ^2 # pitch spring constant

## reduced velocity
V = 65 # = U/(b*ωθ)

## dimensionalized velocity
U = V*b*ωθ

## define submodels
submodels = (Peters{6}(), Section())

## define parameter vector
p = [a, b, a0, α0, cd0, cm0, kh, kθ, m, Sθ, Iθ, U, ρ, c]

## define coupled model
model = CoupledModel(submodels, p)

## construct function
f = ODEFunction(model, p)

## initial states
λ = zeros(6)
h = 1.6
theta = 0.8
hdot = 0.32
thetadot = 0.16
x0 = vcat(λ, h, theta, hdot, thetadot)

## simulate for 100 seconds
tspan = (0.0, 100.0)

## construct problem
prob = ODEProblem(f, x0, tspan, p)

## solve problem
sol = solve(prob, Rodas4())


using Plots
pyplot()

plot(sol,
    idxs = [7,8,9,10],
    xlabel = "t",
    ylabel = permutedims([
        "\$h\$",
        "\$\\theta\$",
        "\$\\dot{h}\$",
        "\$\\dot{\\theta}\$",
        ]),
    label = "",
    layout = (4, 1),
    size = (600,1200)
    )

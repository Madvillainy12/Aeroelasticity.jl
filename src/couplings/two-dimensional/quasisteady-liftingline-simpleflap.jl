"""
    couple_models(aero::QuasiSteady, stru::LiftingLineSection, flap::SimpleFlap,
        ctrl::LiftingLineSectionControl)

Create an aerostructural model using a quasi-steady aerodynamics model, a
lifting line aerodynamic model, and a linear, steady-state control surface model.
The existence of this coupling allows the [`QuasiSteady`](@ref) and [`SimpleFlap`](@ref)
to be used with [`LiftingLine`](@ref) and [`LiftingLineFlaps`](@ref). This model
introduces the freestream air density ``\\rho`` as an additional parameter.
"""
function couple_models(aero::QuasiSteady, stru::LiftingLineSection, flap::SimpleFlap,
    ctrl::LiftingLineSectionControl)

    return (aero, stru, flap, ctrl)
end

# --- Traits --- #

# steady
function number_of_additional_parameters(::Type{QuasiSteady{0}}, ::Type{LiftingLineSection},
    ::Type{SimpleFlap}, ::Type{LiftingLineSectionControl})

    return 1
end

function coupling_inplaceness(::Type{QuasiSteady{0}}, ::Type{LiftingLineSection},
    ::Type{SimpleFlap}, ::Type{LiftingLineSectionControl})

    return OutOfPlace()
end

function coupling_rate_jacobian_type(::Type{QuasiSteady{0}}, ::Type{LiftingLineSection},
    ::Type{SimpleFlap}, ::Type{LiftingLineSectionControl})

    return Zeros()
end

function coupling_state_jacobian_type(::Type{QuasiSteady{0}}, ::Type{LiftingLineSection},
    ::Type{SimpleFlap}, ::Type{LiftingLineSectionControl})

    return Nonlinear()
end

function coupling_parameter_jacobian_type(::Type{QuasiSteady{0}}, ::Type{LiftingLineSection},
    ::Type{SimpleFlap}, ::Type{LiftingLineSectionControl})

    return Nonlinear()
end

function coupling_time_gradient_type(::Type{QuasiSteady{0}}, ::Type{LiftingLineSection},
    ::Type{SimpleFlap}, ::Type{LiftingLineSectionControl})

    return Zeros()
end

# quasisteady without acceleration terms
function number_of_additional_parameters(::Type{QuasiSteady{1}}, ::Type{LiftingLineSection},
    ::Type{SimpleFlap}, ::Type{LiftingLineSectionControl})

    return 1
end

function coupling_inplaceness(::Type{QuasiSteady{1}}, ::Type{LiftingLineSection},
    ::Type{SimpleFlap}, ::Type{LiftingLineSectionControl})

    return OutOfPlace()
end

function coupling_rate_jacobian_type(::Type{QuasiSteady{1}}, ::Type{LiftingLineSection},
    ::Type{SimpleFlap}, ::Type{LiftingLineSectionControl})

    return Zeros()
end

function coupling_state_jacobian_type(::Type{QuasiSteady{1}}, ::Type{LiftingLineSection},
    ::Type{SimpleFlap}, ::Type{LiftingLineSectionControl})

    return Nonlinear()
end

function coupling_parameter_jacobian_type(::Type{QuasiSteady{1}}, ::Type{LiftingLineSection},
    ::Type{SimpleFlap}, ::Type{LiftingLineSectionControl})

    return Nonlinear()
end

function coupling_time_gradient_type(::Type{QuasiSteady{1}}, ::Type{LiftingLineSection},
    ::Type{SimpleFlap}, ::Type{LiftingLineSectionControl})

    return Zeros()
end

# quasisteady with acceleration terms
function number_of_additional_parameters(::Type{QuasiSteady{2}}, ::Type{LiftingLineSection},
    ::Type{SimpleFlap}, ::Type{LiftingLineSectionControl})
    return 1
end

function coupling_inplaceness(::Type{QuasiSteady{2}}, ::Type{LiftingLineSection},
    ::Type{SimpleFlap}, ::Type{LiftingLineSectionControl})

    return OutOfPlace()
end

function coupling_rate_jacobian_type(::Type{QuasiSteady{2}}, ::Type{LiftingLineSection},
    ::Type{SimpleFlap}, ::Type{LiftingLineSectionControl})

    return Linear()
end

function coupling_state_jacobian_type(::Type{QuasiSteady{2}}, ::Type{LiftingLineSection},
    ::Type{SimpleFlap}, ::Type{LiftingLineSectionControl})

    return Nonlinear()
end

function coupling_parameter_jacobian_type(::Type{QuasiSteady{2}}, ::Type{LiftingLineSection},
    ::Type{SimpleFlap}, ::Type{LiftingLineSectionControl})

    return Nonlinear()
end

function coupling_time_gradient_type(::Type{QuasiSteady{2}}, ::Type{LiftingLineSection},
    ::Type{SimpleFlap}, ::Type{LiftingLineSectionControl})

    return Zeros()
end

# --- Methods --- #

# steady
function get_coupling_inputs(aero::QuasiSteady{0}, stru::LiftingLineSection,
    flap::SimpleFlap, ctrl::LiftingLineSectionControl, dx, x, p, t)
    # extract rate variables
    dvx, dvy, dvz, dωx, dωy, dωz, dδ = dx
    # extract state variables
    vx, vy, vz, ωx, ωy, ωz, δ = x
    # extract parameters
    a, b, a0, α0, cd0, cm0, cnδ, caδ, cmδ, ρ = p
    # local freestream velocity components
    u, v, ω = liftingline_velocities(vx, vz, ωy)
    # calculate loads
    Na, Aa, Ma = quasisteady0_loads(a, b, ρ, a0, α0, cd0, cm0, u, v)
    # loads due to flap deflections
    Nf, Af, Mf = simpleflap_loads(b, u, ρ, cnδ, caδ, cmδ, δ)
    # loads per unit span
    f = SVector(Aa+Af, 0, Na+Nf)
    m = SVector(0, Ma+Mf, 0)
    # return inputs
    return SVector(f..., m...)
end

# quasisteady without acceleration terms
function get_coupling_inputs(aero::QuasiSteady{1}, stru::LiftingLineSection,
    flap::SimpleFlap, ctrl::LiftingLineSectionControl, dx, x, p, t)
    # extract rate variables
    dvx, dvy, dvz, dωx, dωy, dωz, dδ = dx
    # extract state variables
    vx, vy, vz, ωx, ωy, ωz, δ = x
    # extract parameters
    a, b, a0, α0, cd0, cm0, cnδ, caδ, cmδ, ρ = p
    # local freestream velocity components
    u, v, ω = liftingline_velocities(vx, vz, ωy)
    # calculate aerodynamic loads
    Na, Aa, Ma = quasisteady1_loads(a, b, ρ, a0, α0, cd0, cm0, u, v, ω)
    # loads due to flap deflections
    Nf, Af, Mf = simpleflap_loads(b, u, ρ, cnδ, caδ, cmδ, δ)
    # loads per unit span
    f = SVector(Aa+Af, 0, Na+Nf)
    m = SVector(0, Ma+Mf, 0)
    # return inputs
    return SVector(f..., m...)
end

# quasisteady with acceleration terms
function get_coupling_inputs(aero::QuasiSteady{2}, stru::LiftingLineSection,
    flap::SimpleFlap, ctrl::LiftingLineSectionControl, dx, x, p, t)
    # extract rate variables
    dvx, dvy, dvz, dωx, dωy, dωz, dδ = dx
    # extract state variables
    vx, vy, vz, ωx, ωy, ωz, δ = x
    # extract parameters
    a, b, a0, α0, cd0, cm0, cnδ, caδ, cmδ, ρ = p
    # local freestream velocity components
    u, v, ω = liftingline_velocities(vx, vz, ωy)
    udot, vdot, ωdot = liftingline_accelerations(dvx, dvz, dωy)
    # calculate aerodynamic loads
    Na, Aa, Ma = quasisteady2_loads(a, b, ρ, a0, α0, cd0, cm0, u, v, ω, vdot, ωdot)
    # loads due to flap deflections
    Nf, Af, Mf = simpleflap_loads(b, u, ρ, cnδ, caδ, cmδ, δ)
    # loads per unit span
    f = SVector(Aa+Af, 0, Na+Nf)
    m = SVector(0, Ma+Mf, 0)
    # return inputs
    return SVector(f..., m...)
end

# --- Performance Overloads --- #

# TODO: State rate, state, and parameter jacobians

# --- Convenience Methods --- #

function set_additional_parameters!(padd, aero::QuasiSteady, stru::LiftingLineSection,
    flap::SimpleFlap, ctrl::LiftingLineSectionControl; rho)

    padd[1] = rho

    return padd
end

function separate_additional_parameters(aero::QuasiSteady, stru::LiftingLineSection,
    flap::SimpleFlap, ctrl::LiftingLineSectionControl, padd)

    return (rho = padd[1],)
end
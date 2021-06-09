"""
    Peters{N,TF,SV,SA} <: AbstractModel

Peter's finite state model with `N` state variables, inputs ``d = \\begin{bmatrix}
\\dot{\\theta} & \\ddot{h} & \\ddot{\\theta}\\end{bmatrix}^T`` and parameters
``p_a = \\begin{bmatrix} a & b & U & \\rho \\end{bmatrix}^T``
"""
struct Peters{N,TF,TV<:SVector{N,TF},TA<:SMatrix{N,N,TF}} <: AbstractModel
    A::TA
    b::TV
    c::TV
end

# --- Constructors --- #

"""
    Peters{N,TF=Float64}()

Initialize an object of type `Peters` which has `N` aerodynamic
degrees of freedom.
"""
Peters{N}() where N = Peters{N,Float64}()

function Peters{N,TF}() where {N,TF}

    b = zeros(TF, N)
    for n = 1:N-1
        b[n] = (-1)^(n-1)*factorial(big(N + n - 1))/factorial(big(N - n - 1))*
            1/factorial(big(n))^2
    end
    b[N] = (-1)^(N-1)

    c = zeros(TF, N)
    for n = 1:N
        c[n] = 2/n
    end

    d = zeros(TF, N)
    d[1] = 1/2

    D = zeros(TF, N, N)
    for m in 1:N-1
        n = m + 1
        D[n, m] = 1/(2*n)
    end
    for m in 2:N
        n = m - 1
        D[n, m] = -1/(2*n)
    end

    A = D + d*b' + c*d' + 1/2*c*b'

    return Peters(SMatrix{N,N,TF}(A), SVector{N,TF}(b), SVector{N,TF}(c))
end

# --- Traits --- #

number_of_states(::Type{Peters{N,TF,SV,SA}}) where {N,TF,SV,SA} = N
number_of_inputs(::Type{<:Peters}) = 3
number_of_parameters(::Type{<:Peters}) = 6
inplaceness(::Type{<:Peters}) = OutOfPlace()
mass_matrix_type(::Type{<:Peters}) = Constant()
state_jacobian_type(::Type{<:Peters}) = Varying()
input_jacobian_type(::Type{<:Peters}) = Varying()
input_dependence_type(::Type{<:Peters}) = Linear()

# --- Methods --- #

get_mass_matrix(model::Peters) = model.A

function get_rates(model::Peters{N,TF,SV,SA}, λ, d, p, t) where {N,TF,SV,SA}
    # extract aerodynamic states as statically sized vector
    λ = SVector{N}(λ)
    # extract structural deflections
    θdot, hddot, θddot = d
    # extract parameters
    a, b, U, ρ, a0, α0 = p
    # extract model constants
    cbar = model.c
    # calculate rates
    return peters_rhs(a, b, U, cbar, θdot, hddot, θddot, λ)
end

function get_state_jacobian(model::Peters, λ, d, p, t)
    # extract parameters
    a, b, U, ρ, a0, α0 = p
    # extract model constants
    cbar = model.c
    # jacobian with respect to aerodynamic states
    return peters_state_jacobian(b, U, cbar)
end

function get_input_jacobian(model::Peters, λ, d, p, t)
    # extract parameters
    a, b, U, ρ, a0, α0 = p
    # extract model constants
    cbar = model.c
    # return jacobian
    return peters_input_jacobian(a, b, U, cbar)
end

# TODO: Add parameter jacobian

# --- Coupled Model Properties --- #

# traits
inplaceness(::Type{<:Peters}, ::Type{TypicalSection}) = OutOfPlace()
mass_matrix_type(::Type{<:Peters}, ::Type{TypicalSection}) = Varying()
state_jacobian_type(::Type{<:Peters}, ::Type{TypicalSection}) = Varying()

# interface methods
function get_input_mass_matrix(aero::Peters{N,TF,SV,SA},
    stru::TypicalSection, u, p, t) where {N,TF,SV,SA}
    # extract aerodynamic parameters
    a, b, U, ρ, a0, α0 = p
    # construct submatrices
    Mda = zeros(SMatrix{3,N,TF})
    Mds = @SMatrix [0 0 0 0; 0 0 -1 0; 0 0 0 -1]
    Mra = zeros(SMatrix{2,N,TF})
    Mrs = hcat(
        zeros(SVector{2,TF}),
        zeros(SVector{2,TF}),
        -peters_loads_hddot(a, b, ρ),
        -peters_loads_θddot(a, b, ρ))
    # assemble mass matrix
    return [Mda Mds; Mra Mrs]
end

function get_inputs(aero::Peters{N,TF,SV,SA}, stru::TypicalSection,
    u, p, t) where {N,TF,SV,SA}
    # indices for extracting state variables
    iλ = SVector{N}(1:N)
    iq = SVector{4}(N+1:N+4)
    # separate aerodynamic and structural states
    λ = u[iλ]
    q = u[iq]
    # extract structural state variables
    h, θ, hdot, θdot = q
    # extract aerodynamic parameters
    a, b, U, ρ, a0, α0 = p
    # extract model constants
    bbar = aero.b
    # calculate aerodynamic loads
    L, M = peters_loads(a, b, U, ρ, a0, α0, bbar, θ, hdot, θdot, 0, 0, λ)
    # return portion of inputs that is not dependent on the state rates
    return SVector(θdot, 0, 0, L, M)
end

function get_input_state_jacobian(aero::Peters{N,TF,SV,SA},
    stru::TypicalSection, u, p, t) where {N,TF,SV,SA}
    # extract aerodynamic parameters
    a, b, U, ρ, a0, α0 = p
    # extract model constants
    bbar = aero.b
    # compute jacobian sub-matrices
    Jda = zeros(SMatrix{3,N,TF}) # d(d)/d(dλ)
    Jds = @SMatrix [0 0 0 1; 0 0 0 0; 0 0 0 0]
    Jra = peters_loads_λ(a, b, U, ρ, a0, bbar)
    Jrs = hcat(
        peters_loads_h(),
        peters_loads_θ(a, b, U, ρ, a0),
        peters_loads_hdot(a, b, U, ρ, a0),
        peters_loads_θdot(a, b, U, ρ, a0)
        )
    # return jacobian
    return [Jda Jds; Jra Jrs]
end

# TODO: Parameter jacobian

# --- Internal Methods --- #
peters_lhs(Abar, dλ) = Abar*dλ
peters_rhs(a, b, U, cbar, θdot, hddot, θddot,  λ) = cbar*(hddot + U*θdot + (b/2-a*b)*θddot) - U/b*λ
peters_mass_matrix(Abar) = Abar
peters_state_jacobian(b, U, cbar) = -U/b*Diagonal(one.(cbar))
peters_input_jacobian(a, b, U, cbar) = hcat(U*cbar, cbar, (b/2-a*b)*cbar)

function peters_loads(a, b, U, ρ, a0, α0, bbar, θ, hdot, θdot, hddot, θddot, λ)
    # circulatory load factor
    tmp1 = a0*ρ*U*b
    # non-circulatory load factor
    tmp2 = pi*ρ*b^3
    # constant based on geometry
    d = b/2 - a*b
    # induced flow velocity
    λ0 = 1/2 * bbar'*λ
    # lift at reference point
    L = tmp1*(U*θ + hdot + d*θdot - λ0 - U*α0) + tmp2*(hddot/b + U/b*θdot - a*θddot)
    # moment at reference point
    M = -tmp2*(hddot/2 + U*θdot + b*(1/8 - a/2)*θddot) + (b/2 + a*b)*L

    return SVector(L, M)
end

function peters_loads_λ(a, b, U, ρ, a0, bbar)
    tmp = a0*ρ*U*b
    L_λ = -tmp/2*bbar'
    M_λ = (b/2 + a*b)*L_λ
    return vcat(L_λ, M_λ)
end

function peters_loads_λdot(cbar)
    tmp = zero(cbar)'
    return vcat(tmp, tmp)
end

peters_loads_h() = SVector(0, 0)

function peters_loads_θ(a, b, U, ρ, a0)
    L_θ = a0*ρ*U^2*b
    M_θ = (b/2 + a*b)*L_θ
    return SVector(L_θ, M_θ)
end

function peters_loads_hdot(a, b, U, ρ, a0)
    L_hdot = a0*ρ*U*b
    M_hdot = (b/2 + a*b)*L_hdot
    return SVector(L_hdot, M_hdot)
end

function peters_loads_θdot(a, b, U, ρ, a0)
    tmp = pi*ρ*b^3
    L_θdot = a0*ρ*U*b*(b/2 - a*b) + tmp*U/b
    M_θdot = -tmp*U + (b/2 + a*b)*L_θdot
    return SVector(L_θdot, M_θdot)
end

function peters_loads_hddot(a, b, ρ)
    tmp1 = pi*ρ*b^3
    tmp2 = b/2 + a*b
    L_hddot = tmp1/b
    M_hddot = -tmp1/2 + tmp2*L_hddot
    return SVector(L_hddot, M_hddot)
end

function peters_loads_θddot(a, b, ρ)
    tmp1 = pi*ρ*b^3
    tmp2 = b/2 + a*b
    L_θddot = -tmp1*a
    M_θddot = -tmp1*(b/8 - a*b/2) + tmp2*L_θddot
    return SVector(L_θddot, M_θddot)
end

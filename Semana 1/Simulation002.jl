using LinearAlgebra, Distributions, Random, Plots

# 1. Model parameters (fixed-effects drift ⨯ mixed-effects diffusion)
α_d, β_d      = 35.3462, 0.0280        # fixed-effects diameter intercept & pull-back (Table 2 in the paper)
α_h, β_h      = 25.3301, 0.0396        # fixed-effects height intercept & pull-back (Table 2 in the paper)
σ11, σ12, σ22 = 2.3736, 0.6334, 0.3615  # mixed-effects diffusion entries (PD!). For the fixed ones,the matrix is singular

# 2. Simulation settings
Δt      = 1.0     # time step
T       = 100     # number of steps
n_paths = 5       # how many sample paths
t0      = 5.0     # initial “age”
x0      = [0.1, 1.3]  # start (diameter, height) at t0

# 3. Build OU matrices
Γ = Diagonal([β_d, β_h])
F = exp(-Γ * Δt)
Q = [
    σ11*(1 - exp(-2β_d*Δt))/(2β_d)         σ12*(1 - exp(-(β_d+β_h)*Δt))/(β_d+β_h);
    σ12*(1 - exp(-(β_d+β_h)*Δt))/(β_d+β_h) σ22*(1 - exp(-2β_h*Δt))/(2β_h)
]

# 4. Equilibrium level μ (fixed-effects only)
μ = [α_d, α_h]

# 5. Allocate & initialize
Xs = Array{Float64,3}(undef, 2, T+1, n_paths)
for j in 1:n_paths
    Xs[:,1,j] .= x0
end

# 6. Simulate all paths
rng  = MersenneTwister(123)
dist = MvNormal(zeros(2), Q)
for j in 1:n_paths, k in 1:T
    η = rand(rng, dist)
    Xs[:,k+1,j] = μ .+ F*(Xs[:,k,j] .- μ) .+ η
end

# 7. Plot
t = t0:Δt:(t0 + T*Δt)

p1 = plot(title="Diameter D(t)", xlabel="Age", ylabel="D (cm)", legend=false)
for j in 1:n_paths
    plot!(p1, t, Xs[1,:,j])
end

p2 = plot(title="Height H(t)", xlabel="Age", ylabel="H (m)", legend=false)
for j in 1:n_paths
    plot!(p2, t, Xs[2,:,j])
end

plot(p1, p2, layout=(2,1), size=(700,500),
     title="Fixed-Effects OU Simulation (PD Diffusion)")

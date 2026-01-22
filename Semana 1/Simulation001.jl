# single-path simulation of diameter and height processes

using LinearAlgebra, Distributions, Random, Plots

# ─────────────────────────────────────────────
# 1. Paper’s fixed‐effects parameter estimates. With these, the matrix is not positive definite.
# α_d, β_d, σ11  = 35.3462, 0.0280, 3.5968
# α_h, β_h, σ22  = 25.3301, 0.0396, 1.1968
# σ12           = 3.9187

# 1. Model parameters (fixed-effects drift ⨯ mixed-effects diffusion)
α_d, β_d      = 35.3462, 0.0280        # fixed-effects diameter intercept & pull-back (Table 2 in the paper)
α_h, β_h      = 25.3301, 0.0396        # fixed-effects height intercept & pull-back (Table 2 in the paper)
σ11, σ12, σ22 = 2.3736, 0.6334, 0.3615  # mixed-effects diffusion entries (PD!). For the fixed ones,the matrix is singular

Δt      = 1.0     # time step


# 2. Build drift + diffusion
Γ = Diagonal([β_d, β_h])
F = exp(-Γ * Δt)
Q = [ σ11*(1-exp(-2β_d*Δt))/(2β_d)     σ12*(1-exp(-(β_d+β_h)*Δt))/(β_d+β_h);
      σ12*(1-exp(-(β_d+β_h)*Δt))/(β_d+β_h)   σ22*(1-exp(-2β_h*Δt))/(2β_h) ]

# 3. Equilibrium μ (fixed-effects, φ=0)
μ = [α_d, α_h]   # ≈ [35.35, 25.33]

# 4. Initial condition at t0=5, below μ
t0 = 5
d0, h0 = 0.0, 1.3
X = zeros(2, T+1)
X[:,1] = [d0, h0]

# 5. Simulate forward from k=1..T
rng = MersenneTwister(39)
dist = MvNormal(zeros(2), Q)
for k in 1:T
    η = rand(rng, dist)
    X[:,k+1] = μ .+ F*(X[:,k] .- μ) .+ η
end

# 6. Plot the single trajectory
t = 0:Δt:T*Δt

p1 = plot(t, X[1,:], label="Diameter D(t)", xlabel="Time", ylabel="D", legend=:bottomright)
p2 = plot(t, X[2,:], label="Height   H(t)", xlabel="Time", ylabel="H", legend=:bottomright)

plot(p1, p2, layout=(2,1), size=(700,500), title="Single‐Tree Growth Simulation")

# #-D growth trajectory: time, diameter, height
using LinearAlgebra, Distributions, Random, Plots

# 1. Model parameters
β_d, β_h     = 0.1, 0.05           # pull-back rates
σ11, σ12, σ22 = 0.2, 0.1, 0.3      # diffusion covariance entries
α_d, φ_d     = 1.0, 10.0           # intercept & optimum for diameter
α_h, φ_h     = 0.5, 15.0           # intercept & optimum for height
Δt           = 1.0                 # time step
T            = 100                 # number of steps

# 2. Matrices
Γ = Diagonal([β_d, β_h])           # drift matrix
F = exp(-Γ * Δt)                   # discrete‐time multiplier
Q = [
    σ11*(1 - exp(-2β_d*Δt))/(2β_d)   σ12*(1 - exp(-(β_d+β_h)*Δt))/(β_d+β_h);
    σ12*(1 - exp(-(β_d+β_h)*Δt))/(β_d+β_h)   σ22*(1 - exp(-2β_h*Δt))/(2β_h)
]

# 3. Equilibrium
μ = [α_d + φ_d,  α_h + φ_h]

# 4. Simulate single path
X = zeros(2, T+1)
X[:,1] .= μ
rng = MersenneTwister(42)
dist = MvNormal(zeros(2), Q)
for k in 1:T
    η = rand(rng, dist)
    X[:,k+1] = μ .+ F*(X[:,k] .- μ) .+ η
end

# 5. Prepare data
t = 0:Δt:T*Δt
D = X[1, :]
H = X[2, :]

# 6. 3D plot: time on x, diameter on y, height on z
plt = plot(t, D, H,
    st = :path,          # connect points
    lw = 2,              # line width
    xlabel = "Time",
    ylabel = "Diameter",
    zlabel = "Height",
    title = "3D Growth Trajectory of a Single Tree",
    legend = false
)
display(plt)

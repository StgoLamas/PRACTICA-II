using Random, Distributions, LinearAlgebra
using Plots
gr()  # use GR backend

# ─── MODEL PARAMETERS ───────────────────────────────────────────────────────────
β_d, β_h      = 0.1, 0.05
σ11, σ12, σ22 = 0.2, 0.1, 0.3
α_d, φ_d      = 1.0, 10.0
α_h, φ_h      = 0.5, 15.0
Δt, T         = .01, 1000

function A(x)
    d, h = x
    return [β_d*(α_d-d), β_h*(α_h-h)] 
end

B = [σ11 σ12; σ12 σ22]
Ch = cholesky(B).U

# OU matrices
Γ = Diagonal([β_d, β_h])
F  = exp(-Γ * Δt)
Q  = [
  σ11*(1-exp(-2β_d*Δt))/(2β_d)   σ12*(1-exp(-(β_d+β_h)*Δt))/(β_d+β_h);
  σ12*(1-exp(-(β_d+β_h)*Δt))/(β_d+β_h)   σ22*(1-exp(-2β_h*Δt))/(2β_h)
]
μ = [α_d+φ_d, α_h+φ_h]

# ─── INITIAL CONDITION: start small! ───────────────────────────────────────────
D0, H0 = 0.1, 0.5    # sapling diameter & height at time 0
X = zeros(2, T+1)
X[:,1] = [D0, H0]    # << here is the key change

# 5. Simulate forward from k=1..T
rng = MersenneTwister(39)
dist = MvNormal(zeros(2), Δt)
for k in 1:T
    η = rand(rng, dist)
    X[:,k+1] = X[:,k] + A(X[:,k])*Δt .+ Ch * η
end

# ─── ANIMATE TAPERED “TREE” ────────────────────────────────────────────────────
anim = @animate for k in 1:(T+1)
  base_r = X[1,k]/2
  h      = X[2,k]
  taper_frac = 0.2

  θ = range(0, 2π; length=60)
  z = range(0, h;  length=30)
  x = [ (base_r*(1 - taper_frac*zj/h) * cos(tt)) for tt in θ, zj in z ]
  y = [ (base_r*(1 - taper_frac*zj/h) * sin(tt)) for tt in θ, zj in z ]
  zz= [ zj for _ in θ, zj in z ]

  surface(
    x, y, zz;
    color = :saddlebrown,
    xlims = (-φ_d, φ_d),
    ylims = (-φ_d, φ_d),
    zlims = (0, maximum(X[2,:])),
    legend = false,
    camera = (30, 30),
    xlabel = "X", ylabel = "Y", zlabel = "Z",
    title = "Tree Growth (t=$(k-1))"
  )
end

gif(anim, "tree_growth.gif", fps=60)
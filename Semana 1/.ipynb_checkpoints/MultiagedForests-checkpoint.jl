### Single-species, multiaged forest model, taken from "A Discussion of Vintage Optimization Models in Forest Economics",
# A. Piazza, Forest Science, Vol 66(4) pp.469-477, 2020. 

using JuMP
using SDDP
using Ipopt
using Gurobi
import Plots

const GRB_ENV = Gurobi.Env()

UB = 1000 # upper bound for the objective function
b = 0.975 # discount factor used by Buongiorno
α = 0.5 # risk aversion parameter for the utility function
n = 5 # number of age classes
f = [0.8, 4.2, 16.99, 68.1, 84] # biomass coeffcients
x0 = [1.0, 0.0, 0.0, 0.0, 0.0] # initial state

graph = SDDP.UnicyclicGraph(b; num_nodes = 1)

model = SDDP.PolicyGraph(
    graph;
    sense = :Max,
    upper_bound = UB,
    #optimizer = Ipopt.Optimizer,
    optimizer = () -> Gurobi.Optimizer(GRB_ENV),
) do sp, t
    @variable(sp, x[s = 1:n] >=0, SDDP.State, initial_value = x0[s]) # state variables, multidimensional. Time to figure this out: 3 hours!!!!
    @variable(sp, c[1:n] >= 0) # control variables
    @variable(sp, TV >= 0) # Timber volume
    @constraints(sp, begin
        c[n] == x[n].in # we harvest everything for the oldest age class
        [s = 1:n-1], x[s+1].out == x[s].in - c[s] # balance equations
        sum(f[s] * c[s] for s in 1:n) == TV # timber volume
        [s = 1:n], x[s].out <= 1 # state space constraints
        sum(x[s].out for s in 1:n) == 1 # x[1].out will be determined by this equation because the other x[s].out are determined by the balance equations
    end)
    @stageobjective(sp, TV^(1-α)/(1-α)) # CRRA objective function
    #@stageobjective(sp, log(TV)) # Logarithmic objective function
end

n_train_iterations = 50
n_simulations = 1

SDDP.train(model, iteration_limit = n_train_iterations)
# Simulations to compute bounds. These will have different lengths.
simulations = SDDP.simulate(model, n_simulations, [:x])
objective_values = [
    sum(stage[:stage_objective] for stage in sim) for sim in simulations
]
μ, ci = round.(SDDP.confidence_interval(objective_values, 1.96); digits = 2)
lower_bound = SDDP.calculate_bound(model)
println("Confidence interval: ", μ, " ± ", ci)
println("Lower bound: ", round(lower_bound, digits = 2))

# Simulations to compute bounds. These need to have the same length
simulations = SDDP.simulate(
    model,
    n_simulations,
    [:x];
    sampling_scheme = SDDP.InSampleMonteCarlo(;
        max_depth = 20,
        terminate_on_dummy_leaf = false,
    ),
)

Plots.plot(
    SDDP.publication_plot(simulations, title = "x[1]") do data
        return data[:x][1].out # area occupied by trees in  age class 1
    end,
    SDDP.publication_plot(simulations, title = "x[2]") do data
        return data[:x][2].out # area occupied by trees in  age class 2
    end,
    SDDP.publication_plot(simulations, title = "x[3]") do data
        return data[:x][3].out # area occupied by trees in  age class 3
    end,
    SDDP.publication_plot(simulations, title = "x[4]") do data
        return data[:x][4].out # area occupied by trees in  age class 4
    end,
    SDDP.publication_plot(simulations, title = "x[5]") do data
        return data[:x][5].out # area occupied by trees in  age class 5
    end,
    xlabel = "Stage",
    ylims = (0, 1),
    layout = (2, 3),
)

Plots.savefig("ForestState.pdf")

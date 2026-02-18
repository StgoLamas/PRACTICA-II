# SOLO PARA INICIALIZAR

samples = 20000

# Allocate & initialize
Xs = Array{Float64, 3}(undef, 2, n+1, samples)
X0 = reshape(repeat(x0, samples), (2, samples))
Xs[:,1,:] .= X0

# Simulate forward from k=1..n
rng = MersenneTwister(39)
dist = MvNormal(zeros(2), Δt*I2)
for k in 1:n
    for j in 1:samples
        η = rand(rng, dist)
        Xs[:,k+1,j] = Xs[:,k,j] .+ gamma*(mu .- Xs[:,k,j])*Δt .+ sqB * η
    end
end

# Plot the trajectories
t = 0:Δt:T
p1 = plot(t, Xs[1,:,1:20], xlabel="Time", ylabel="D", legend=false)
p2 = plot(t, Xs[2,:,1:20], xlabel="Time", ylabel="H", legend=false)
plot(p1, p2, layout=(2,1), size=(700,500), title="Single‐Tree Growth Simulation")

# Espaciado de las muestras para samplear cada 5 años
t_muestras = floor(Int,5/Δt)
# Clasificación de los árboles
Ys = max.(floor.(Int, Xs[:,1:t_muestras:end,:]/5),0)
 
# Se registran todas las posibles transiciones.
key_list = [([d1,h1], [d2,h2]) for d1 in 0:12, d2 in 0:12, h1 in 0:6, h2 in 0:6]
# Diccionario que contará todas las transiciones.
total_transiciones = Dict(key => 0 for key in key_list)
 
# Se cuentan todas las transiciones dada la simplificación.
for i in 1:samples
    for k in 1:size(Ys)[2]-1
        total_transiciones[Ys[:,k,i],Ys[:,k+1,i]] += 1
    end
end
 
sumas_transiciones = Dict()
valid_idx = []
for d1 in 0:9
    for h1 in 0:6
        s = sum(total_transiciones[[d1,h1],[d2,h2]] for d2 in 0:9, h2 in 0:6)
        # Omitir estados visitados menos del 1% del número de trayectorias
        if s >= samples/100 
            sumas_transiciones[[d1,h1]] = s
            push!(valid_idx, [d1,h1])
        end
    end
end

touch("markovian.txt")
touch("markovian_sum.txt")
io = open("markovian.txt", "w")
io_sum = open("markovian_sum.txt", "w")

mc = Dict()
mc_sum = Dict()
for ik in valid_idx
    s = sumas_transiciones[ik]
    print("Estado de partida: ")
    println(ik)
    println(" ")
    for jl in valid_idx
        mc[ik,jl] = trunc(round(total_transiciones[ik,jl]/s; digits=4); digits = 2)
        mc_str = string(ik,jl,mc[ik,jl],"\n") 
        write(io, mc_str)
        if mc[ik,jl] > 0
            print("Estado de llegada: ")
            print(jl)
            print(", prob. de transición = ")
            println(mc[ik,jl])
        end
    end
    mc_sum[ik] = sum(mc[ik,jl] for jl in valid_idx)
    mc_sum_str = string(ik,mc_sum[ik],"\n") 
    write(io_sum, mc_sum_str)
    print("Suma de probs. de transición = ")
    println(mc_sum[ik])
    println(" ")
    println("------------------")
    println(" ")
end

close(io)
close(io_sum)
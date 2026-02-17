# SOLO PARA INICIALIZAR
 
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
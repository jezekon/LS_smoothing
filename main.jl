# add packages if you need:
# add PlotlyJS, DelimitedFiles, LinearAlgebra, Statistics, BenchmarkTools, FLoops, SparseArrays
using PlotlyJS,
    DelimitedFiles,
    LinearAlgebra,
    Statistics,
    SparseArrays,
    BenchmarkTools,
    FLoops

include("subscripts/Subscripts_LS_2D.jl")
include("subscripts/Grafy2D.jl")

data = open("data/geometry.dat")
desities_raw = readlines(data)

# matrix dimension
Lx = Int((length(desities_raw[1]) - 1) / 2 + 1)
Ly = length(desities_raw)

# Create matrix
function Density_Matrix(Lx::Int, Ly::Int, desities_raw::Vector{String})
    densities = zeros(Int, (Lx, Ly))
    for i = 1:Lx, j = 1:Ly
        clen = desities_raw[i][(j-1)*2+1]
        if clen == '1'
            densities[i, j] = 1
        elseif clen == '0'
            densities[i, j] = 0
        end
    end
    return densities
end
dense_mat = Density_Matrix(Ly, Lx, desities_raw)

# Element density to nodal desity
function Nodal_desity(dense_mat::Matrix{Int64})
    aa, bb = size(dense_mat)
    # #Pomocné matice (rozsirene)
    zero_ver_shot = zeros(aa)
    zero_leng_long = zeros(bb + 1)'
    # Leva horni
    LH = vcat(zero_leng_long, hcat(zero_ver_shot, dense_mat))
    # Prava horni
    PH = vcat(zero_leng_long, hcat(dense_mat, zero_ver_shot))
    # Leva spotní
    LS = vcat(hcat(zero_ver_shot, dense_mat), zero_leng_long)
    # Prava spotní
    PS = vcat(hcat(dense_mat, zero_ver_shot), zero_leng_long)
    # # Matice na dělení:
    hrana = [1 fill(2.0, (1, bb - 1)) 1]
    stred = [2 fill(4.0, (1, bb - 1)) 2]
    stred_mat = stred
    for k = 1:(aa-2)
        stred_mat = vcat(stred_mat, stred)
    end
    delitel = vcat(vcat(hrana, stred_mat), hrana)
    # Složení:
    Hustota_uzly = (LH + PH + LS + PS) ./ delitel
    return Hustota_uzly
end
Hustota_uzly = Nodal_desity(dense_mat)


scale = 1
B = 1.
obl = 3
Obalka = range(start = -obl * B, stop = obl * B, length = 2 * obl + 1)

function build_A(Lx::Integer, Ly::Integer, obl::Real, B::Real)
    nnod = Lx * Ly
    A = spzeros(nnod, nnod)
    r_max = 2 * obl * B
    for i = 1:nnod
        Xᵢ = computeGridCoordinates(i, Ly, Lx)
        for j = 1:nnod
            Xⱼ = computeGridCoordinates(j, Ly, Lx)
            r = norm(Xⱼ - Xᵢ) * B
            if r <= r_max
                A[i, j] = radialBasisFunction(r, B)
            end
        end
    end
    return A
end

@time A = build_A(Lx, Ly, obl, B)
s = A \ vec(dense_mat)
scale = 1

x_r = range(start = first(0), stop = last(Lx), length = Lx * scale + 1)
y_r = range(start = first(0), stop = last(Ly), length = Ly * scale + 1)
okno = range(start = first(Obalka), stop = last(Obalka), length = (length(Obalka) - 1) * scale + 1)

Grid_y, Grid_x = mgrid(y_r,x_r)


A_matr = zeros(Int.(size(Grid_x) .+ (2 * obl) * scale))
        Grid_oy, Grid_ox = mgrid(okno, okno)
        rr = sqrt.((Grid_ox .^ 2) + (Grid_oy .^ 2))
        A_kop = exp.(-(rr ./ B) .^ 2)

        for i = 1:(Lx), j = 1:(Ly)
            Pos_i = RBF_Position(i, scale)
            Pos_j = RBF_Position(j, scale)
            Ar = Hustota_uzly[j, i] * A_kop
            A_matr[Pos_j[1]:Pos_j[2], Pos_i[1]:Pos_i[2]] += Ar
        end

# Rudukce matice na počáteční velikost
A_mat = A_matr[Int(obl/scale)+1:Int(size(A_matr,1)-obl/scale),Int(obl/scale)+1:Int(size(A_matr,2)-obl/scale)]

mean(dense_mat)
(Th_mat, hladina) = LS_Threshold_2D(A_mat, mean(dense_mat), 4.0)
HeatmapYellowBlack(A_mat)
HeatmapYellowBlack(Th_mat)
HeatmapYellowBlack(dense_mat.*1.0)

# writedlm("output.txt", Int.(Th_mat))

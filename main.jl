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
    LH = vcat(zero_leng_long, hcat(zero_ver_shot, dense_mat))  # Leva horni
    PH = vcat(zero_leng_long, hcat(dense_mat, zero_ver_shot))  # Prava horni
    LS = vcat(hcat(zero_ver_shot, dense_mat), zero_leng_long)  # Leva spotní
    PS = vcat(hcat(dense_mat, zero_ver_shot), zero_leng_long)  # Prava spotní
    # # Matice na dělení:
    hrana = [1 fill(2.0, (1, bb - 1)) 1]
    stred = [2 fill(4.0, (1, bb - 1)) 2]
    stred_mat = zeros(aa-1, bb+1)
    for k = 1:(aa-1)
        stred_mat[k,:] = stred
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
# @time A = build_A(Lx+1, Ly+1, obl, B)
# s = A \ vec(Hustota_uzly)
# s_mat = reshape(s, (Ly+1, Lx+1))  # from system of equations
s_mat = Hustota_uzly            # form nodal density


x_r = range(start = first(0), stop = last(Lx), length = Lx * scale + 1)
y_r = range(start = first(0), stop = last(Ly), length = Ly * scale + 1)
okno = range(start = first(Obalka), stop = last(Obalka), length = (length(Obalka) - 1) * scale + 1)
Grid_y, Grid_x = mgrid(y_r,x_r)

function BuildLevelSetFunction(
    okno::StepRangeLen,
    B::Float64,
    Lx::Int,#Grid::BoxOfStructuredGrid,
    Ly::Int,
    scale::Int,
    obl::Int,
    vyska_mat::Matrix,
    FGx::Matrix,
)
    @floop begin
        A_matr = zeros(Int.(size(FGx) .+ (2 * obl) * scale))
        Grid_oy, Grid_ox = mgrid(okno, okno)
        rr = sqrt.((Grid_ox .^ 2) + (Grid_oy .^ 2))
        A_kop = exp.(-(rr ./ B) .^ 2)
        # for i = 1:(Grid.NoW[1]), j = 1:(Grid.NoW[2])
        for i = 1:(Lx), j = 1:(Ly)
            Pos_i = RBF_Position(i, scale)
            Pos_j = RBF_Position(j, scale)
            Ar = vyska_mat[j, i] * A_kop
            A_matr[Pos_j[1]:Pos_j[2], Pos_i[1]:Pos_i[2]] += Ar
        end
    end
    return A_mat = A_matr[
        Int(obl * scale)+1:Int(size(A_matr, 1) - obl * scale),
        Int(obl * scale)+1:Int(size(A_matr, 2) - obl * scale),
    ]
end
A_mat = @time(BuildLevelSetFunction(okno, B, Lx, Ly, scale, obl, s_mat, Grid_x))

# Element density
El_dense = zeros(size(dense_mat))
for i = 1:Lx, j = 1:Ly
    El_dense[j,i] =  mean([A_mat[j,i], A_mat[j,i+1], A_mat[j+1,i], A_mat[j+1,i+1]])
end

# Finding the level to maintain the volume ratio:
(Th_mat, hladina) = LS_Threshold_2D(El_dense, mean(dense_mat), 4.0)

# Plot:
# HeatmapYellowBlack(A_mat)
HeatmapYellowBlack(Th_mat)
# HeatmapYellowBlack(El_dense)
# HeatmapYellowBlack(dense_mat.*1.0)

# Write data:
# writedlm("output.txt", Int.(Th_mat))

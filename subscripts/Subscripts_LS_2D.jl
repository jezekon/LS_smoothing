# restart Julia 
function restart()
    startup = """
        Base.ACTIVE_PROJECT[]=$(repr(Base.ACTIVE_PROJECT[]))
        Base.HOME_PROJECT[]=$(repr(Base.HOME_PROJECT[]))
        cd($(repr(pwd())))
        """
    cmd = `$(Base.julia_cmd()) -ie $startup`
    atexit(()->run(cmd))
    exit(0)
end

# Point value of RBF
function radialBasisFunction(r::Float64, B::Float64)::Float64
    return exp(-(r / B)^2)
end

# Grid coordinates 3D
function computeGridCoordinates(i::Int, Lx::Int, Ly::Int)::Vector{Float64}
    xᵢ = (mod(i - 1, Lx))
    yᵢ = (Int(floor((i - 1) / Lx)))
    return [xᵢ, yᵢ]
end

function RBF_Position(i::Int, scale::Int)
    L = Int((i - 1) * scale + 1)
    R = Int((i - 1 + 2 * obl) * scale + 1)
    return [L, R]
end

# Level-set threshold for specific volume ratio
function LS_Threshold_2D(s_mat::Matrix, V_frac::Float64, Exp::Float64)
    Eps = 1.0
    # Hledání "levelu" pro dosažení pořebného objemového poměru:
    Th_mat = zeros(size(s_mat))
    th = maximum(s_mat) / 2 # -> odhad
    krok = maximum(s_mat) / 10 # počáteční krok
    smer = 0
    n = 0 # počáteční směr není znám

    while n < 40 && Eps > 10^(-Exp)
        # global Volume, Eps, krok, th, smer, n
        # Výpočet objemu pro aktuální threshold (odchylky)
        for i = 1:size(s_mat, 1), j = 1:size(s_mat, 2)
            if s_mat[i, j] >= th
                Th_mat[i, j] = 1
            else
                Th_mat[i, j] = 0
            end
        end
        Volume = sum(Th_mat) / (size(s_mat, 1) * size(s_mat, 2))
        Eps = abs(V_frac - Volume) # --> odchylka
        # Změna kroku v závislosti na směru
        if (Volume - V_frac) > 0 && smer == 2
            krok = krok / 2
        elseif (Volume - V_frac) < 0 && smer == 1
            krok = krok / 2
        end
        # Posouvání levelu
        if (Volume - V_frac) > 0
            th = th + krok
            smer = 1
        elseif (Volume - V_frac) < 0
            th = th - krok
            smer = 2
        end
        n = n + 1
    end
    return [Th_mat, th]
end


## A structured grid covering the area in which an unstructured grid is located.
# function Structured_mesh2D(X::Matrix, vs::Float64, obl::Int)
#     Box_step = nothing
#     for i = 1:length(X[:, 1])
#         range_X = [minimum(X[i, :]), maximum(X[i, :])]
#         st = range((range_X[1]-obl*vs), step = vs, stop = (range_X[2]+obl*vs))
#         # st = range(range_X[1], step = vs, stop = range_X[2] + vs) # rozšíření oblasti o jeden element
#         if Box_step == nothing
#             Box_step = [st]
#         else
#             Box_step = vcat(Box_step, [st])
#         end
#     end
#     L_box = abs.(last.(Box_step) .- first.(Box_step))
#     N_box = length.(Box_step)
#     Grid_y, Grid_x = mgrid(Box_step[2], Box_step[1])
#     Grid = BoxOfStructuredGrid(
#         Box_step[1],
#         Box_step[2],
#         Grid_x,
#         Grid_y,
#         L_box,
#         N_box,
#     )
#     return Grid
# end



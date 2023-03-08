## Heatmap for Level-set:
function HeatmapYellowBlack(A_mat::Matrix{Float64})
    a, b = size(A_mat)
    ratio = b/a
    A_mat = A_mat ./ maximum(A_mat)
    data = heatmap(;
        z = A_mat,#dense_mat, A_mat, HF
        colorbar = attr(;
            title = "Range",
            titleside = "right",
            titlefont = attr(; size = 14, family = "Arial, sans-serif"),
        ),
        colorscale = [
            [0.0, "rgb(0,0,0)"],
            [0.2, "rgb(65, 10, 104)"],
            [0.4, "rgb(147, 38, 103)"],
            [0.6, "rgb(221, 81, 58)"],
            [0.8, "rgb(252, 165, 10)"],
            [1.0, "rgb(249, 252, 158)"],
        ],
    )
    layout = Layout(
        width =  1000,
        height = 1000/ratio,
        title = "LSF of the MBB beam",
        titlefont = attr(; size = 22),#, family="Arial, sans-serif"),
        xaxis = attr(
            title = "",
            showgrid = false,
            zeroline = false,
            showticklabels = false,
        ),
        yaxis = attr(title = "", zeroline = false, showticklabels = false),
    )
    return plot(data, layout)
end

## Heatmap for threshold:
function HeatmapWhiteBlue(Th_mat::Matrix{Float64})
    a, b = size(Th_mat)
    ratio = b/a
    data = heatmap(;
        z = Th_mat,
        showscale = false,
        colorscale = [
            [0.0, "rgb(250,250,250)"], [1.0, "rgb(27,75,121)"],
        ],
    )

    layout = Layout(
        width =  1000,
        height = 1000/ratio,
        xaxis = attr(
            title = "",
            showgrid = false,
            zeroline = false,
            showticklabels = false,
        ),
        yaxis = attr(title = "", showgrid = false, zeroline = false, showticklabels = false),
    )
    return plot(data, layout)
end

## Obálka (LS = 0):
function ScatterObalka(souradnice::Matrix{Float64})
    trace =
        scatter(; x = souradnice[:, 1], y = souradnice[:, 2], mode = "markers")

    layout = Layout(
        width = 1000,
        height = 400,
        title = "Contour",
        titlefont = attr(; size = 22),
        xaxis = attr(title = "x-axis", showgrid = true, zeroline = false),
        yaxis = attr(title = "y-axis", zeroline = false),
    )

    return plot([trace], layout)
end

## Mesh - porovnání strukturované a nestrukturované sítě.
function ComparisonOfGrids(
    X::Matrix,
    IEN::Matrix,
    color_map::Vector{Float64},
    Grid_x::Matrix{Float64},
    Grid_y::Matrix{Float64},
    dense_npr::Matrix{Float64},
    Xe::Matrix{Float64},
    Ye::Matrix{Float64},
)
    pos_x = X[1, :]
    pos_y = X[2, :]
    # Plot
    edges_trace = [
        scatter(
            mode = "lines",
            x = Xe[i, :],
            y = Ye[i, :],
            line = attr(width = 0.5, color = "#888"),
        ) for i in [1:1:length(Xe[:, 1]);]
    ]

    # Create nodes
    nodes_trace = scatter(
        x = pos_x,
        y = pos_y,
        mode = "markers",
        text = [string("density: ", dens) for dens in vec(color_map)],
        marker = attr(
            showscale = true,
            # colorscale=colors.imola,
            colorscale = [
                [0.0, "rgb(240,249,33)"],
                [0.2, "rgb(251,165,54)"],
                [0.4, "rgb(224,100,97)"],
                [0.6, "rgb(176,43,143)"],
                [0.8, "rgb(106,1,166)"],
                [1.0, "rgb(13,7,134)"],
            ],
            color = color_map,
            size = 10,
            line_width = 2,
            colorbar = attr(
                thickness = 15,
                title = "Nodal density",
                xanchor = "left",
                titleside = "right",
            ),
        ),
    )

    # Structured grid
    structured_grid = scatter(
        x = vec(Grid_x),
        y = vec(Grid_y),
        mode = "markers",
        text = [string("density: ", dens) for dens in vec(dense_npr)],
        marker = attr(
            showscale = false,
            colorscale = [
                [0.0, "rgb(240,249,33)"],
                [0.2, "rgb(251,165,54)"],
                [0.4, "rgb(224,100,97)"],
                [0.6, "rgb(176,43,143)"],
                [0.8, "rgb(106,1,166)"],
                [1.0, "rgb(13,7,134)"],
            ],
            color = vec(dense_npr),
            size = 10,
            symbol = "square",
            # line_width=2,
            line = attr(width = 2, color = "#2b2930"),
        ),
    )

    # Create Plot
    # return plot([edges_trace, structured_grid, nodes_trace],
    return plot(
        vcat(edges_trace, nodes_trace, structured_grid),
        Layout(
            # width=1200, height=450,
            hovermode = "closest",
            title = "Nodal density of an unstructured grid",
            titlefont_size = 21,
            showlegend = false,
            showarrow = false,
            xaxis = attr(
                showgrid = false,
                zeroline = false,
                showticklabels = false,
            ),
            yaxis = attr(
                showgrid = false,
                zeroline = false,
                showticklabels = false,
            ),
            plot_bgcolor = "white",
        ),
    )
end

## Grid with nodal density
function GridWithNodes(
    X::Matrix,
    IEN::Matrix,
    color_map::Vector{Float64},
    Xe::Matrix{Float64},
    Ye::Matrix{Float64},
)
    delka = 1200
    pos_x = X[1, :]
    pos_y = X[2, :]
    X_range = maximum(Xe) - minimum(Xe)
    Y_range = maximum(Ye) - minimum(Ye)
    # Plot
    edges_trace = [
        scatter(
            mode = "lines",
            x = Xe[:, i],
            y = Ye[:, i],
            line = attr(width = 1, color = "#888"),
        ) for i in [1:1:length(Xe[1, :]);]
    ]

    # Create nodes
    nodes_trace = scatter(
        x = pos_x,
        y = pos_y,
        mode = "markers",
        marker = attr(
            showscale = true,
            # colorscale=colors.imola,
            colorscale = [
                [0.0, "rgb(240,249,33)"],
                [0.2, "rgb(251,165,54)"],
                [0.4, "rgb(224,100,97)"],
                [0.6, "rgb(176,43,143)"],
                [0.8, "rgb(106,1,166)"],
                [1.0, "rgb(13,7,134)"],
            ],
            color = color_map,
            size = 10,
            line_width = 1.2,
            colorbar = attr(
                thickness = 15,
                title = "Nodal density",
                xanchor = "left",
                titleside = "right",
            ),
        ),
    )

    # Create Plot
    return plot(
        vcat(edges_trace, nodes_trace),
        Layout(
            width = delka,
            height = delka/(X_range/Y_range),
            hovermode = "closest",
            title = "Nodal density of an unstructured grid",
            titlefont_size = 21,
            showlegend = false,
            showarrow = false,
            xaxis = attr(
                showgrid = false,
                zeroline = false,
                showticklabels = false,
            ),
            yaxis = attr(
                showgrid = false,
                zeroline = false,
                showticklabels = false,
            ),
            plot_bgcolor = "white",
        ),
    )
end

function DenseInGPsGrid(
    X::Matrix,
    Xe::Matrix{Float64},
    Ye::Matrix{Float64},
    Gauss::Matrix{Float64},
    rho::Vector{Float64},
)
    pos_x = X[1, :]
    pos_y = X[2, :]
    # Plot
    edges_trace = [
        scatter(
            mode = "lines",
            x = Xe[i, :],
            y = Ye[i, :],
            line = attr(width = 1, color = "#888"),
        ) for i in [1:1:length(Xe[:, 1]);]
    ]

    # Create nodes
    nodes_trace = scatter(
        x = Gauss[:, 1],
        y = Gauss[:, 2],
        mode = "markers",
        text = [string("density: ", dens) for dens in rho],
        marker = attr(
            showscale = true,
            # colorscale=colors.imola,
            colorscale = [
                [0.0, "rgb(240,249,33)"],
                [0.2, "rgb(251,165,54)"],
                [0.4, "rgb(224,100,97)"],
                [0.6, "rgb(176,43,143)"],
                [0.8, "rgb(106,1,166)"],
                [1.0, "rgb(13,7,134)"],
            ],
            color = rho,
            size = 10,
            line_width = 1.2,
            colorbar = attr(
                thickness = 15,
                title = "Nodal density",
                xanchor = "left",
                titleside = "right",
            ),
        ),
    )

    return plot(
        vcat(edges_trace, nodes_trace),
        Layout(
            width = 1200,
            height = 450,
            hovermode = "closest",
            title = "Nodal density of an unstructured grid",
            titlefont_size = 21,
            showlegend = false,
            showarrow = false,
            xaxis = attr(
                showgrid = false,
                zeroline = false,
                showticklabels = false,
            ),
            yaxis = attr(
                showgrid = false,
                zeroline = false,
                showticklabels = false,
            ),
            plot_bgcolor = "white",
        ),
    )
end

function LeastSquareNodeFit(
    Xe::Matrix,
    Ye::Matrix,
    Gauss::Matrix,
    ρ::Matrix,
    mid_x::Matrix,
    mid_y::Matrix,
    Plocha_x::StepRangeLen,
    Plocha_y::StepRangeLen,
    Fρ::Matrix,
    bod::Float64,
    eye::Vector,
)
    # Nepravidelná síť (Xe, Ye, Ze)
    # GP - šedé (Gauss)
    # čárkovaně midsides (mid_x, mid_y)
    # plocha (Gauss, Fρ)
    # promítnutí
    edges_trace = [
        scatter(
            mode = "lines",
            x = Xe[i, :],
            y = Ye[i, :],
            z = zeros(length(Ye[i, :])),
            line = attr(width = 4.0, color = "#888"),
            type = "scatter3d",
        ) for i in [1:1:length(Xe[:, 1]);]
    ]
    # edges_trace = [
    #     scatter(
    #         mode = "lines",
    #         x = Xe[:, i],
    #         y = Ye[:, i],
    #         z = zeros(length(Ye[:, i])),
    #         line = attr(width = 4.0, color = "#888"),
    #         type = "scatter3d",
    #     ) for i in [1:1:length(Xe[1, :]);]
    # ]

    trace1 = scatter(
        x = Gauss[:, 1],
        y = Gauss[:, 2],
        z = zeros(length(Gauss[:, 1])),
        mode = "markers",
        marker = attr(
            opacity = 0.8,
            size = 8,
            colorscale = [
                [0.0, "rgb(240,249,33)"],
                [0.2, "rgb(251,165,54)"],
                [0.4, "rgb(224,100,97)"],
                [0.6, "rgb(176,43,143)"],
                [0.8, "rgb(106,1,166)"],
                [1.0, "rgb(13,7,134)"],
            ],
            color = vec(ρ),
        ),
        type = "scatter3d",
    )

    trace2 = scatter(
        x = Gauss[:, 1],
        y = Gauss[:, 2],
        z = vec(ρ),
        mode = "markers",
        marker = attr(
            # opacity=0.8,
            size = 8,
            colorscale = [
                [0.0, "rgb(240,249,33)"],
                [0.2, "rgb(251,165,54)"],
                [0.4, "rgb(224,100,97)"],
                [0.6, "rgb(176,43,143)"],
                [0.8, "rgb(106,1,166)"],
                [1.0, "rgb(13,7,134)"],
            ],
            color = vec(ρ),
        ),
        type = "scatter3d",
    )

    trace3 = scatter(
        x = vcat(Xe[1, 3], Xe[1, 3]),
        y = vcat(Ye[1, 3], Ye[1, 3]),
        z = vcat(bod, 0),
        mode = "markers",
        marker = attr(
            opacity = 0.8,
            size = 6,
            symbol = 'x',
            # color = "rgb(0,152,0)",
            colorscale = [
                [0.0, "rgb(240,249,33)"],
                [0.2, "rgb(251,165,54)"],
                [0.4, "rgb(224,100,97)"],
                [0.6, "rgb(176,43,143)"],
                [0.8, "rgb(106,1,166)"],
                [1.0, "rgb(13,7,134)"],
            ],
            color = bod,
        ),
        type = "scatter3d",
    )

    midsides = [
        scatter(
            mode = "lines",
            x = mid_x[i, :],
            y = mid_y[i, :],
            z = zeros(length(mid_x[i, :])),
            line = attr(width = 2.0, color = "#888"),
            dash = "dashdot",
            type = "scatter3d",
        ) for i in [1:1:length(mid_y[:, 1]);]
    ]

    Grid_y, Grid_x = mgrid(Plocha_y, Plocha_x)
    plocha = surface(
        x = Grid_x,
        y = Grid_y,
        z = Fρ',
        # colorscale = [[0.0, "rgb(152,152,196)"],[1.0, "rgb(152,152,196)"]],
        colorscale = [
            [0.0, "rgb(240,249,33)"],
            [0.2, "rgb(251,165,54)"],
            [0.4, "rgb(224,100,97)"],
            [0.6, "rgb(176,43,143)"],
            [0.8, "rgb(106,1,166)"],
            [1.0, "rgb(13,7,134)"],
        ],
        color = Fρ,
        opacity = 0.5,
        showscale = false,
    )
    Promitnuti = [
        scatter(
            mode = "lines",
            x = (vcat(Gauss[i, 1], Gauss[i, 1])),
            y = (vcat(Gauss[i, 2], Gauss[i, 2])),
            z = (vcat(0.0, ρ[i])),
            line = attr(width = 2.0, color = "#888"),
            type = "scatter3d",
        ) for i in [1:1:length(Gauss[:, 1]);]
    ]

    Promitnuti2 = scatter(
        mode = "lines",
        x = vcat(Xe[1, 3], Xe[1, 3]),
        y = vcat(Ye[1, 3], Ye[1, 3]),
        z = vcat(bod, 0),
        line = attr(width = 2.0, color = "#888"),
        type = "scatter3d",
    )

    layout = Layout(
        width = 900,
        height = 800,
        showlegend = false,
        scene = attr(
            camera = attr(eye = attr(x = eye[1], y = eye[1], z = eye[1])),
            xaxis = attr(
                title = "x",
                titlefont = attr(; size = 20),
                zeroline = false,
                showticklabels = false,
            ),
            yaxis = attr(
                title = "y",
                titlefont = attr(; size = 20),
                zeroline = false,
                showticklabels = false,
            ),
            zaxis = attr(
                title = "ᵨ",
                titlefont = attr(; size = 25),
                zeroline = false,
            ),
            aspectratio = attr(x = 1, y = 1, z = 0.7),
            aspectmode = "manual",
        ),
        margin = attr(l = 0, r = 0, b = 0, t = 0),
    )

    return plot(
        vcat(
            edges_trace,
            trace1,
            trace2,
            trace3,
            midsides,
            plocha,
            Promitnuti,
            Promitnuti2,
        ),
        layout,
    )
end

## Mesh - porovnání strukturované a nestrukturované sítě.
function StructuredGrids(
    Grid_x::Matrix{Float64},
    Grid_y::Matrix{Float64},
    dense_npr::Matrix{Float64},
    Xe::Matrix{Float64},
    Ye::Matrix{Float64},
)
    a, b = size(A_mat)
    ratio = b/a
    # Plot
    edges_trace = [
        scatter(
            mode = "lines",
            x = Xe[i, :],
            y = Ye[i, :],
            line = attr(width = 1.2, color = "#888"),
        ) for i in [1:1:length(Xe[:, 1]);]
    ]

    # Structured grid
    structured_grid = scatter(
        x = vec(Grid_x),
        y = vec(Grid_y),
        mode = "markers",
        text = [string("density: ", dens) for dens in vec(dense_npr)],
        marker = attr(
            showscale = true,
            colorscale = [
                [0.0, "rgb(240,249,33)"],
                [0.2, "rgb(251,165,54)"],
                [0.4, "rgb(224,100,97)"],
                [0.6, "rgb(176,43,143)"],
                [0.8, "rgb(106,1,166)"],
                [1.0, "rgb(13,7,134)"],
            ],
            color = vec(dense_npr),
            size = 9,
            symbol = "square",
            # line_width=2,
            line = attr(width = 1.2, color = "#2b2930"),
            colorbar = attr(
                thickness = 15,
                title = "Nodal density",
                xanchor = "left",
                titleside = "right",
            ),
        ),
    )

    # Create Plot
    return plot(
        vcat(edges_trace, structured_grid),
        Layout(
            width =  1000,
            height = 1000/ratio,
            hovermode = "closest",
            title = "Nodal density of an structured grid",
            titlefont_size = 21,
            showlegend = false,
            showarrow = false,
            xaxis = attr(
                showgrid = false,
                zeroline = false,
                showticklabels = false,
            ),
            yaxis = attr(
                showgrid = false,
                zeroline = false,
                showticklabels = false,
            ),
            plot_bgcolor = "white",
        ),
    )
end

## Level-set projection
function LSprojection(
    s_mat::Matrix,
    Th_mat::Matrix,
    eye::Vector,
    velikost::Vector,
    hladina::Float64
)
    x_eye = eye[1]
    y_eye = eye[2]
    z_eye = eye[3]

    a, b = size(A_mat)
    ratio = b/a

    for i in 1:size(Th_mat,1), j in 1:size(Th_mat,2)
        if Th_mat[i,j] == 0
            Th_mat[i,j] = NaN
        end
    end

    s_mat = (s_mat./maximum(s_mat)).*2 .-1
    hladina = (hladina./maximum(s_mat)).*2 .-1

    for i in 1:size(s_mat,1), j in 1:size(s_mat,2)
        if s_mat[i,j] < - 0.95
            s_mat[i,j] = NaN
        end
    end

    trace1 = surface(;
        # z = s_mat'.+hladina,
        z = s_mat'.+0.35,
        # z = s_mat' .- 0.469,
        colorscale = [[0.0, "rgb(107,27,94)"], [1.0, "rgb(49,54,149)"]],
        showscale = false,
    )
    trace2 = surface(;
        z = Th_mat' .- 5,
        colorscale = [[0.0, "rgb(37,83,125)"], [1.0, "rgb(37,83,125)"]],
        showscale = false,
    )

    trace3 = surface(;
        z = zeros(size(Th_mat))',
        colorscale = [[0.0, "rgb(152,152,196)"], [1.0, "rgb(152,152,196)"]],
        opacity = 0.5,
        showscale = false,
    )

    layout = Layout(
        width = velikost[1],
        height = velikost[2],
        scene = attr(
            camera = attr(eye = attr(x = x_eye, y = y_eye, z = z_eye)),
            xaxis = attr(title = "", zeroline = false, showticklabels = false),
            yaxis = attr(title = "", zeroline = false, showticklabels = false),
            zaxis = attr(
                title = "ᵨ",
                titlefont = attr(; size = 30),
                zeroline = false,
            ),
            aspectratio = attr(x = 1, y = 1/ratio, z = 0.4),
            aspectmode = "manual",
        ),
        margin = attr(l = 0, r = 0, b = 0, t = 0),
    )
    return p = plot([trace1, trace2, trace3], layout)
end

## RBFs in nodes
function RBFsNodes(
    Grid_x::Matrix,
    Grid_y::Matrix,
    A::Matrix,
    s::Matrix,
    Lx::Int,
    Ly::Int,
)
    x_eye = eye[1]
    y_eye = eye[2]
    z_eye = eye[3]

    data = [
        surface(
            x = Grid_x .+ (mod(i - 1, Lx) + 1),
            y = Grid_y .+ (Int(floor((i - 1) / Lx)) + 1),
            z = A .* s[(Int(floor((i - 1) / Lx))+1), (mod(i - 1, Lx)+1)],
            # colorscale = [[0.0, "rgb(235,235,235)"],[1.0, "rgb(100,13,13)"]],
            colorscale = [
                [0.0, "rgb(235,235,235)"],
                [
                    1.0,
                    "rgb($(235-135*s[Int(floor((i-1)/Lx))+1, mod(i-1,Lx)+1]),$(235-222*s[Int(floor((i-1)/Lx))+1, mod(i-1,Lx)+1]),$(235-222*s[Int(floor((i-1)/Lx))+1, mod(i-1,Lx)+1]))",
                ],
            ],
            showscale = false,
        ) for i in [1:1:length(s);]
    ]

    layout = Layout(;
        title = "Gaussian Radial Basis Function",
        width = 800,
        height = 800,
        scene = attr(
            camera = attr(eye = attr(x = x_eye, y = y_eye, z = z_eye)),
            xaxis = attr(title = "", zeroline = false, showticklabels = false),
            yaxis = attr(title = "", zeroline = false, showticklabels = false),
            zaxis = attr(title = "", zeroline = false),
            aspectratio = attr(x = 1, y = (Ly / Lx), z = 0.2),
            aspectmode = "manual",
        ),
        margin = attr(l = 0, r = 0, b = 0, t = 0),
    )
    return plot(data, layout)
end

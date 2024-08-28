using CSV, DataFrames

N = 30
gr = (1 + 5^0.5) / 2

i = 0:(N-1)
phi = 2 * pi / gr .* i
theta = acos.(1 .- 2 / N .* i)

# df = DataFrame(phi=phi[:], theta=theta[:])
# println(df)
x = sin.(theta) .* cos.(phi)
y = sin.(theta) .* sin.(phi)
z = cos.(theta)
df = DataFrame(x=x, y=y, z=z)
CSV.write("test/data/points.csv", df)

using PlotlyJS
plot(
    scatter3d(x=x, y=y, z=z, mode="markers"),
    Layout(
        scene=attr(
            xaxis=attr(title="x"),
            yaxis=attr(title="y"),
            zaxis=attr(title="z")
        ),
        yaxis=attr(scaleanchor="x", scaleratio=1),
        zaxis=attr(scaleanchor="x", scaleratio=1),
    )
)

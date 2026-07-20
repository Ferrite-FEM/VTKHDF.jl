# # ImageData: Mandelbrot volume
#
# Port of the `mandelbrot-vti.hdf` example from the
# [VTKHDF specification](https://docs.vtk.org/en/latest/vtk_file_formats/vtkhdf_file_format/index.html):
# an ImageData volume sampling the Mandelbrot iteration count. The x/y axes
# span the complex c-plane and the z axis varies the real part of the
# starting value, so every slice is a different member of the Mandelbrot
# family. Coordinate ranges define the uniform grid (origin and spacing).
# `Iterations` volume-rendered in ParaView:
#
# ![Volume rendering of the Mandelbrot iteration count](../assets/examples/image_data-light.png)
# ![Volume rendering of the Mandelbrot iteration count](../assets/examples/image_data-dark.png)

using VTKHDF

x = range(-1.75, 0.75; length = 20)
y = range(-1.25, 1.25; length = 21)
z = range(0.0, 2.0; length = 22)
nothing #hide

function iterations(c, w)
    for n in 1:100
        abs2(w) > 4 && return Float32(n)
        w = w * w + c
    end
    return Float32(100)
end

iters = [iterations(complex(xi, yi), complex(zi, 0)) for xi in x, yi in y, zi in z];

# A central-difference gradient of the iteration count, as a
# `(3, nx, ny, nz)` vector field:

function gradient(f, spacing)
    g = zeros(3, size(f)...)
    R = CartesianIndices(f)
    for d in 1:3, I in R
        δ = CartesianIndex(ntuple(==(d), 3))
        hi, lo = min(I + δ, last(R)), max(I - δ, first(R))
        g[d, I] = (f[hi] - f[lo]) / ((hi[d] - lo[d]) * spacing[d])
    end
    return g
end
nothing #hide

# Point data shaped like the grid is detected automatically; `attribute`
# marks the active scalars/vectors.

vtkhdf_grid("mandelbrot", x, y, z) do vtk
    vtk["Iterations", VTKPointData(), attribute = :Scalars] = iters
    vtk["IterationsGradient", VTKPointData(), attribute = :Vectors] =
        gradient(iters, step.((x, y, z)))
end
nothing #hide

# ## Reading it back
#
# The grid metadata comes back through `grid_info`:

r_mandel = vtkhdf_open("mandelbrot")
VTKHDF.grid_info(r_mandel)

# Arrays are read by indexing (image-like data keeps its full 3-D shape),
# and the marked active attributes are queryable:

size(r_mandel["Iterations"]), size(r_mandel["IterationsGradient"])

#-

VTKHDF.active_attributes(r_mandel, VTKPointData())

#-

close(r_mandel)

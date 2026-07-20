# # OverlappingAMR: Gaussian pulse
#
# Port of the `amr_gaussian_pulse.vtkhdf` example from the VTKHDF
# specification: a two-level overlapping AMR of a Gaussian pulse with the
# same box layout as the reference file. Box extents are inclusive
# cell-index ranges in each level's spacing; here the refinement ratio
# between the levels is 2.

using WriteVTKHDF

origin = (-2.0, -2.0, 0.0)
center = (-0.75, -0.75, 1.25)

# Cell centroids of a box, as a `(3, ncells)` matrix (x fastest, VTK order):

function centroids(extents, spacing)
    (i0, i1, j0, j1, k0, k1) = extents
    cells = vec(collect(Iterators.product(i0:i1, j0:j1, k0:k1)))
    return [origin[d] + spacing[d] * (c[d] + 0.5) for d in 1:3, c in cells]
end

function add_pulse_box(lvl, extents, spacing)
    c = centroids(extents, spacing)
    pulse = [2 * exp(-sum(abs2, c[:, i] .- center) / 0.5) for i in axes(c, 2)]
    return add_box(lvl, extents; celldata = ("Gaussian-Pulse" => pulse, "Centroid" => c))
end

vtkhdf_amr("gaussian_pulse"; origin) do amr
    level0 = add_level(amr; spacing = (0.5, 0.5, 0.5))
    add_pulse_box(level0, (0, 4, 0, 4, 0, 4), (0.5, 0.5, 0.5))
    level1 = add_level(amr; spacing = (0.25, 0.25, 0.25))
    add_pulse_box(level1, (0, 3, 0, 5, 0, 9), (0.25, 0.25, 0.25))
    add_pulse_box(level1, (6, 9, 4, 9, 0, 9), (0.25, 0.25, 0.25))
end

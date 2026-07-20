# # PolyData: warped torus surface
#
# Port of the `test_poly_data.vtkhdf` example from the VTKHDF specification:
# a closed polygonal surface with `Normals` and `Warping` point vectors and
# a `Materials` cell array. PolyData files always store the four cell
# categories Vertices, Lines, Polygons and Strips (in that order, which is
# also the cell-data order); only Polygons is populated here. The surface
# warped by `Warping` and colored by `Materials` in ParaView:
#
# ![Warped torus colored by material](../assets/examples/poly_data-light.png)
# ![Warped torus colored by material](../assets/examples/poly_data-dark.png)

using WriteVTKHDF

nu, nv = 48, 24
R, r = 1.0f0, 0.4f0

points = Matrix{Float32}(undef, 3, nu * nv)
normals = similar(points)
warping = similar(points)
id(iu, iv) = mod(iu, nu) + 1 + nu * mod(iv, nv)
for iv in 0:(nv - 1), iu in 0:(nu - 1)
    u, v = 2π * iu / nu, 2π * iv / nv
    normal = (cos(v) * cos(u), cos(v) * sin(u), sin(v))
    points[:, id(iu, iv)] .= (R * cos(u), R * sin(u), 0) .+ r .* normal
    normals[:, id(iu, iv)] .= normal
    warping[:, id(iu, iv)] .= 0.2 * sin(3u) .* normal
end

# One quad per grid cell, wrapping around in both directions:

quads = vec(
    [
        MeshCell(PolyData.Polys(), [id(iu, iv), id(iu + 1, iv), id(iu + 1, iv + 1), id(iu, iv + 1)])
            for iu in 0:(nu - 1), iv in 0:(nv - 1)
    ]
)
materials = vec([iu < nu ÷ 2 ? 1 : 2 for iu in 0:(nu - 1), iv in 0:(nv - 1)])

vtkhdf_grid("torus", points, quads) do vtk
    vtk["Normals", VTKPointData(), attribute = :Normals] = normals
    vtk["Warping", VTKPointData(), attribute = :Vectors] = warping
    vtk["Materials", VTKCellData()] = materials
end

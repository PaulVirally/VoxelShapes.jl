using VoxelShapes
using Plots

vacuum = 0.0
medium = 1.0
sphere = FillableSphere((1//4, 1//4, 3//4), 1//8, medium)
cube = FillableCube((2//4, 2//4, 1//4), 1//4, medium)
shapes = [sphere, cube]
num_cells = 32
scale = 1//(num_cells + 1)
world = World((num_cells, num_cells, num_cells), (scale, scale, scale), shapes, vacuum, SuperResolutionAntiAliasing())

arr = zeros(typeof(vacuum), size(world))
fill!(arr, world)

anim = @animate for i in 1:size(arr, 3)
    heatmap(arr[:, :, i], title="Filled World", xlabel="x", ylabel="y", color=:bone, aspect_ratio=:equal)
end
gif(anim, "anim.gif", fps=30)

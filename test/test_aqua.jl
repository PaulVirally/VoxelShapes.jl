using Aqua

@testset "Aqua" begin
    Aqua.test_all(VoxelShapes; ambiguities=false)
end

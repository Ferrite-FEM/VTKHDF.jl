# Test.detect_closure_boxes only exists on Julia nightly (1.13-DEV).
@testset "closure boxes" begin
    if isdefined(Test, :detect_closure_boxes)
        @test isempty(Test.detect_closure_boxes(VTKHDF))
    else
        @test true
    end
end

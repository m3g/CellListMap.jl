@testitem "performance" begin
    let cmd = `$(Base.julia_cmd()) --startup-file=no --check-bounds=auto ../modules/performance.jl`
        success(pipeline(cmd; stdout=stdout, stderr=stderr)) || error("performance check, cmd : $cmd")
    end
end
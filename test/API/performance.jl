@testitem "performance" tags=[:performance] begin
    let cmd = `$(Base.julia_cmd()[1]) --startup-file=no --check-bounds=auto ../modules/performance.jl`
        println(cmd)
        success(pipeline(cmd; stdout=stdout, stderr=stderr)) || error("performance check, cmd : $cmd")
    end
end
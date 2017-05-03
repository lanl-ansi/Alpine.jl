# Used for Windows system
# https://github.com/JuliaOpt/Ipopt.jl/issues/77
using WinRPM
pkgs = String[]
for pkg in eachline(WinRPM.installedlist)
    name = match(Regex("$(WinRPM.OS_ARCH)-(.*)"), last(split(pkg)))
    if name !== nothing && !isempty(WinRPM.lookup(name[1]))
        push!(pkgs, name[1])
    end
end
WinRPM.install(pkgs)
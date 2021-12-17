function read_problem(filename)
    caps = Array{Int64, 1}()
    dems = Array{Int64, 1}()
    ass = zeros(Int64, 0, 0)
    open(filename) do f
        tok = split(readline(f))
        ndeps = parse(Int64, tok[1])
        ncusts = parse(Int64, tok[2])
        caps = zeros(Int64, ndeps)
        dems = zeros(Int64, ncusts)
        ass = zeros(Int64, ndeps, ncusts)
        while !eof(f)
            tok = split(readline(f))
            if tok[1] == "cap"
                i = parse(Int64, tok[2])
                c = parse(Int64, tok[3])
                caps[i] = c
            elseif tok[1] == "dem"
                j = parse(Int64, tok[2])
                d = parse(Int64, tok[3])
                dems[j] = d
            elseif tok[1] == "cost"
                i = parse(Int64, tok[2])
                j = parse(Int64, tok[3])
                c = parse(Int64, tok[4])
                ass[i, j] = c
            end
        end
    end
    TranspProblemData(caps, dems, ass)
end


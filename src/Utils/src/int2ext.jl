"""
    int2ext(i2e::Vector{Int}, bus::Matrix{Float64}, gen::Matrix{Float64}, 
            branch::Matrix{Float64}, load::Matrix{Float64}, pvarray, 
            areas::Union{Matrix{Float64},Nothing}=nothing)

Convert internal bus numbering to external bus numbering for power system data.

This function maps internal bus indices to their corresponding external bus numbers
using the provided mapping vector.

# Arguments
- `i2e::Vector{Int}`: Mapping from internal to external bus numbers
- `bus::Matrix{Float64}`: Bus data matrix
- `gen::Matrix{Float64}`: Generator data matrix
- `branch::Matrix{Float64}`: Branch data matrix
- `load::Matrix{Float64}`: Load data matrix
- `pvarray`: PV array data
- `areas::Union{Matrix{Float64},Nothing}`: Area data matrix (optional)

# Returns
- Tuple of converted matrices with external bus numbering
"""
function int2ext(i2e::Vector{Int}, bus::Matrix{Float64}, gen::Matrix{Float64}, 
    branch::Matrix{Float64}, load::Matrix{Float64}, pvarray, areas::Union{Matrix{Float64},Nothing}=nothing)

# Convert bus numbers
bus[:, BUS_I] = i2e[Int.(bus[:, BUS_I])]
gen[:, GEN_BUS] = i2e[Int.(gen[:, GEN_BUS])]
branch[:, F_BUS] = i2e[Int.(branch[:, F_BUS])]
branch[:, T_BUS] = i2e[Int.(branch[:, T_BUS])]
load[:, LOAD_CND] = i2e[Int.(load[:, LOAD_CND])]
if size(pvarray, 1) > 0
    pvarray[:, PV_BUS] = i2e[Int.(pvarray[:, PV_BUS])]
end

return bus, gen, branch, load, pvarray, areas
end

function int2ext(mpc::Dict)
"""
    int2ext(mpc::Dict)

Convert a power system case from internal to external bus numbering.

This function converts all data in the MATPOWER case from internal to external
bus numbering format, restoring the original bus numbers.

# Arguments
- `mpc::Dict`: The MATPOWER case dictionary with internal numbering

# Returns
- `Dict`: The MATPOWER case with external bus numbering
"""
if !haskey(mpc, "order")
error("int2ext: mpc does not have the 'order' field required for conversion back to external numbering.")
end

o = mpc["order"]

if o["state"] == "i"
# Define constant indices
BUS_I = 1
GEN_BUS = 1
F_BUS = 1
T_BUS = 2

# Save internal numbered data and restore original data
o["int"]["bus"] = copy(mpc["bus"])
o["int"]["branch"] = copy(mpc["branch"])
o["int"]["gen"] = copy(mpc["gen"])
mpc["bus"] = copy(o["ext"]["bus"])
mpc["branch"] = copy(o["ext"]["branch"])
mpc["gen"] = copy(o["ext"]["gen"])

# If needed, pad with zeros on the right
nci = size(o["int"]["bus"], 2)
nr, nc = size(mpc["bus"])
if nc < nci
mpc["bus"] = hcat(mpc["bus"], zeros(nr, nci-nc))
end

nci = size(o["int"]["branch"], 2)
nr, nc = size(mpc["branch"])
if nc < nci
mpc["branch"] = hcat(mpc["branch"], zeros(nr, nci-nc))
end

nci = size(o["int"]["gen"], 2)
nr, nc = size(mpc["gen"])
if nc < nci
mpc["gen"] = hcat(mpc["gen"], zeros(nr, nci-nc))
end

# Update data
bus_on = o["bus"]["status"]["on"]
branch_on = o["branch"]["status"]["on"]
gen_on = o["gen"]["status"]["on"]

mpc["bus"][bus_on, :] = o["int"]["bus"]
mpc["branch"][branch_on, :] = o["int"]["branch"]
mpc["gen"][gen_on, :] = o["int"]["gen"][o["gen"]["e2i"], :]

# Restore original bus numbers
mpc["bus"][bus_on, BUS_I] = o["bus"]["i2e"][Int.(mpc["bus"][bus_on, BUS_I])]
mpc["branch"][branch_on, F_BUS] = o["bus"]["i2e"][Int.(mpc["branch"][branch_on, F_BUS])]
mpc["branch"][branch_on, T_BUS] = o["bus"]["i2e"][Int.(mpc["branch"][branch_on, T_BUS])]
mpc["gen"][gen_on, GEN_BUS] = o["bus"]["i2e"][Int.(mpc["gen"][gen_on, GEN_BUS])]

# Process additional fields
if haskey(mpc, "gencost")
if size(mpc["gencost"], 1) == 2 * size(mpc["gen"], 1) && size(mpc["gencost"], 1) != 0
    ordering = ["gen", "gen"]  # Includes Qg cost
else
    ordering = ["gen"]         # Only Pg cost
end
mpc = i2e_field(mpc, "gencost", ordering)
end

if haskey(mpc, "bus_name")
mpc = i2e_field(mpc, "bus_name", ["bus"])
end

# Update state
if haskey(o, "ext")
delete!(o, "ext")
end
o["state"] = "e"
mpc["order"] = o

else
error("int2ext: mpc claims it is already using external numbering.")
end

return mpc
end

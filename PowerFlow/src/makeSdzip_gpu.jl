mutable struct Sd_gpu
    z::CuVector{ComplexF64}
    i::CuVector{ComplexF64}
    p::CuVector{ComplexF64}
end

function makeSdzip_gpu(baseMVA, bus_gpu, pw_1, pw_2, pw_3)
    (PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM,VA, 
    BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN, PER_CONSUMER) = PowerFlow.idx_bus();
    qw_1 = pw_1
    qw_2 = pw_2
    qw_3 = pw_3
    z = (bus_gpu[:, PD] .* pw_3  + 1im * bus_gpu[:, QD] .* qw_3) / baseMVA
    i = (bus_gpu[:, PD] .* pw_2  + 1im * bus_gpu[:, QD] .* qw_2) / baseMVA
    p = (bus_gpu[:, PD] .* pw_1  + 1im * bus_gpu[:, QD] .* qw_1) / baseMVA
    sd = Sd_gpu(z[:,1], i[:,1], p[:,1])
    return sd
end
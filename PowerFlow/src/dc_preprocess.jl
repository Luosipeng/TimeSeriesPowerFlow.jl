function dc_preprocess(mpc,opt)
    mpc_list, isolated = PowerFlow.extract_islands(mpc)
    n_islands = length(mpc_list)
    if(opt["PF"]["DC_PREPROCESS"]==1)   
        preconditioned_list = Vector{Any}(undef, n_islands)
        
        @threads for i in 1:n_islands
            preconditioned_list[i] = PowerFlow.runprepf(mpc_list[i], opt)
        end
        
        # 更新相角也使用多线程
        @threads for i in 1:n_islands
            mpc_list[i]["bus"][:, 9] = preconditioned_list[i]["bus"][:, 9]
        end
    end
    return mpc_list, isolated
end
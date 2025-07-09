using SparseArrays

"""
   Definition of the Microgrid Planning Problem structure.
   This structure contains all parameters and variables for microgrid planning optimization.
"""
mutable struct MicrogridPlanningProblem
    # Constants
    bigM::Float64
    
    # Planning parameters
    nPV::Int
    nBESS::Int
    DeltaPV::Float64
    DeltaBESS::Float64
    
    # Technical parameters
    soc0::Float64
    soc_max::Float64
    soc_min::Float64
    eta_ug::Float64
    eff_ch::Float64
    eff_dc::Float64
    
    # Economic parameters
    pv_unit_cost::Float64
    ess_unit_cost_e::Float64
    ess_unit_cost_p::Float64
    ug_unit_cost::Float64
    carbon_permit_price::Float64
    carbon_emission_coeff::Float64
    npv::Float64

    # Reliability parameters
    ug_prob::Float64
    mg_prob::Float64
    power_requirements::Float64
    energy_requirements::Float64
    pv_power_factor::Float64
    pv_energy_factor::Float64
    
    # Carbon emission limit
    carbon_limit::Float64
    # Uncertainty parameters
    ns::Int
    ws::Vector{Float64}
    T::Int
    ppv_un::Int
    pl::Int
    fl::Int
    ug_status::Int
    ele_price::Int
    carbon_price::Int
    carbon_emission_ug::Int
    nu::Int
    NU::Int
    
    # Indices
    IPV::Int
    IESS::Int
    IPESS::Int
    PUG::Int
    CARBON_PERMIT::Int
    CARBON::Int
    NX::Int
    
    # First stage constraints
    lx::Vector{Float64}
    ux::Vector{Float64}
    vtypex::Vector{Char}
    Aeq::SparseMatrixCSC{Float64, Int}
    beq::Vector{Float64}
    A::SparseMatrixCSC{Float64, Int}
    b::Vector{Float64}
    
    # Objective function
    c::Vector{Float64}

    # Second stage indices
    pug::Int
    pug_sold::Int
    ppv::Int
    pess_ch::Int
    pess_dc::Int
    eess::Int
    iess_ch::Int
    iess_dc::Int
    wpv::Int
    wpess_dc::Int
    wpess_ch::Int
    plc::Int
    pf::Int
    ny::Int
    NY::Int
    
    """
    Constructor for MicrogridPlanningProblem.
    
    # Arguments
    - `data`: Dictionary containing planning parameters, technical parameters, economic parameters, and reliability parameters
    
    # Returns
    - A new instance of MicrogridPlanningProblem with initialized fields
    """
    function MicrogridPlanningProblem(data::Dict{Symbol, Any})
        # Create instance
        self = new()
    
        # Initialize constants
        self.bigM = 1e6

        # Read planning parameters
        planning_params = data[:planning_params]
        self.nPV = planning_params[1, "PV units"]
        self.nBESS = planning_params[1, "Storage units"]
        self.DeltaPV = planning_params[1, "PV unit capacity"]
        self.DeltaBESS = planning_params[1, "Storage unit capacity"]
    
        # Read technical parameters
        tech_params = data[:tech_params]
        self.soc0 = tech_params[1, "Initial SOC"]
        self.soc_max = tech_params[1, "Maximum SOC"]
        self.soc_min = tech_params[1, "Minimum SOC"]
        self.eta_ug = tech_params[1, "Grid efficiency"]
        self.eff_ch = tech_params[1, "Charging efficiency"]
        self.eff_dc = tech_params[1, "Discharging efficiency"]
    
        # Read economic parameters
        econ_params = data[:econ_params]
        self.pv_unit_cost = econ_params[1, "PV unit cost (\$/kW)"]
        self.ess_unit_cost_e = econ_params[1, "Storage energy cost (\$/kWh)"]
        self.ess_unit_cost_p = econ_params[1, "Storage power cost (\$/kW)"]
        self.ug_unit_cost = econ_params[1, "Capacity cost (\$/year/kW)"]
        self.carbon_permit_price = econ_params[1, "Carbon permit price (\$/kWh)"]
        self.carbon_emission_coeff = econ_params[1, "Carbon emission coefficient (kg/kWh)"]
        self.npv = capital_recovery_factor(econ_params[1, "Discount rate"], econ_params[1, "Planning years"])

        # Read reliability parameters
        reliability_params = data[:reliability_params]
        self.ug_prob = reliability_params[1, "Grid reliability"]
        self.mg_prob = reliability_params[1, "Microgrid reliability"]
    
        # Define uncertainty parameters (placeholders to be filled from scenario generation)
        self.ns = 0
        self.ws = Float64[]
        self.T = planning_params[1, "Time periods"]
        self.ppv_un = 0  # PV output in per unit
        self.pl = self.ppv_un + 1  # Active load consumption
        self.fl = self.pl + 1  # Flexible load consumption
        self.ug_status = self.fl + 1  # Utility grid operating status
        self.ele_price = self.ug_status + 1  # Electricity price
        self.carbon_price = self.ele_price + 1  # Carbon price
        self.carbon_emission_ug = self.carbon_price + 1  # Carbon price
        self.nu = self.carbon_emission_ug + 1 # Uncertainties of loads
        self.NU = self.nu * self.T  # Total uncertainties
    
        # Return instance
        return self
    end
end

"""
    Definition of the power flow problem structure.
    Aligned with JPC structure for consistency.
"""
mutable struct JuliaPowerCase
    version::String
    baseMVA::Float32
    basef::Float32
    
    # AC Network Components
    busesAC::Vector{Bus}                # Matches busAC in JPC
    branchesAC::Vector{Line}            # Matches branchAC in JPC
    loadsAC::Vector{Load}               # Matches loadAC in JPC
    loadsAC_flex::Vector{FlexLoad}      # Matches loadAC_flex in JPC
    loadsAC_asymm::Vector{AsymmetricLoad}  # Matches loadAC_asymm in JPC
    # branches3ph::Vector{ThreePhaseBranch}  # Matches branch3ph in JPC
    gensAC::Vector{Generator}            # Matches genAC in JPC
    
    # DC Network Components
    busesDC::Vector{BusDC}              # Matches busDC in JPC
    branchesDC::Vector{LineDC}          # Matches branchDC in JPC
    loadsDC::Vector{LoadDC}          # Matches loadDC in JPC
    
    # Distributed Energy Resources
    sgensAC::Vector{StaticGenerator}     # Matches sgenAC in JPC
    storages::Vector{Storage}            # Matches storage in JPC
    storageetap::Vector{Storageetap}      # Matches storageetap in JPC
    sgensDC::Vector{StaticGeneratorDC}   # Matches sgenDC in JPC
    pvarray::Vector{PVArray}          # Matches pv_array in JPC
    ACPVSystems::Vector{ACPVSystem}  # Matches AC PV systems in JPC
    
    # Special Components
    converters::Vector{Converter}        # Matches converter in JPC
    ext_grids::Vector{ExternalGrid}      # Matches ext_grid in JPC
    hvcbs::Vector{HighVoltageCircuitBreaker}  # Matches hvcb in JPC
    microgrids::Vector{Microgrid}        # Matches microgrid in JPC
    
    # Original components without direct JPC equivalents (preserved)
    transformers_2w::Vector{Transformer2W}
    transformers_3w::Vector{Transformer3W}
    transformers_2w_etap::Vector{Transformer2Wetap}
    charging_stations::Vector{ChargingStation}
    chargers::Vector{Charger}
    ev_aggregators::Vector{EVAggregator}
    virtual_power_plants::Vector{VirtualPowerPlant}
    carbon_time_series::Vector{CarbonTimeSeries}
    equipment_carbon::Vector{EquipmentCarbon}
    
    # Lookup dictionaries
    bus_name_to_id::Dict{String, Int}
    busdc_name_to_id::Dict{String, Int}
    zone_to_id::Dict{String, Int}
    area_to_id::Dict{String, Int}
    
    """
    Constructor for JuliaPowerCase.
    
    # Arguments
    - `version`: Version string (default: "2.0")
    - `baseMVA`: Base MVA for the system (default: 100.0)
    - `basef`: Base frequency in Hz (default: 50.0)
    
    # Returns
    - A new instance of JuliaPowerCase with empty component vectors
    """
    function JuliaPowerCase(version::String = "2.0", baseMVA::Float32 = 100.0f0, basef::Float32 = 50.0f0)
        return new(
            version, baseMVA, basef,
            # AC Network Components
            Bus[], Line[], Load[], FlexLoad[], AsymmetricLoad[], Generator[],
            # DC Network Components
            BusDC[], LineDC[],LoadDC[],
            # Distributed Energy Resources
            StaticGenerator[], Storage[], Storageetap[], StaticGeneratorDC[], PVArray[], ACPVSystem[],
            # Special Components
            Converter[], ExternalGrid[], HighVoltageCircuitBreaker[], Microgrid[],
            # Original components without direct JPC equivalents
            Transformer2W[], Transformer3W[],Transformer2Wetap[], ChargingStation[], Charger[], 
            EVAggregator[], VirtualPowerPlant[], CarbonTimeSeries[], EquipmentCarbon[],
            # Lookup dictionaries
            Dict{String, Int}(), Dict{String, Int}(), Dict{String, Int}(), Dict{String, Int}(),
        )
    end
end

"""
    Definition of the JPC structure as one matrix format.
    This structure represents a power system case in a one matrix format,
    including buses, generators, branches, loads, and other power system components.
    Each matrix is sized according to the number of attributes defined in idx.jl.
"""
mutable struct JPC
    version::String
    baseMVA::Float64
    success::Bool
    iterationsAC::Int
    iterationsDC::Int
    
    # AC Network Components
    busAC::Array{Float64,2}         # Bus data (indexed by idx_bus(), dimensions: n_bus × N_BUS_ATTR)
    genAC::Array{Float64,2}         # Generator data (indexed by idx_gen(), dimensions: n_gen × N_GEN_ATTR)
    branchAC::Array{Float64,2}      # Branch data (indexed by idx_brch(), dimensions: n_branch × N_BRANCH_ATTR)
    loadAC::Array{Float64,2}        # Load data (indexed by idx_ld(), dimensions: n_load × N_LOAD_ATTR)
    loadAC_flex::Array{Float64,2}   # Flexible load data (indexed by idx_ld_flex(), dimensions: n_flex_load × N_FLEX_LOAD_ATTR)
    loadAC_asymm::Array{Float64,2}  # Asymmetric load data (indexed by idx_ld_asymmetric(), dimensions: n_asym_load × N_ASYM_LOAD_ATTR)
    branch3ph::Array{Float64,2}     # Three-phase branch data (indexed by idx_3ph_brch(), dimensions: n_3ph_branch × N_3PH_BRANCH_ATTR)
    
    # DC Network Components
    busDC::Array{Float64,2}         # DC bus data (indexed by idx_dcbus(), dimensions: n_dcbus × N_DCBUS_ATTR)
    branchDC::Array{Float64,2}      # DC branch data (indexed by idx_brch(), dimensions: n_dc_branch × N_BRANCH_ATTR)
    genDC::Array{Float64,2}         # DC generator data (indexed by idx_gen(), dimensions: n_dc_gen × N_GEN_ATTR)
    loadDC::Array{Float64,2}        # DC load data (indexed by idx_ld(), dimensions: n_dc_load × N_LOAD_ATTR)
    
    # Distributed Energy Resources
    sgenAC::Array{Float64,2}        # Solar/PV generation data (indexed by idx_pv(), dimensions: n_pv × N_PV_ATTR)
    storageetap::Array{Float64,2}  # Energy storage system data (indexed by idx_storageetap(), dimensions: n_storageetap × N_STORAGEETAP_ATTR)
    storage::Array{Float64,2}       # Energy storage system data (indexed by idx_ess(), dimensions: n_ess × N_ESS_ATTR)
    sgenDC::Array{Float64,2}        # DC-connected PV data (indexed by idx_pv(), dimensions: n_dc_pv × N_PV_ATTR)
    pv::Array{Float64,2}          # PV array data (indexed by idx_pv_array(), dimensions: n_pv_array × N_PV_ARRAY_ATTR)
    pv_acsystem::Array{Float64,2}  # AC PV systems data (indexed by idx_ac_pv_system(), dimensions: n_ac_pv_system × N_AC_PV_SYSTEM_ATTR)
    
    # Special Components
    converter::Array{Float64,2}     # AC/DC converter data (indexed by idx_conv(), dimensions: n_conv × N_CONV_ATTR)
    ext_grid::Array{Float64,2}      # External grid data (indexed by idx_ext_grid(), dimensions: n_ext_grid × N_EXT_GRID_ATTR)
    hvcb::Array{Float64,2}          # High voltage circuit breaker data (indexed by idx_hvcb(), dimensions: n_hvcb × N_HVCB_ATTR)
    microgrid::Array{Float64,2}     # Microgrid data (indexed by idx_microgrid(), dimensions: n_mg × N_MG_ATTR)
    
    """
    Constructor for JPC with default empty arrays.
    
    # Arguments
    - `version`: Version string (default: "2.0")
    - `baseMVA`: Base MVA for the system (default: 100.0)
    - `success`: Boolean indicating if power flow converged (default: false)
    - `iterationsAC`: Number of AC power flow iterations (default: 0)
    - `iterationsDC`: Number of DC power flow iterations (default: 0)
    
    # Returns
    - A new instance of JPC with empty component arrays
    """
    function JPC(version::String = "2.0", baseMVA::Float64 = 100.0, success::Bool = false, iterationsAC::Int = 0, iterationsDC::Int = 0)
        # Initialize with empty arrays - they will be properly sized when data is loaded
        new(version, baseMVA, success, iterationsAC, iterationsDC,
            Array{Float64}(undef, 0, 22),  # busAC (N_BUS_ATTR = 22)
            Array{Float64}(undef, 0, 32),  # genAC (N_GEN_ATTR = 32)
            Array{Float64}(undef, 0, 28),  # branchAC (N_BRANCH_ATTR = 28)
            Array{Float64}(undef, 0, 8),   # loadAC (N_LOAD_ATTR = 8)
            Array{Float64}(undef, 0, 25),  # loadAC_flex (N_FLEX_LOAD_ATTR = 25)
            Array{Float64}(undef, 0, 12),  # loadAC_asymm (N_ASYM_LOAD_ATTR = 12)
            Array{Float64}(undef, 0, 30),  # branch3ph (N_3PH_BRANCH_ATTR = 30)
            Array{Float64}(undef, 0, 21),  # busDC (N_DCBUS_ATTR = 21)
            Array{Float64}(undef, 0, 28),  # branchDC (uses same structure as branchAC)
            Array{Float64}(undef, 0, 32),  # genDC (N_GEN_ATTR = 32)
            Array{Float64}(undef, 0, 8),   # loadDC (N_LOAD_ATTR = 8)
            Array{Float64}(undef, 0, 3),   # sgenAC (N_PV_ATTR = 3)
            Array{Float64}(undef, 0, 15),  # storageetap (N_STORAGEETAP_ATTR = 15)
            Array{Float64}(undef, 0, 15),  # storage (N_ESS_ATTR = 15)
            Array{Float64}(undef, 0, 3),   # sgenDC (uses same structure as sgenAC)
            Array{Float64}(undef, 0, 9),  # pv (N_PV_ARRAY_ATTR = 10)
            Array{Float64}(undef, 0, 15),  # pv_acsystem (N_AC_PV_SYSTEM_ATTR = 18)
            Array{Float64}(undef, 0, 18),  # converter (N_CONV_ATTR = 18)
            Array{Float64}(undef, 0, 13),  # ext_grid (N_EXT_GRID_ATTR = 13)
            Array{Float64}(undef, 0, 5),   # hvcb (N_HVCB_ATTR = 5)
            Array{Float64}(undef, 0, 5)    # microgrid (N_MG_ATTR = 5)
        )
    end
end

import Base: getindex, setindex!

"""
    getindex(jpc::JPC, key::String)

Access JPC components by string key.

# Arguments
- `jpc`: JPC structure
- `key`: String key representing the component to access

# Returns
- The requested component data
"""
function getindex(jpc::JPC, key::String)
    # Return the field corresponding to the string key
    if key == "version"
        return jpc.version
    elseif key == "baseMVA"
        return jpc.baseMVA
    elseif key == "success"
        return jpc.success
    elseif key == "iterationsAC"
        return jpc.iterationsAC
    elseif key == "iterationsDC"
        return jpc.iterationsDC
    elseif key == "busAC"
        return jpc.busAC
    elseif key == "genAC"
        return jpc.genAC
    elseif key == "branchAC"
        return jpc.branchAC
    elseif key == "loadAC"
        return jpc.loadAC
    elseif key == "loadAC_flex"
        return jpc.loadAC_flex
    elseif key == "loadAC_asymm"
        return jpc.loadAC_asymm
    elseif key == "branch3ph"
        return jpc.branch3ph
    elseif key == "busDC"
        return jpc.busDC
    elseif key == "branchDC"
        return jpc.branchDC
    elseif key == "genDC"
        return jpc.genDC
    elseif key == "loadDC"
        return jpc.loadDC
    elseif key == "sgenAC"
        return jpc.sgenAC
    elseif key == "storageetap"
        return jpc.storageetap
    elseif key == "storage"
        return jpc.storage
    elseif key == "sgenDC"
        return jpc.sgenDC
    elseif key == "pv"
        return jpc.pv
    elseif key == "pv_acsystem"
        return jpc.pv_acsystem
    elseif key == "converter"
        return jpc.converter
    elseif key == "ext_grid"
        return jpc.ext_grid
    elseif key == "hvcb"
        return jpc.hvcb
    elseif key == "microgrid"
        return jpc.microgrid
    else
        error("Key does not exist in JPC structure: $key")
    end
end

"""
    setindex!(jpc::JPC, value, key::String)

Set JPC components by string key.

# Arguments
- `jpc`: JPC structure
- `value`: Value to assign
- `key`: String key representing the component to modify
"""
function setindex!(jpc::JPC, value, key::String)
    # Set the field corresponding to the string key
    if key == "version"
        jpc.version = value
    elseif key == "baseMVA"
        jpc.baseMVA = value
    elseif key == "success"
        jpc.success = value
    elseif key == "iterationsAC"
        jpc.iterationsAC = value
    elseif key == "iterationsDC"
        jpc.iterationsDC = value
    elseif key == "busAC"
        jpc.busAC = value
    elseif key == "genAC"
        jpc.genAC = value
    elseif key == "branchAC"
        jpc.branchAC = value
    elseif key == "loadAC"
        jpc.loadAC = value
    elseif key == "loadAC_flex"
        jpc.loadAC_flex = value
    elseif key == "loadAC_asymm"
        jpc.loadAC_asymm = value
    elseif key == "branch3ph"
        jpc.branch3ph = value
    elseif key == "busDC"
        jpc.busDC = value
    elseif key == "branchDC"
        jpc.branchDC = value
    elseif key == "genDC"
        jpc.genDC = value
    elseif key == "loadDC"
        jpc.loadDC = value
    elseif key == "sgenAC"
        jpc.sgenAC = value
    elseif key == "storageetap"
        jpc.storageetap = value
    elseif key == "storage"
        jpc.storage = value
    elseif key == "sgenDC"
        jpc.sgenDC = value
    elseif key == "pv"
        jpc.pv = value
    elseif key == "pv_acsystem"
        jpc.pv_acsystem = value
    elseif key == "converter"
        jpc.converter = value
    elseif key == "ext_grid"
        jpc.ext_grid = value
    elseif key == "hvcb"
        jpc.hvcb = value
    elseif key == "microgrid"
        jpc.microgrid = value
    else
        error("Key does not exist in JPC structure: $key")
    end
end

"""
    Definition of the JPC_3ph structure for three-phase power flow analysis.
    This structure extends the JPC format to support three-phase power flow calculations.
"""
mutable struct JPC_3ph
    version::String
    baseMVA::Float32
    basef::Float32
    mode::String
    success::Bool  # Flag indicating if three-phase power flow converged
    iterations::Int  # Number of iterations for three-phase power flow
    
    # AC Network Components - Matrix representations of JuliaPowerCase components
    busAC_0::Array{Float64,2}  # Phase A bus data
    busAC_1::Array{Float64,2}  # Phase B bus data
    busAC_2::Array{Float64,2}  # Phase C bus data

    branchAC_0::Array{Float64,2}  # Phase A branch data
    branchAC_1::Array{Float64,2}  # Phase B branch data
    branchAC_2::Array{Float64,2}  # Phase C branch data

    loadAC_0::Array{Float64,2}  # Phase A load data
    loadAC_1::Array{Float64,2}  # Phase B load data
    loadAC_2::Array{Float64,2}  # Phase C load data

    genAC_0::Array{Float64,2}  # Phase A generator data
    genAC_1::Array{Float64,2}  # Phase B generator data
    genAC_2::Array{Float64,2}  # Phase C generator data

    storageAC::Array{Float64,2}  # Storage data
    
    # DC Network Components
    busDC::Array{Float64,2}  # DC bus data
    branchDC::Array{Float64,2}  # DC branch data
    loadDC::Array{Float64,2}  # DC load data
    genDC::Array{Float64,2}  # DC generator data
    storageDC::Array{Float64,2}  # DC storage data
    
    # Special Components
    ext_grid::Array{Float64,2}  # External grid data
    switche::Array{Float64,2}  # Switch data

    # Three-phase power flow results
    res_bus_3ph::Array{Float64,2}  # Bus results from three-phase power flow
    res_loadsAC_3ph::Array{Float64,2}  # Load results from three-phase power flow
    res_ext_grid_3ph::Array{Float64,2}  # External grid results from three-phase power flow
    
    # Lookup dictionaries
    bus_name_to_id::Dict{String, Int}
    zone_to_id::Dict{String, Int}
    area_to_id::Dict{String, Int}
    
    """
    Constructor for JPC_3ph.
    
    # Arguments
    - `version`: Version string (default: "2.0")
    - `baseMVA`: Base MVA for the system (default: 100.0)
    - `basef`: Base frequency in Hz (default: 50.0)
    - `mode`: Mode of operation (default: "etap")
    - `success`: Boolean indicating if power flow converged (default: false)
    - `iterations`: Number of power flow iterations (default: 0)
    
    # Returns
    - A new instance of JPC_3ph with empty component arrays
    """
    function JPC_3ph(version::String = "2.0", baseMVA::Float32 = 100.0f0, basef::Float32 = 50.0f0, mode::String = "etap",
                     success::Bool = false, iterations::Int = 0)
        # Initialize with empty matrices - they will be properly sized when data is loaded
        return new(
            version, baseMVA, basef, mode, success, iterations,
            # Initialize empty matrices - column counts match attribute counts
            Array{Float64}(undef, 0, 22),  # busAC_0 (N_BUS_ATTR = 22)
            Array{Float64}(undef, 0, 22),  # busAC_1 (N_BUS_ATTR = 22)
            Array{Float64}(undef, 0, 22),  # busAC_2 (N_BUS_ATTR = 22)

            Array{Float64}(undef, 0, 28),  # branchAC_0 (N_BRANCH_ATTR = 28)
            Array{Float64}(undef, 0, 28),  # branchAC_1 (N_BRANCH_ATTR = 28)
            Array{Float64}(undef, 0, 28),  # branchAC_2 (N_BRANCH_ATTR = 28)
            Array{Float64}(undef, 0, 8),   # loadAC_0 (N_LOAD_ATTR = 8)
            Array{Float64}(undef, 0, 8),   # loadAC_1 (N_LOAD_ATTR = 8)
            Array{Float64}(undef, 0, 8),   # loadAC_2 (N_LOAD_ATTR = 8)

            Array{Float64}(undef, 0, 32),  # genAC_0 (N_GEN_ATTR = 32)
            Array{Float64}(undef, 0, 32),  # genAC_1 (N_GEN_ATTR = 32)
            Array{Float64}(undef, 0, 32),  # genAC_2 (N_GEN_ATTR = 32)

            Array{Float64}(undef, 0, 15),  # storageAC (N_ESS_ATTR = 15)

            Array{Float64}(undef, 0, 21),  # busDC (N_DCBUS_ATTR = 21)
            Array{Float64}(undef, 0, 28),  # branchDC (uses same structure as branchAC)
            Array{Float64}(undef, 0, 8),   # loadDC (N_LOAD_ATTR = 8)
            Array{Float64}(undef, 0, 32),  # genDC (N_GEN_ATTR = 32)

            Array{Float64}(undef, 0, 15),  # storageDC (N_ESS_ATTR = 15)
            Array{Float64}(undef, 0, 13),  # ext_grid (N_EXT_GRID_ATTR = 13)
            Array{Float64}(undef, 0, 5),   # hvcb (N_HVCB_ATTR = 5)

                        Array{Float64}(undef, 0, 15),   # res_bus_3ph (N_RES_BUS_3PH_ATTR = 15)
            Array{Float64}(undef, 0, 18),  # res_loadsAC_3ph (N_RES_LOADS_AC_3PH_ATTR = 18)
            Array{Float64}(undef, 0, 13),  # res_ext_grid_3ph (N_RES_EXT_GRID_3PH_ATTR = 13)

            Dict{String, Int}(),  # bus_name_to_id
            Dict{String, Int}(),  # zone_to_id
            Dict{String, Int}(),  # area_to_id
        )
    end
end

"""
    getindex(jpc::JPC_3ph, key::String)

Access JPC_3ph components by string key.

# Arguments
- `jpc`: JPC_3ph structure
- `key`: String key representing the component to access

# Returns
- The requested component data
"""
function getindex(jpc::JPC_3ph, key::String)
    # Return the field corresponding to the string key
    if key == "version"
        return jpc.version
    elseif key == "baseMVA"
        return jpc.baseMVA
    elseif key == "basef"
        return jpc.basef
    elseif key == "mode"
        return jpc.mode
    elseif key == "success"
        return jpc.success
    elseif key == "iterations"
        return jpc.iterations
    elseif key == "busAC_0"
        return jpc.busAC_0
    elseif key == "busAC_1"
        return jpc.busAC_1
    elseif key == "busAC_2"
        return jpc.busAC_2
    elseif key == "branchAC_0"
        return jpc.branchAC_0
    elseif key == "branchAC_1"
        return jpc.branchAC_1
    elseif key == "branchAC_2"
        return jpc.branchAC_2
    elseif key == "loadAC_0"
        return jpc.loadAC_0
    elseif key == "loadAC_1"
        return jpc.loadAC_1
    elseif key == "loadAC_2"
        return jpc.loadAC_2
    elseif key == "genAC_0"
        return jpc.genAC_0
    elseif key == "genAC_1"
        return jpc.genAC_1
    elseif key == "genAC_2"
        return jpc.genAC_2
    elseif key == "storageAC"
        return jpc.storageAC
    elseif key == "busDC"
        return jpc.busDC
    elseif key == "branchDC"
        return jpc.branchDC
    elseif key == "loadDC"
        return jpc.loadDC
    elseif key == "genDC"
        return jpc.genDC
    elseif key == "storageDC"
        return jpc.storageDC
    elseif key == "ext_grids"
        return jpc.ext_grids
    elseif key == "switches"
        return jpc.switches
    elseif key == "res_bus_3ph"
        return jpc.res_bus_3ph
    elseif key == "res_loadsAC_3ph"
        return jpc.res_loadsAC_3ph
    elseif key == "res_ext_grid_3ph"
        return jpc.res_ext_grid_3ph
    else
        error("Key does not exist in JPC_3ph structure: $key")
    end
end

"""
    setindex!(jpc::JPC_3ph, value, key::String)

Set JPC_3ph components by string key.

# Arguments
- `jpc`: JPC_3ph structure
- `value`: Value to assign
- `key`: String key representing the component to modify
"""
function setindex!(jpc::JPC_3ph, value, key::String)
    # Set the field corresponding to the string key
    if key == "version"
        jpc.version = value
    elseif key == "baseMVA"
        jpc.baseMVA = value
    elseif key == "basef"
        jpc.basef = value
    elseif key == "mode"
        jpc.mode = value
    elseif key == "success"
        jpc.success = value
    elseif key == "iterations"
        jpc.iterations = value
    elseif key == "busAC_0"
        jpc.busAC_0 = value
    elseif key == "busAC_1"
        jpc.busAC_1 = value
    elseif key == "busAC_2"
        jpc.busAC_2 = value
    elseif key == "branchAC_0"
        jpc.branchAC_0 = value
    elseif key == "branchAC_1"
        jpc.branchAC_1 = value
    elseif key == "branchAC_2"
        jpc.branchAC_2 = value
    elseif key == "loadAC_0"
        jpc.loadAC_0 = value
    elseif key == "loadAC_1"
        jpc.loadAC_1 = value
    elseif key == "loadAC_2"
        jpc.loadAC_2 = value
    elseif key == "genAC_0"
        jpc.genAC_0 = value
    elseif key == "genAC_1"
        jpc.genAC_1 = value
    elseif key == "genAC_2"
        jpc.genAC_2 = value
    elseif key == "storageAC"
        jpc.storageAC = value
    elseif key == "busDC"
        jpc.busDC = value
    elseif key == "branchDC"
        jpc.branchDC = value
    elseif key == "loadDC"
        jpc.loadDC = value
    elseif key == "genDC"
        jpc.genDC = value
    elseif key == "storageDC"
        jpc.storageDC = value
    elseif key == "ext_grid"
        jpc.ext_grid = value
    elseif key == "switches"
        jpc.switches = value
    elseif key == "res_bus_3ph"
        jpc.res_bus_3ph = value
    elseif key == "res_loadsAC_3ph"
        jpc.res_loadsAC_3ph = value
    elseif key == "res_ext_grid_3ph"
        jpc.res_ext_grid_3ph = value
    else
        error("Key does not exist in JPC_3ph structure: $key")
    end
end


# ComponentModel Module

The ComponentModel module provides a comprehensive framework for modeling and simulating power systems, including AC, DC, and hybrid AC/DC components, as well as carbon emission models. This document outlines the key components and their usage within the module.

## AC Power System Components

AC components serve as the fundamental building blocks for constructing power system models, including buses, lines, transformers, generators, and loads.

### Bus Components

Buses are nodes in power systems that connect various devices such as generators, loads, lines, and transformers.

```julia
# Create a new bus
bus1 = Bus(
    index = 1,
    name = "Bus 1",
    vn_kv = 110.0,  # Nominal voltage (kV)
    type = "b",     # Slack bus
    zone = 1,
    v_pu = 1.0,     # Initial voltage (per unit)
    theta_degree = 0.0  # Initial angle (degrees)
)
```

### Line Components

Line components model AC transmission and distribution lines through impedance parameters.

```julia
# Create a new line
line1 = Line(
    index = 1,
    name = "Line 1-2",
    from_bus = 1,
    to_bus = 2,
    length_km = 10.0,
    r_ohm_per_km = 0.1,  # Resistance per km (Ω/km)
    x_ohm_per_km = 0.3,  # Reactance per km (Ω/km)
    c_nf_per_km = 10.0,  # Capacitance per km (nF/km)
    max_i_ka = 0.5,      # Maximum current capacity (kA)
    df = 1.0,            # Derating factor
    parallel = 1         # Number of parallel lines
)
```

### Transformer Components

Transformer components model voltage conversion devices connecting buses at different voltage levels.

```julia
# Create a two-winding transformer
transformer = Transformer2W(
    index = 1,
    name = "T1",
    std_type = "160 MVA 110/20 kV",
    hv_bus = 1,       # High voltage bus
    lv_bus = 2,       # Low voltage bus
    sn_mva = 160.0,   # Rated power (MVA)
    vn_hv_kv = 110.0, # HV rated voltage (kV)
    vn_lv_kv = 20.0,  # LV rated voltage (kV)
    vk_percent = 10.0,# Short-circuit impedance (%)
    vkr_percent = 0.3,# Short-circuit losses (%)
    pfe_kw = 60.0,    # Iron losses (kW)
    i0_percent = 0.1, # No-load current (%)
    shift_degree = 0.0,# Phase shift angle (degrees)
    tap_side = "hv",  # Tap side (high voltage)
    tap_pos = 0,      # Tap position
    tap_neutral = 0,  # Neutral tap position
    tap_min = -10,    # Minimum tap position
    tap_max = 10      # Maximum tap position
)
```

### Generator Components

Static generator components provide simplified modeling of generating units.

```julia
# Create a static generator
gen = StaticGenerator(
    index = 1,
    name = "Generator 1",
    bus = 1,           # Connected bus number
    p_mw = 100.0,      # Active power output (MW)
    q_mvar = 20.0,     # Reactive power output (MVar)
    scaling = 1.0,     # Power scaling factor
    max_p_mw = 150.0,  # Maximum active power limit (MW)
    min_p_mw = 50.0,   # Minimum active power limit (MW)
    max_q_mvar = 50.0, # Maximum reactive power limit (MVar)
    min_q_mvar = -50.0 # Minimum reactive power limit (MVar)
)
```

### Load Components

Load components model power consumption in electrical systems.

```julia
# Create a load
load = Load(
    index = 1,
    name = "Load 1",
    bus = 2,          # Connected bus number
    p_mw = 50.0,      # Active power demand (MW)
    q_mvar = 10.0,    # Reactive power demand (MVar)
    scaling = 1.0,    # Power scaling factor
    type = "wye"      # Load type (Y-connected)
)
```

### Network Construction Example

```julia
using ComponentModel
using Utils

# Create power system model
case = JuliaPowerCase()

# Add buses
bus1 = Bus(1, "Bus 1", 110.0, "b", 1, 1.0, 0.0)
bus2 = Bus(2, "Bus 2", 110.0, "pq", 1, 1.0, 0.0)
bus3 = Bus(3, "Bus 3", 20.0, "pq", 1, 1.0, 0.0)
push!(case.busesAC, bus1)
push!(case.busesAC, bus2)
push!(case.busesAC, bus3)

# Add generator
gen = StaticGenerator(1, "Gen 1", 1, 100.0, 20.0, 1.0, 150.0, 50.0, 50.0, -50.0)
push!(case.gensAC, gen)

# Add line
line = Line(1, "Line 1-2", 1, 2, 10.0, 0.1, 0.3, 10.0, 0.5, 1.0, 1)
push!(case.branchesAC, line)

# Add transformer
transformer = Transformer2W(1, "T1", "160 MVA 110/20 kV", 2, 3, 160.0, 110.0, 20.0, 10.0, 0.3, 60.0, 0.1, 0.0, "hv", 0, 0, -10, 10)
push!(case.transformers_2w, transformer)

# Add load
load = Load(1, "Load 1", 3, 50.0, 10.0, 1.0, "wye")
push!(case.loadsAC, load)

```

## Converters

Converters model power electronic interfaces between AC and DC systems in power networks.

### Converter Structure

```julia
# Create a new converter
converter = Converter(
    index = 1,
    name = "VSC1",
    bus_ac = 5,         # Connected AC bus number
    bus_dc = 101,       # Connected DC bus number
    p_mw = 100.0,       # Active power transfer (MW)
    q_mvar = 25.0,      # Reactive power at AC side (MVar)
    vm_ac_pu = 1.0,     # AC voltage magnitude (per unit)
    vm_dc_pu = 1.05,    # DC voltage magnitude (per unit)
    loss_percent = 1.5, # Losses as percentage of active power
    loss_mw = 1.5,      # Losses (MW)
    max_p_mw = 200.0,   # Maximum active power limit (MW)
    min_p_mw = -200.0,  # Minimum active power limit (MW)
    max_q_mvar = 100.0, # Maximum reactive power limit (MVar)
    min_q_mvar = -100.0,# Minimum reactive power limit (MVar)
    control_mode = "P-Q", # Control mode
    droop_kv = 0.0,     # DC voltage droop parameter
    in_service = true,  # Operational status
    controllable = true # Controllability status
)
```

### Control Modes

Converters can operate in different control modes:

1. **P-Q Mode**: The converter controls active and reactive power
2. **Vdc-Q Mode**: The converter controls DC voltage and reactive power
3. **P-Vac Mode**: The converter controls active power and AC voltage
4. **Vdc-Vac Mode**: The converter controls DC voltage and AC voltage

### Integration Example

```julia
using ComponentModel
using Utils

# Create power system model
case = JuliaPowerCase()

# Add AC and DC buses
ac_bus = Bus(1, "AC Bus", 400.0, "pq", 1, 1.0, 0.0)
dc_bus = DCBus(101, "DC Bus", 500.0, 1.0)
push!(case.busesAC, ac_bus)
push!(case.busesDC, dc_bus)

# Add converter
converter = Converter(
    1, "VSC1", 1, 101, 100.0, 25.0, 1.0, 1.05, 
    1.5, 1.5, 200.0, -200.0, 100.0, -100.0, 
    "P-Q", 0.0, true, true
)
push!(case.converters, converter)

```

## DC Components

DC components serve as building blocks for constructing DC power system models.

### DC Bus Components

```julia
# Create a new DC bus
dc_bus = BusDC(
    index = 101,
    name = "DC Bus 1",
    vn_kv = 500.0,  # Nominal voltage (kV)
    type = "b",     # Slack bus
    zone = 1,
    v_pu = 1.0,     # Initial voltage (per unit)
    in_service = true
)
```

### DC Line Components

```julia
# Create a new DC line
dc_line = LineDC(
    index = 1,
    name = "DC Line 1-2",
    from_bus = 101,
    to_bus = 102,
    length_km = 100.0,
    r_ohm_per_km = 0.01,  # Resistance per km (Ω/km)
    max_i_ka = 2.0,       # Maximum current capacity (kA)
    in_service = true
)
```

### DC Generator Components

```julia
# Create a DC static generator
dc_gen = StaticGeneratorDC(
    index = 1,
    name = "DC Generator 1",
    bus = 101,         # Connected DC bus number
    p_mw = 200.0,      # Active power output (MW)
    scaling = 1.0,     # Power scaling factor
    max_p_mw = 300.0,  # Maximum active power limit (MW)
    min_p_mw = 0.0,    # Minimum active power limit (MW)
    in_service = true
)
```

### DC Load Components

```julia
# Create a DC load
dc_load = LoadDC(
    index = 1,
    name = "DC Load 1",
    bus = 102,          # Connected DC bus number
    p_mw = 150.0,       # Active power demand (MW)
    const_z_percent = 0.0,  # Percentage of constant impedance load
    const_i_percent = 0.0,  # Percentage of constant current load
    const_p_percent = 100.0, # Percentage of constant power load
    scaling = 1.0,      # Power scaling factor
    in_service = true
)
```

### DC Storage Components

```julia
# Create a DC storage system
dc_storage = Storage(
    index = 1,
    name = "Battery Storage 1",
    bus = 101,                # Connected DC bus number
    power_capacity_mw = 50.0, # Maximum power capacity (MW)
    energy_capacity_mwh = 200.0, # Energy storage capacity (MWh)
    soc_init = 0.5,           # Initial state of charge (50%)
    min_soc = 0.1,            # Minimum allowed state of charge (10%)
    max_soc = 0.9,            # Maximum allowed state of charge (90%)
    efficiency = 0.95,        # Round-trip efficiency
    self_discharge = 0.001,   # Self-discharge rate per hour
    p_mw = 0.0,               # Initial active power (MW, positive for charging)
    in_service = true
)
```

### DC Network Construction Example

```julia
using ComponentModel
using Utils

# Create DC power system model
case = JuliaPowerCase()

# Add DC buses
bus1 = BusDC(101, "DC Bus 1", 500.0, "b", 1, 1.0, true)
bus2 = BusDC(102, "DC Bus 2", 500.0, "pq", 1, 1.0, true)
push!(case.busesDC, bus1)
push!(case.busesDC, bus2)

# Add DC generator
gen = StaticGeneratorDC(1, "DC Gen 1", 101, 200.0, 1.0, 300.0, 0.0, true)
push!(case.sgensDC, gen)

# Add DC line
line = LineDC(1, "DC Line 1-2", 101, 102, 100.0, 0.01, 2.0, true)
push!(case.branchesDC, line)

# Add DC load
load = LoadDC(1, "DC Load 1", 102, 150.0, 0.0, 0.0, 100.0, 1.0, true)
push!(case.loadsDC, load)

# Add DC storage
storage = Storage(1, "Battery 1", 101, 50.0, 200.0, 0.5, 0.1, 0.9, 0.95, 0.001, 0.0, true)
push!(case.storages, storage)

```

## Hybrid AC/DC Systems

DC components can be integrated with AC components to form hybrid AC/DC power systems through converter components.

```julia
using ComponentModel
using Utils

# Create hybrid power system model
case = JuliaPowerCase()

# Add AC buses
ac_bus1 = Bus(1, "AC Bus 1", 400.0, "b", 1, 1.0, 0.0)
ac_bus2 = Bus(2, "AC Bus 2", 400.0, "pq", 1, 1.0, 0.0)
push!(case.busesAC, ac_bus1)
push!(case.busesAC, ac_bus2)

# Add DC buses
dc_bus1 = BusDC(101, "DC Bus 1", 500.0, "b", 1, 1.0, true)
dc_bus2 = BusDC(102, "DC Bus 2", 500.0, "pq", 1, 1.0, true)
push!(case.busesDC, dc_bus1)
push!(case.busesDC, dc_bus2)

# Add AC generator
ac_gen = StaticGenerator(1, "AC Gen", 1, 300.0, 100.0, 1.0, 500.0, 0.0, 200.0, -200.0)
push!(case.sgensAC, ac_gen)

# Add DC generator
dc_gen = StaticGeneratorDC(1, "DC Gen", 102, 100.0, 1.0, 200.0, 0.0, true)
push!(case.sgensDC, dc_gen)

# Add converter between AC and DC
converter = Converter(1, "VSC", 2, 101, 150.0, 50.0, 1.0, 1.0, 1.0, 1.5, 200.0, -200.0, 100.0, -100.0, "P-Q", 0.0, true, true)
push!(case.converters, converter)


```

## Carbon Emission Models

Carbon emission modeling components enable users to incorporate carbon emissions data into power system analyses.

### Carbon Time Series

```julia
# Create a carbon time series entry
carbon_ts = CarbonTimeSeries(
    index = 1,
    timestamp = DateTime(2025, 7, 8, 12, 0, 0),
    grid_carbon_intensity_kgCO2e_per_MWh = 350.0,
    renewable_generation_carbon_intensity_kgCO2e_per_MWh = 25.0,
    storage_carbon_intensity_kgCO2e_per_MWh = 120.0
)
```

### Carbon Scenarios

```julia
# Create a carbon scenario
carbon_scenario = CarbonScenario(
    index = 1,
    name = "Net Zero 2050",
    description = "Scenario aligned with net zero emissions by 2050",
    year = 2050,
    grid_carbon_intensity_kgCO2e_per_MWh = 50.0,
    renewable_penetration_percent = 85.0,
    ev_adoption_percent = 95.0,
    storage_capacity_mwh = 500000.0
)
```

### Equipment Carbon Models

```julia
# Create an equipment carbon record
equipment_carbon = EquipmentCarbon(
    index = 1,
    element_type = "transformer",
    element_id = 101,
    carbon_embodied_kgCO2e = 25000.0,
    carbon_operational_kgCO2e_per_year = 500.0,
    lifetime_years = 30,
    manufacturing_date = Date(2023, 3, 15),
    installation_date = Date(2023, 6, 10),
    recycling_rate_percent = 75.0
)
```

### Carbon Reduction Strategies

The carbon emission models support various strategies for reducing carbon emissions in power systems:

1. **Equipment Replacement**: Calculate the carbon payback period for replacing high-emission equipment.
2. **Scenario Analysis**: Compare system emissions under different future scenarios.
3. **Temporal Optimization**: Schedule operations to minimize carbon emissions based on time-varying carbon intensities.

## Advanced Applications

ComponentModel components can be used for various power system analyses, including:

1. **Power Flow Analysis**: Determining voltage and power distribution under steady-state conditions
2. **Fault Analysis**: Assessing the impact of short-circuit faults on the system
3. **Stability Analysis**: Studying the dynamic behavior of the system following disturbances
4. **Optimal Dispatch**: Determining the optimal output of generators to meet load demands
5. **HVDC Transmission**: Long-distance bulk power transfer with minimal losses
6. **Renewable Integration**: Connection of remote renewable energy sources to the grid
7. **Microgrids**: Formation of DC or hybrid AC/DC microgrids
8. **Energy Storage**: Integration of large-scale energy storage systems
9. **Environmental Impact Assessment**: Evaluating carbon emissions and reduction strategies

For more advanced application examples, refer to the relevant sections in the user manual.
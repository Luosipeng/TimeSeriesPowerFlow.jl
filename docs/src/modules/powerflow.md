# PowerFlow Module

## Overview

The PowerFlow module provides comprehensive power flow analysis capabilities for electrical power systems. It includes implementations of various power flow algorithms including AC Newton-Raphson, DC power flow, and hybrid power flow for integrated AC/DC systems. The module supports both balanced and unbalanced power flow analysis with various features for handling different system components and control modes.

## Flow Chart

![PowerFlow flow chart](../assets/powerflow.png)

## Core Features

- AC power flow using Newton-Raphson method
- DC power flow analysis
- Hybrid power flow for integrated AC/DC systems
- Unbalanced power flow analysis
- GPU-accelerated computations
- Support for ZIP load models
- PV system integration
- Comprehensive results processing and reporting

## Main Functions

### Power Flow Execution

#### `runpf(case, opt)`

Main function to execute AC power flow calculations.

**Arguments:**
- `case`: Power system case data
- `opt`: Options dictionary containing power flow parameters

**Returns:**
- Results structure containing power flow solution

#### `rundcpf(case, opt)`

Executes DC power flow analysis.

**Arguments:**
- `case`: Power system case data
- `opt`: Options dictionary containing DC power flow parameters

**Returns:**
- Results structure containing DC power flow solution

#### `runhpf(jpc, opt)`

Runs hybrid power flow calculation for integrated AC/DC systems.

**Arguments:**
- `jpc`: MATPOWER-style power flow case structure
- `opt`: Options dictionary containing hybrid power flow parameters

**Returns:**
- Results structure containing hybrid power flow solution

#### `runupf(case, jpc_3ph, gs_eg, bs_eg, opt)`

Executes unbalanced power flow analysis on unsymmetrical load nodes.

**Arguments:**
- `case`: Power system case data
- `jpc_3ph`: Three-phase power flow case structure
- `gs_eg`, `bs_eg`: Shunt conductance and susceptance
- `opt`: Options dictionary

**Returns:**
- Results structure containing unbalanced power flow solution

### Solution Algorithms

#### `newtonpf(baseMVA, bus, gen, branch, Ybus, Yf, Yt, V0, ref, pv, pq, mpopt)`

Solves power flow using Newton-Raphson method.

**Arguments:**
- `baseMVA`: Base MVA
- `bus`, `gen`, `branch`: System component data
- `Ybus`, `Yf`, `Yt`: Admittance matrices
- `V0`: Initial voltage vector
- `ref`, `pv`, `pq`: Bus type indices
- `mpopt`: MATPOWER options structure

**Returns:**
- Voltage solution and convergence information

#### `newtondcpf(baseMVA, bus, branch, Bbus, Bf, Pbusinj, Pfinj, Va0, ref, pv, pq, mpopt)`

Solves DC power flow using Newton method.

#### `adaptive_damped_newton(baseMVA, bus, gen, branch, Ybus, Yf, Yt, V0, ref, pv, pq, mpopt)`

Implements adaptive damped Newton method for improved convergence.

#### `currentinjectionpf(baseMVA, bus, gen, branch, Ybus, Yf, Yt, V0, ref, pv, pq, mpopt)`

Solves power flow using current injection method.

### System Construction

#### `makeYbus(baseMVA, bus, branch)`

Builds the bus admittance matrix (Ybus) and branch admittance matrices (Yf, Yt).

**Arguments:**
- `baseMVA`: Base MVA
- `bus`: Bus data matrix
- `branch`: Branch data matrix

**Returns:**
- `Ybus`: Bus admittance matrix
- `Yf`: From-bus branch admittance matrix
- `Yt`: To-bus branch admittance matrix

#### `makeBdc(baseMVA, bus, branch)`

Builds the DC power flow matrices.

**Arguments:**
- `baseMVA`: Base MVA
- `bus`: Bus data matrix
- `branch`: Branch data matrix

**Returns:**
- `Bbus`: Nodal susceptance matrix
- `Bf`: Branch susceptance matrix for line flows
- Other DC model matrices

#### `makeSbus(baseMVA, bus, gen, Vm, load, pvarray; dc=false, Sg=nothing, return_derivative=false)`

Builds the vector of complex bus power injections.

**Arguments:**
- `baseMVA`: Base MVA
- `bus`: Bus data matrix
- `gen`: Generator data matrix
- `Vm`: Voltage magnitude vector
- `load`: Load data matrix
- `pvarray`: PV array data
- Additional optional parameters

**Returns:**
- Vector of complex bus power injections

#### `calculate_line_parameter(net, jpc, sequence, opt)`

Calculates line parameters for the branch matrix based on the specified sequence.

**Arguments:**
- `net`: Network data structure
- `jpc`: MATPOWER-style power flow case structure
- `sequence`: Sequence type (1 for positive, 2 for negative, 0 for zero)
- `opt`: Options dictionary

**Returns:**
- Updated jpc structure with branch parameters

### Bus and Generator Construction

#### `calculate_bus(net, jpc, sequence, slack_bus, opt)`

Builds the bus matrix for power flow calculations.

#### `build_gen(net, jpc)`

Builds the generator matrix for power flow calculations.

### Bus Type Identification

#### `bustypes(bus)`

Identifies bus types (reference, PV, PQ) for AC power flow.

#### `dcbustypes(bus)`

Identifies bus types for DC power flow.

### Jacobian Calculation

#### `dSbus_dV(Ybus, V, vcart)`

Computes partial derivatives of power injection w.r.t. voltage.

**Arguments:**
- `Ybus`: Bus admittance matrix
- `V`: Complex bus voltage vector
- `vcart`: Boolean flag for cartesian coordinates

**Returns:**
- Partial derivatives of power injection

### Solution Processing

#### `pfsoln(baseMVA, bus, gen, branch, Ybus, Yf, Yt, V, ref, pv, pq, opt)`

Updates power flow solution with computed voltage values and calculates generator outputs and branch flows.

#### `dcpfsoln(baseMVA, bus0, gen0, branch0, Bbus, Bf, Pbusinj, Pfinj, Va, success, et, opt)`

Updates DC power flow solution.

#### `merge_results(results)`

Merges power flow results from multiple calculations.

#### `process_result(results, case, opt)`

Processes power flow calculation results and generates summary.

#### `process_pv_acsystem(case, results, opt)`

Processes PV AC systems and integrates them into the power flow solution.

### Utility Functions

#### `julinsolve(A, b, solver)`

Solves linear system Ax = b using specified solver.

#### `total_load(bus, load)`

Calculates total load in the system.

#### `eliminate_element(mpc, element_type, element_id)`

Eliminates a specific element from the power system model.

### Results Reporting

#### `pf_summary`

Generates comprehensive power flow summary reports including:
- System summary
- Bus data
- Branch data
- Generator data
- Voltage violations
- Branch violations
- System losses

#### `compare_voltage_results(results, case, reference_file; tolerance_mag=1e-4, tolerance_ang=1e-3)`

Compares power flow calculation results with reference values from external tools like ETAP.

## GPU Support

The module includes GPU-accelerated versions of key functions for improved performance on compatible hardware:

- `newtonpf_gpu`
- `makeSbus_gpu`
- `makeSdzip_gpu`

## Settings and Options

#### `settings()`

Defines default settings for power flow calculations including:
- Convergence tolerance
- Maximum iterations
- Algorithm selection
- Output verbosity
- Enforcement options
- Initialization methods

## Usage Example

```julia
using HyDistFlow.PowerFlow

# Load case data
case = case3()

# Set options
opt = PowerFlow.settings()
opt[:pf_alg] = :NR  # Newton-Raphson algorithm
opt[:verbose] = 1

# Run power flow
results = PowerFlow.runpf(case, opt)

# Process results
PowerFlow.process_result(results, case, opt)
```

## Advanced Features

### Hybrid Power Flow

The module supports hybrid power flow for integrated AC/DC systems with various converter control modes:
- Vθ control
- PQ control
- PV control
- Vdc-Q control
- Vdc-Vac control
- Droop control

### Unbalanced Power Flow

Unbalanced power flow analysis is supported for unsymmetrical load nodes:
1. Identifies unbalanced nodes in the system
2. Finds interface branches between balanced and unbalanced nodes
3. Solves the unbalanced power flow

### ZIP Load Models

The module supports ZIP (constant impedance, constant current, constant power) load models through the `makeSdzip` function.
"""
    dSbus_dV(Ybus, V, vcart=0)

Compute partial derivatives of power injection with respect to voltage.

# Arguments
- `Ybus`: Bus admittance matrix
- `V`: Complex bus voltage vector
- `vcart`: Flag that determines output format:
  - 0: Derivatives with respect to voltage magnitude and angle (polar form, default)
  - 1: Derivatives with respect to real and imaginary voltage components (Cartesian form)

# Returns
- `dSbus_dV1`: Partial derivatives of power injection with respect to:
  - Voltage angle (∂S/∂θ) if vcart=0
  - Real part of voltage (∂S/∂Vr) if vcart=1
- `dSbus_dV2`: Partial derivatives of power injection with respect to:
  - Voltage magnitude (∂S/∂|V|) if vcart=0
  - Imaginary part of voltage (∂S/∂Vi) if vcart=1

# Description
This function computes the partial derivatives of the complex power injection
at each bus with respect to voltage components. These Jacobian matrices are
essential for power flow analysis, particularly in Newton-Raphson methods.

For polar coordinates (vcart=0), the function returns:
- ∂S/∂θ: Derivatives with respect to voltage angle
- ∂S/∂|V|: Derivatives with respect to voltage magnitude

For Cartesian coordinates (vcart=1), the function returns:
- ∂S/∂Vr: Derivatives with respect to real part of voltage
- ∂S/∂Vi: Derivatives with respect to imaginary part of voltage

The calculations use the relationships:
- S = V ⋅ conj(I) = V ⋅ conj(Y⋅V)
- I = Y⋅V

Note: All derivatives are complex matrices of size n×n, where n is the number of buses.
"""
function dSbus_dV(Ybus, V, vcart=0)
    n = length(V)
    Ibus = Ybus * V

    diagV = Diagonal(V)
    diagIbus = Diagonal(Ibus)

    if vcart == 0
        diagVnorm = Diagonal(V ./ abs.(V))
    end

    if vcart == 1
        dSbus_dV1 = conj(diagIbus) + diagV * conj(Ybus)  # dSbus/dVr
        dSbus_dV2 = 1im * (conj(diagIbus) - diagV * conj(Ybus))  # dSbus/dVi
    else
        dSbus_dV1 = 1im * diagV * conj(diagIbus - Ybus * diagV)  # dSbus/dVa
        dSbus_dV2 = diagV * conj(Ybus * diagVnorm) + conj(diagIbus) * diagVnorm  # dSbus/dVm
    end

    return dSbus_dV1, dSbus_dV2
end

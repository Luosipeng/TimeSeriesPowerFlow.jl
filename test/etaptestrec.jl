#This function is used to test the etap hybrid power flow result
using ModelingToolkit
using NonlinearSolve
const baseMVA = 10.0
#base KV in AC side is 0.4KV
#base V in DC side is 0.824KV
@variables Vm_ac2,Va_ac2
@variables Vm_dc1,Vm_dc2
@parameters P_ac2,Q_ac2
@parameters P_dc1,P_dc2
@parameters Vm_ac1,Va_ac1
@parameters Vm_dc3
@parameters G_ac11,B_ac11,G_ac12,B_ac12
@parameters G_ac21,B_ac21,G_ac22,B_ac22
@parameters G_dc11,G_dc12,G_dc13
@parameters G_dc21,G_dc22,G_dc23
@parameters G_dc31,G_dc32,G_dc33
@parameters eta,Ps,Qs,Ps_dc

eqs1=[
    Vm_ac2*Vm_ac1*(G_ac21*cos(Va_ac2-Va_ac1)+B_ac21*sin(Va_ac2-Va_ac1))+Vm_ac2^2*(G_ac22)-Ps-P_ac2 ~ 0,
    Vm_ac2*Vm_ac1*(G_ac21*sin(Va_ac2-Va_ac1)-B_ac21*cos(Va_ac2-Va_ac1))+Vm_ac2^2*(-B_ac22)-Qs-Q_ac2 ~ 0,
]

eqs2=[
    Vm_dc1*Vm_dc2*(G_dc12)+Vm_dc1^2*(G_dc11)+Ps_dc-P_dc1~ 0,
    Vm_dc2*Vm_dc1*(G_dc21)+Vm_dc2^2*(G_dc22)+Vm_dc2*Vm_dc3*(G_dc23)-P_dc2 ~ 0,
]
eqs=vcat(eqs1,eqs2)
parameters=[P_ac2,Q_ac2,Vm_ac1,Va_ac1,G_ac11,B_ac11,G_ac12,B_ac12,G_ac21,B_ac21,G_ac22,B_ac22,
    Vm_dc3,P_dc1,P_dc2,G_dc11,G_dc12,G_dc13,G_dc21,G_dc22,G_dc23,G_dc31,G_dc32,G_dc33,eta,Ps,Qs,Ps_dc]
@mtkbuild ns=NonlinearSystem(eqs,[Vm_dc1,Vm_dc2,Vm_ac2,Va_ac2],parameters)

u0=[
    Vm_ac2 => 1.0,
    Va_ac2 => 0.0,
    Vm_dc1 => 1.0,
    Vm_dc2 => 1.0,
]
para=[
    P_ac2 => -0.17/baseMVA,
    Q_ac2 => -0.1054/baseMVA,
    Vm_ac1 => 1.0,
    Va_ac1 => 0.0,
    G_ac11 => 1,
    B_ac11 => -1,
    G_ac12 => -1,
    B_ac12 => 1,
    G_ac21 => -1,
    B_ac21 => 1,
    G_ac22 => 1,
    B_ac22 => -1,
    Vm_dc3 => 1.0,
    P_dc1 => -0.2/baseMVA,
    P_dc2 => 0.0/baseMVA,
    G_dc11 => 0.3395/2.0,
    G_dc12 => -0.3395/2.0,
    G_dc13 => 0,
    G_dc21 => -0.3395/2.0,
    G_dc22 => 5.9049/2.0,
    G_dc23 => -5.5654/2.0,
    G_dc31 => 0,
    G_dc32 => -5.5654/2.0,
    G_dc33 => 5.5654/2.0,
    eta => 0.9,
    Ps => -0.013/baseMVA,
    Qs => -0.05/baseMVA,
    Ps_dc => -0.0144/baseMVA,
]
prob=NonlinearProblem(ns,u0,para,jac=true)
sol = solve(prob, NewtonRaphson(); abstol=1e-8, reltol=1e-8, maxiters=100)
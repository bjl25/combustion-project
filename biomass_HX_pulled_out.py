# -*- coding: utf-8 -*-
"""
Created on Mon Sep  4 14:18:30 2023

@author: bjl25
"""
#Importing required pyomo and idaes components
from pyomo.environ import (
    Constraint,
    Var,
    ConcreteModel,
    Expression,
    Objective,
    SolverFactory,
    TransformationFactory,
    value,
    units as pyunits
)
from pyomo.network import Arc, SequentialDecomposition

#Todo add the four other unit operations
from idaes.models.unit_models import (
Mixer,
StoichiometricReactor,
Heater,
HeatExchanger
)

from idaes.models.unit_models.pressure_changer import ThermodynamicAssumption
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core import FlowsheetBlock
# Import idaes logger to set output levels
import idaes.logger as idaeslog
from idaes.models.properties.modular_properties import GenericParameterBlock
from  single_comp_biomass_comb_pp import configuration 
from  biomass_combustion_rp import BMCombReactionParameterBlock

#helmholtz import for water
from idaes.models.properties.general_helmholtz import (
        HelmholtzParameterBlock,
        HelmholtzThermoExpressions,
        AmountBasis,
        PhaseType,
    )
from idaes.models.unit_models.heat_exchanger import delta_temperature_amtd_callback, HX0DInitializer

from idaes.core.util.model_statistics import degrees_of_freedom

m = ConcreteModel()

m.fs = FlowsheetBlock(dynamic=False)

m.fs.biomass_properties = GenericParameterBlock(**configuration)
m.fs.reaction_params = BMCombReactionParameterBlock(
    property_package=m.fs.biomass_properties
)
m.fs.reaction_params.h.fix(0.06)
m.fs.reaction_params.w.fix(0.09)

m.fs.steam_properties = HelmholtzParameterBlock(
        pure_component="h2o", amount_basis=AmountBasis.MOLE,
        phase_presentation=PhaseType.LG,
    )


#flue-water heat exchanger
m.fs.E101 = HeatExchanger(
    delta_temperature_callback=delta_temperature_amtd_callback,
    hot_side_name="shell",
    cold_side_name="tube",
    shell={"property_package": m.fs.biomass_properties},
    tube={"property_package": m.fs.steam_properties}
)

#reactor flow sheet feed via mixer -> reactor -> product via separator


flowTotal = 1
extentR1 = 0.005*flowTotal #absolute extent of reaction

#Inlet
m.fs.E101.shell_inlet.flow_mol.fix(1.0200)
m.fs.E101.shell_inlet.mole_frac_comp[0,"N2"].fix(0.68627)
m.fs.E101.shell_inlet.mole_frac_comp[0,"O2"].fix(0.25980)
m.fs.E101.shell_inlet.mole_frac_comp[0,"CO2"].fix(0.029412)
m.fs.E101.shell_inlet.mole_frac_comp[0,"H2O"].fix(0.024510) 
m.fs.E101.shell_inlet.mole_frac_comp[0,"CO"].fix(1.1869e-12) 
m.fs.E101.shell_inlet.mole_frac_comp[0,"biomass"].fix(1.19e-12) 
m.fs.E101.shell_inlet.temperature.fix(912.65)
m.fs.E101.shell_inlet.pressure.fix(101325)

#specifying heat exchanger
m.fs.E101.area.fix(0.25)
m.fs.E101.overall_heat_transfer_coefficient[0].fix(100)
m.fs.E101.tube_inlet.flow_mol.fix(0.1)
m.fs.E101.tube_inlet.pressure.fix(101325)
m.fs.E101.tube_inlet.enth_mol.fix(m.fs.steam_properties.htpx(p=101325*pyunits.Pa,T=290*pyunits.K))

initializer = HX0DInitializer()
initializer.initialize(m.fs.E101)

m.fs.E101.shell_outlet.temperature.fix(600)
#m.fs.E101.area.unfix()
m.fs.E101.tube_inlet.flow_mol.unfix()
m.fs.E101.tube_outlet.enth_mol.fix(m.fs.steam_properties.htpx(p=101325*pyunits.Pa,T=400*pyunits.K))
m.fs.E101.tube_outlet.enth_mol.unfix()
#m.fs.E101.tube_outlet.enth_mol.fix(m.fs.steam_properties.htpx(p=101325*pyunits.Pa,T=400*pyunits.K))

print(degrees_of_freedom(m))

solver=SolverFactory("ipopt")
status=solver.solve(m,tee=True)

print(degrees_of_freedom(m))
assert degrees_of_freedom(m) == 0

# print(value(m.fs.R101.outlet.temperature[0]))
m.fs.E101.report()

import inspect
import cantera_adaptive_testing.models as models
import cantera_adaptive_testing.surfaces as surfs

mods = {k:v for k,v in inspect.getmembers(models, inspect.isclass)}
ms = ('Hydrogen', 'DME', 'JetA', 'Butane', 'TwoButonane', 'IsoButene', 'NHeptane','IsoOctane', 'ThreeMethylHeptane', 'NHexadecane', 'MethylFiveDeconate', 'MethylDeconateNHeptane', 'TwoMethylnonadecane')

mod_instances = [ mods[m](preconditioned=True, runtype="steady", verbose=True) for m in ms]

# add surfaces
surf = surfs.PlatinumLarge()
for m in mod_instances:
    m.add_surface(surf)

for mi in mod_instances:
    mi.well_stirred_reactor()

for mi in mod_instances:
    mi.plug_flow_reactor()

for mi in mod_instances:
    mi.network_combustor_exhaust()

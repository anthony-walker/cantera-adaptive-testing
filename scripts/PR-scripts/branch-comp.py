import inspect
import cantera_adaptive_testing.models as models
import cantera_adaptive_testing.surfaces as surfs

mods = {k:v for k,v in inspect.getmembers(models, inspect.isclass)}
ms = ('Hydrogen', 'A2', 'A3', 'C1', 'C5', 'Butane', 'TwoButonane', 'IsoButene', 'NHeptane','IsoOctane', 'ThreeMethylHeptane', 'NHexadecane', 'MethylFiveDeconate', 'MethylDeconateNHeptane', 'TwoMethylnonadecane')

mod_instances = [ mods[m](preconditioned=True, endtime=0.001, verbose=True) for m in ms]

# for m in ms:
#     mods[m].print_model_information()
# add surfaces
# surf = surfs.PlatinumLarge()
# for m in mod_instances:
#     m.add_surface(surf)

for mi in mod_instances:
    mi.well_stirred_reactor()

for mi in mod_instances:
    mi.plug_flow_reactor()

for mi in mod_instances:
    mi.network_combustor_exhaust()

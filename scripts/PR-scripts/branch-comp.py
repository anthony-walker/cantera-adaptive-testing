import inspect
import cantera_adaptive_testing.models as models

mods = {k:v for k,v in inspect.getmembers(models, inspect.isclass)}
ms = ('Hydrogen', 'DME', 'JetA', 'Butane', 'TwoButonane', 'IsoButene', 'NHeptane','IsoOctane', 'ThreeMethylHeptane', 'NHexadecane', 'MethylFiveDeconate', 'MethylDeconateNHeptane', 'TwoMethylnonadecane')
mod_instances = [ mods[m](preconditioned=True) for m in ms]
for mi in mod_instances:
    mi.well_stirred_reactor()
    mi.plug_flow_reactor()

from cantera_adaptive_testing.models import *

m = PlatinumMediumHydrogen(runtype="performance", preconditioned=True, threshold=1e-12, remove_falloff=True, remove_thirdbody=True)
m.network_combustor_exhaust()
# m.plug_flow_reactor()

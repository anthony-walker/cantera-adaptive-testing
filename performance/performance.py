import cantera_adaptive_testing.mechanisms as mechanisms
import importlib
mechs = inspect.getmembers(mechanisms, inspect.isclass)
mechs = {element[0]: element[1] for element in mechs}

# ja = mech.GRI(verbose=True, plot=True, tight=True)
# ja.advance()

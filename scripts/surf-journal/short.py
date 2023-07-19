import os
import cantera_adaptive_testing.models as models
from cantera_adaptive_testing.surfaces import *
import cantera as ct
import ruamel.yaml
import inspect

name = "short_symp"
direc = f"{name}_data"
yml = f"{name}.yml"

def run_cases():
    mods = [models.Hydrogen, models.GRI, models.A2, models.A3, models.C1, models.C5, models.Butane, models.TwoButonane, models.IsoButene, models.NHeptane, models.IsoOctane, models.Toluene, models.NHexadecane, models.MethylFiveDeconate, models.MethylDeconateNHeptane, models.TwoMethylnonadecane]
    # create steady state times
    T = 1000
    surf=True
    for m in mods:
        apsf = False #True if "Hydrogen" not in m.__name__ else False
        print(m, apsf)
        cm = m(log=False, runtype="steady", max_time_step=1e-4, preconditioned=True, append_surface_fuel=apsf, T=T, P=ct.one_atm)
        if surf:
            cs = PlatinumLarge()
            cm.add_surface(cs)
        cm.plug_flow_reactor()
        # performance runs
        cm = m(log=True, out_dir=direc, runtype="performance", preconditioned=False, append_surface_fuel=apsf, T=T, P=ct.one_atm)
        if surf:
            cs = PlatinumLarge()
            cm.add_surface(cs)
        cm.plug_flow_reactor()
        # performance runs
        cm = m(log=True, out_dir=direc, runtype="performance", preconditioned=True, append_surface_fuel=apsf, T=T, P=ct.one_atm)
        if surf:
            cs = PlatinumLarge()
            cm.add_surface(cs)
        cm.plug_flow_reactor()


def combine_surf_yamls(direc="performance_data", yml_name="performance.yaml"):
    yaml = ruamel.yaml.YAML()
    files = os.listdir(direc)
    files = list(filter(lambda x: ".yaml" in x, files))
    mods, ___ = zip(*inspect.getmembers(models, inspect.isclass))
    # sort all files into cases
    fp = lambda x: os.path.join(direc, x)
    cases = {}
    for m in mods:
        curr_files = list(filter(lambda x: m in x, files))
        for cf in curr_files:
            with open(fp(cf), "r") as f:
                data = yaml.load(f)
            case_key = list(data.keys())[0]
            if case_key in cases:
                cases[case_key].append(cf)
            else:
                cases[case_key] = [cf]
    # combine cases
    case_data = {}
    for case in cases:
        case_files = cases[case]
        case_data[case] = {}
        for cf in case_files:
            with open(fp(cf), "r") as f:
                data = yaml.load(f)
            if data:
                for k1, v1 in data.items():
                    for k2, v2 in v1.items():
                        # do nothing if exception is a key
                        if "exception" in v1[k2]:
                            pass
                        elif k2 not in case_data[case]:
                            case_data[case][k2] = v1[k2]
                        else:
                            case_data[case][k2]["runtime"] += v2["runtime"]
                            case_data[case][k2]["nruns"] += 1
    # average the run times
    for c in case_data:
        for p in case_data[c]:
            case_data[c][p]["runtime"] /= case_data[c][p]["nruns"]
    with open(yml_name, 'w') as f:
        data = yaml.dump(case_data, f)

if __name__ == "__main__":
    # run_cases()
    combine_surf_yamls(direc=direc, yml_name=yml)

import re
import random
from sympy.parsing.sympy_parser import parse_expr


def extract_from_fac(facfile, replace_generics=True, replace_complex=True):
    prefix = facfile.split(".")[0]
    # break down file into sections
    with open(facfile, "r") as f:
        content = f.read()
        delimiters = re.findall(r"[*][;]\n[*].*[;][\n][*][;]\n", content)
        sections = []
        for d in delimiters[::-1]:
            content, app = content.split(d)
            sections.append(app)
        sections = sections[::-1]
    # write species file
    species = re.findall(r"(?:[A-Z]+[a-z]*[0-9]*)+", sections.pop(0))
    species.remove("VARIABLE")
    with open(f"{prefix}-species.txt", "w") as f:
        for sp in species:
            f.write(f"{sp.strip()}\n")
    # substitute in generic rate coefficients where possible
    gcs = {}
    generic_coeffs = re.sub(r"\s*", "", sections.pop(0))
    for gc in generic_coeffs.split(";")[:-1]:
        gc_key, gc_rate = gc.split("=")
        gc_rate = re.sub(r"EXP", r"exp", gc_rate)
        gc_rate = re.sub(r"(\d+)[DE]([+-]?)(\d+)", r"\1e\2\3", gc_rate)
        gc_rate = re.sub(r"[@]([-]?\d+([.]\d+)?([e][+-]?\d+)?)", r"**(\1)", gc_rate)
        gcs[gc_key] = gc_rate
    # complex coefficients
    # complex_coeffs = sections.pop(0)
    # substitute in generic rate coefficients where possible
    ccs = {}
    complex_coeffs = re.sub(r"\s*", "", sections.pop(0))
    complex_coeffs = re.sub(r"[*]+[;]", "", complex_coeffs)
    for cc in complex_coeffs.split(";")[:-1]:
        cc_key, cc_rate = cc.split("=")
        cc_rate = re.sub(r"EXP", r"exp", cc_rate)
        cc_rate = re.sub(r"(\d+)[DE]([+-]?)(\d+)", r"\1e\2\3", cc_rate)
        cc_rate = re.sub(r"[@]([-]?\d+([.]\d+)?([e][+-]?\d+)?)", r"**(\1)", cc_rate)
        ccs[cc_key] = cc_rate
    # get the simplified complex expressions
    ccs_keys = ccs.keys()
    change = True
    temp_data = {x: y for x, y in ccs.items()}
    while change:
        change = False
        crates = sorted(temp_data.items(), key=lambda x: len(x[0]), reverse=True)
        for rname, rexpr in crates:
                found = re.findall(r"(?:[A-Z]+[0-9]*)+", rexpr)
                found.sort(key=lambda x: len(x), reverse=True)
                for fexp in found:
                    if fexp in ccs_keys:
                        rexpr = rexpr.replace(fexp, ccs[fexp])
                        change = True
                temp_data[rname] = rexpr
    for rname, rate in temp_data.items():
        rate = re.sub(r"EXP", r"exp", rate)
        rate = re.sub(r"(\d+)[DE]([+-]?)(\d+)", r"\1e\2\3", rate)
        rate = re.sub(r"[@]", r"**", rate.strip())
        ccs[rname] = str(parse_expr(rate))
    # write ro2 species file
    species = re.findall(r"(?:[A-Z]+[a-z]*[0-9]*)+", sections.pop(0).split("=")[1])
    with open(f"{prefix}-ro2-sum.txt", "w") as f:
        for sp in species:
            f.write(f"{sp.strip()}\n")
    # separate reactions and rates
    equations = re.findall(r"[%].*[;]", sections.pop(0))
    rrs = []
    for eq in equations:
        ceq = re.sub(r"[%](.*)[;]", r"\1", eq)
        rate, reaction = ceq.split(":")
        rate = re.sub(r"EXP", r"exp", rate)
        rate = re.sub(r"(\d+)[DE]([+-]?)(\d+)", r"\1e\2\3", rate)
        rate = re.sub(r"[@]([-]?\d+([.]\d+)?([e][+-]?\d+)?)", r"**(\1)", rate.strip())
        # check for the generic rate
        rrs.append((reaction.strip(), rate))
    rrs.sort()
    reactions, rates = zip(*rrs)
    # replace rates with generic rates maybe
    if replace_generics:
        rate_str = "\n".join(rates)
        sort_gc = list(gcs.items())
        sort_gc.sort(key=lambda x: len(x[0]), reverse=True)
        for gc, gc_exp in sort_gc:
            rate_str = rate_str.replace(gc, gc_exp)
        rates = rate_str.split("\n")
    if replace_complex:
        arrhen_regex = r"\d+([.]\d*)?([e][+-]\d+)?(([*]\d+([.]\d*)?([e][+-]\d+)?)+)?([*][(]TEMP[/]\d+[)][*][*][-]?\d+([.]\d*)?)?([*]exp[(][-]?\d+[/]TEMP[)])?(([*]\d+([.]\d*)?([e][+-]\d+)?)+)?"
        rate_str = "\n".join(rates)
        sort_cc = list(ccs.items())
        sort_cc.sort(key=lambda x: len(x[0]), reverse=True)
        for cc, cc_exp in sort_cc:
            if re.fullmatch(arrhen_regex,cc_exp):
                rate_str = rate_str.replace(cc, cc_exp)
        rates = rate_str.split("\n")
    # write reactions file
    with open(f"{prefix}-reactions.txt", "w") as f:
        f.write("\n".join(reactions))
    # write rates file
    with open(f"{prefix}-rates.txt", "w") as f:
        f.write("\n".join(rates))


if __name__ == "__main__":
    extract_from_fac("full-mcm.fac", replace_generics=True, replace_complex=False)

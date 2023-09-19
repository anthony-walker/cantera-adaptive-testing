import re
import random

def extract_from_fac(facfile, replace_generics=True):
    prefix = "mcm"
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
    complex_coeffs = sections.pop(0)
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
    # write reactions file
    with open(f"{prefix}-reactions.txt", "w") as f:
        f.write("\n".join(reactions))
    # write rates file
    with open(f"{prefix}-rates.txt", "w") as f:
        f.write("\n".join(rates))


if __name__ == "__main__":
    extract_from_fac("full-mcm.fac", replace_generics=True)

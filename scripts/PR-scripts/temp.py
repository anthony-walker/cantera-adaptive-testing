import re

with open("bc.txt") as f:
    lines = f.read()
lines.strip()
lines = list(filter(None, lines.split("\n")))
lines.sort()
lines = [lines[i:i+2] for i in range(0, len(lines), 2)]
lns = []
for l1, l2 in lines:
    print(l1, l2)
    name, prob, bran, t1, _, nspec = l1.split(" ")
    name, prob, bran, t2, _, nspec = l2.split(" ")
    t1 =  float(t1)
    t2 = float(t2)
    bw = t1 - t2
    name = name[:10] if len(name) >= 10 else name + (10 - len(name)) * " "
    prob = prob[:10] if len(prob) >= 10 else prob + (10 - len(prob)) * " "
    nspec = nspec[:4] if len(nspec) >= 4 else nspec + (4 - len(nspec)) * " "
    ostr = f"{name} {prob} {nspec} {t1:3.3f}    |    {t2:3.3f}    {bw:3.3f}\n"
    lns.append(ostr)

with open("comp.txt", "w") as f:
    f.write("Model     problem     nspecies   main    |    sdev    main - sdev\n")
    f.write("IdealGasMoleReactor\n")
    for l in filter(lambda x: "plug" not in x, lns):
        f.write(l)
    f.write("\n")
    f.write("IdealGasConstPressureMoleReactor\n")
    for l in filter(lambda x: "well" not in x, lns):
        f.write(l)

import os, subprocess

outputs = os.listdir("./slurm-output")
scancel_str = "scancel {:s}\n"
cmds = []
for f in outputs:
	flist = f.split("-") 
	type_str = flist[0]
	slurm_id = flist[-1].split(".")[0]
	if type_str == "anprec":
		cmds.append(scancel_str.format(slurm_id))

with open("cancel_script.sh","w") as f:
	for cmd in cmds:
		f.write(cmd)


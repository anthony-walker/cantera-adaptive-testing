import os, subprocess, sys

def cancel_type():
    if len(sys.argv) > 1:
        ftype = sys.argv[1]
        outputs = os.listdir("./slurm-output")
        scancel_str = "scancel {:s}\n"
        cmds = []
        for f in outputs:
            flist = f.split("-")
            type_str = flist[0]
            slurm_id = flist[-1].split(".")[0]
            if type_str == ftype:
                cmds.append(scancel_str.format(slurm_id))

        with open("cancel_script.sh","w") as f:
            for cmd in cmds:
                f.write(cmd)
    else:
        print("No type string passed as command line argument")

if __name__ == "__main__":
    cancel_type()

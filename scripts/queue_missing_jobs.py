import os, subprocess, sys

def queue_missing_jobs():
    if len(sys.argv) > 2:
        dirname = sys.argv[1]
        opts_type = sys.argv[2]
        command = "adaptive-utilities {:s} count_uniq_yamls".format(dirname)

        process = subprocess.Popen(command.split(" "), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = process.communicate()
        out = str(out)
        process.terminate()
        # Fix first element and remove last
        outlist = out.split("\\n")[:-1]
        outlist[0] = outlist[0][2:]
        # opts_file model script_prefix nruns memory(optional)
        launch_cmd = "./one-launch.sh {:s} {:s} {:s} {:d}"
        # define thresh_id so it won't fail for mass/mole cases
        thresh_id=0
        for f in outlist:
            name, nruns = f.split(" ")
            nruns = 10 - int(nruns)//10
            if nruns > 0:
                nlist = name.split("-")
                # preconditioned run
                if len(nlist) > 2:
                    model = nlist[0]
                    script_prefix = "-".join(nlist[1:3])
                    curr_cmd = launch_cmd.format(opts_type, model, script_prefix, nruns)
                    if nlist[-1][:-1] == "0.0e+00":
                        thresh_id = 0
                    else:
                        thresh_id = int(nlist[-1][:-1])
                # mass or mole run
                else:
                    model = nlist[0]
                    script_prefix = nlist[1][:-1]
                    curr_cmd = launch_cmd.format(opts_type, model, script_prefix, nruns)
                # run the command

                process = subprocess.Popen(curr_cmd.split(), stdout=subprocess.PIPE, env={"TEND":str(thresh_id), "TSTART":str(thresh_id)})
                out, err = process.communicate()
                for ot in str(out).split("\\n"):
                    print(ot)
                process.terminate()
    else:
        print("No type string passed as command line argument")

if __name__ == "__main__":
    queue_missing_jobs()

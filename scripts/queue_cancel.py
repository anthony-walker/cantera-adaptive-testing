import subprocess
get_curr_job_file = "rm curr_jobs; squeue -u walkanth > curr_jobs"
subprocess.Popen(get_curr_job_file.split(), stdout=subprocess.PIPE)
with open("curr_jobs", "r") as f:
	line = f.readline()
	while line:
		if "anprec" in line:
			jid = line.split()[0]
			process = subprocess.Popen("scancel {:s}".format(jid).split(), stdout=subprocess.PIPE)
			print(process.communicate())
		line = f.readline()
		

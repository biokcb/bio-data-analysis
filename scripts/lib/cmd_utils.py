###############################################
#
# Command line utility functions
# 
###############################################

import subprocess

def run_cmd(cmd):  
	""" 
	Runs a command line tool/command entered. 
	Returns stdout, stderr 
	
	cmd = command to run
	"""
	print('Now running:')
	print(cmd)
	process = subprocess.Popen(cmd,shell=True, stdout=subprocess.PIPE)
	output, error = process.communicate()
	return output, error
	
def main():
	# test with ls command & simply print stdout/stderr
	cmd = 'ls -l'
	out, err = run_cmd(cmd)
	print(out)
	print(err)
	
if __name__ == "__main__":
	main()


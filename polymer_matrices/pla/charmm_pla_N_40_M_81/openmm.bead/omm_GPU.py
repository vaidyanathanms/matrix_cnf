#!/home/charmm-gui/local/bin/python3
import subprocess as sp

def get_njob():
    args = ["nvidia-smi","--query-compute-apps=pid,process_name,used_memory", "--format=csv"]
    result = sp.run(args, check=True, stderr=sp.STDOUT, stdout=sp.PIPE)
    if result.returncode == 0:
        stdout = result.stdout
        output = stdout.decode('utf-8')
        lines = [x for x in output.split('\n') if x.strip()]
        njob = len(lines)-1
        print('njob: ', njob)
    else:
        print("Warning: nvidia-smi does not respond!")
    return njob

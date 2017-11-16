from dask.distributed import Client
client = Client(scheduler_file='scheduler.json')

import socket
host = client.run_on_scheduler(socket.gethostname)

def start_jlab(dask_scheduler):
    import subprocess
    proc = subprocess.Popen(['jupyter', 'lab', '--ip', host, '--no-browser'])
    #proc = subprocess.Popen(['jupyter', 'lab', '--ip', host, '--no-browser','--port','9999'])
    #proc = subprocess.Popen(['nohup','jupyter', 'lab', '--ip', host, '--no-browser'])
    dask_scheduler.jlab_proc = proc

client.run_on_scheduler(start_jlab)

print("ssh -N -L 8787:%s:8787 -L 8888:%s:8888 datarmor1-10g" % (host, host))
#print("ssh -N -L 9999:%s.ib0.ice.ifremer.fr:9999 -L 8888:%s.ib0.ice.ifremer.fr:8888 aponte@datarmor1-10g" % (host, host))


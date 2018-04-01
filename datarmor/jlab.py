#!/usr/bin/env python
import os, sys
import socket
import subprocess


def parse_qstat(job):
    #
    username = os.environ['USER']
    print(username)
    #
    bashCommand = 'qstat -f '+job
    output = subprocess.check_output(bashCommand, shell=True)
    #
    for line in output.splitlines():
        lined = line.decode('UTF-8')
        if 'exec_host' in lined:
            print(lined)
            host = lined.split('=')[1].split('/')[0].strip()
            print(host)
            return host

if __name__ == '__main__':

    jlab_port = '8877'

    user = os.environ['USER']
    #user = 'aponte'

    hostname = 'datarmor1-10g'

    print(sys.argv)
    #host = parse_qstat(sys.argv[1])
    host = socket.gethostname()
    print(f'host ={host}')

    com = ['jupyter', 'lab', '--ip', host, '--no-browser','--port', jlab_port]
    print(' '.join(com))
    proc = subprocess.Popen(com)
    print(f'ssh -N -L {jlab_port}:{host}:{jlab_port}  {user}@{hostname}')
    #print("ssh -N -L 8888:%s:8888 -l aponte aponte@datarmor1-10g" % (host))


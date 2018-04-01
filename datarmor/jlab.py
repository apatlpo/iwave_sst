#!/usr/bin/env python
import os, sys
import socket
import subprocess

if __name__ == '__main__':

    jlab_port = '8877'

    user = os.environ['USER']

    hostname = 'datarmor1-10g'

    host = socket.gethostname()

    com = ['jupyter', 'lab', '--ip', host, '--no-browser','--port', jlab_port]
    print(' '.join(com))
    proc = subprocess.Popen(com)
    print(f'ssh -N -L {jlab_port}:{host}:{jlab_port}  {user}@{hostname}')


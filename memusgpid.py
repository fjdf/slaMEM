#!/usr/bin/env python
import sys
import os
import time
import subprocess

if ( len(sys.argv) != 2 ) :
    print ( "Usage: %s <pid>" % str(sys.argv[0]) )
    sys.exit(2)

pid = int(sys.argv[1])
pidtime = time.time()
pidmem = -1

while True:
    psproc = subprocess.Popen(['ps', '-o', 'pid=', '-o', 'rssize=', '-o', 'comm=', '--pid', str(pid)], stdout=subprocess.PIPE, stderr=None, shell=False)
    (stdout, stderr) = psproc.communicate()
    line = stdout.strip()
    if ( line == "" ) : break
    line = line.split()
    #print ( "\n{[memusgpid](%d,%d)%s}" % (os.getpid(),psproc.pid,str(line)) )
    memvalue = int(line[1])
    pidmem = max( pidmem , memvalue)
    time.sleep(0.1)

if ( pidmem == -1 ) :
    sys.stderr.write("[memusgpid] Error: pid " + str(pid) + " not found\n")
    sys.exit(1)

pidtime = int( time.time() - pidtime )
pidmem = int( round( float(pidmem)/float(1000) ) )
sys.stderr.write("[memusgpid] Real Time = " + str(pidtime) + " s ; RSS Memory = " + str(pidmem) + " MB\n")
sys.stderr.write("[memusgpid] " + str(pidmem) + " " + str(pidtime) + "\n")
sys.exit(0)

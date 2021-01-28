#!/bin/python3
import fileinput, os
import time

idx = 0
for line in fileinput.input():
    os.system('make clean')

    # recorre los distintos flags que tira $ cat flaglist.txt
    os.system('make {} > out{:02d}.txt'.format(line.rstrip(), idx))

    # header de out{:02d}.txt 
    os.system('echo "" >> out{:02d}.txt'.format(idx))
    os.system('echo "############# runnig with perf stat" >> out{:02d}.txt'.format(idx))
    os.system('echo "" >> out{:02d}.txt'.format(idx))

    # perf stat
    os.system('perf stat -r 10 ./tiny_md >> out{:02d}.txt'.format(idx))
    idx+=1

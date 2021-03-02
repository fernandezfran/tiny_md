#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""  graficos: comparación de optimizaciones secuenciales, flags, SoA vs AoS, etc
  """
import os
import numpy as np
import matplotlib.pyplot as plt

os.system("echo '# tiempo, incerteza' > time-aos.dat")
os.system("cat 01-aos/slurm-45640.out | grep seconds | awk '{print $1,$3}' >> time-aos.dat")
os.system("echo '# tiempo, incerteza' > time-soa.dat")
os.system("cat 02-soa/slurm-45641.out | grep seconds | awk '{print $1,$3}' >> time-soa.dat")
os.system("echo '# tiempo, incerteza' > time-mix.dat")
os.system("cat 03-naiv-loops/slurm-45642.out | grep seconds | awk '{print $1,$3}' >> time-mix.dat")

aos_time, aos_time_des = np.loadtxt("time-aos.dat", unpack=True)
soa_time, soa_time_des = np.loadtxt("time-soa.dat", unpack=True)
mix_time, mix_time_des = np.loadtxt("time-mix.dat", unpack=True)

x_aos = np.arange(0,len(aos_time))
x_soa = np.arange(0,len(soa_time))

plt.scatter(x_aos, aos_time, marker="o", s=30, color="tab:red", label="AoS")
plt.errorbar(x_aos, aos_time, aos_time_des, ls='none', color="tab:red")
plt.scatter(x_soa, soa_time, marker="s", s=30, color="tab:green", label="SoA")
plt.errorbar(x_soa, soa_time, soa_time_des, ls='none', color="tab:green")
plt.scatter(x_soa, mix_time, marker="^", s=30, color="tab:blue", label="naive loops")
plt.errorbar(x_soa, mix_time, mix_time_des, ls='none', color="tab:blue")
plt.xticks(x_aos)
plt.xlim(np.min(x_aos)-0.5, np.max(x_aos)+0.5)
plt.xlabel("flags")
plt.ylabel("tiempo [s]")
plt.grid(axis="y")
plt.legend()
plt.savefig("time.png", dpi=600)
plt.show()

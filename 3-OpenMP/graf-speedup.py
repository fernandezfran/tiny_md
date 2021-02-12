#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""  graficos: comparación de eficiencia para distintos tamaños de problema y 
  "            cantidad de hilos en OpenMP
  """
import os
import numpy as np
import matplotlib.pyplot as plt

os.system("cat perf0256bind.txt | grep seconds | awk '{print $1,$3}' > 0256.dat")
os.system("cat perf0500bind.txt | grep seconds | awk '{print $1,$3}' > 0500.dat")
os.system("cat perf0864bind.txt | grep seconds | awk '{print $1,$3}' > 0864.dat")
os.system("cat perf1372bind.txt | grep seconds | awk '{print $1,$3}' > 1372.dat")
os.system("cat perf2048bind.txt | grep seconds | awk '{print $1,$3}' > 2048.dat")

tim_0256 = np.loadtxt("0256.dat", unpack=True)
tim_0500 = np.loadtxt("0500.dat", unpack=True)
tim_0864 = np.loadtxt("0864.dat", unpack=True)
tim_1372 = np.loadtxt("1372.dat", unpack=True)
tim_2048 = np.loadtxt("2048.dat", unpack=True)

su_0256 = (tim_0256[0][0])/tim_0256[0]
su_0500 = (tim_0500[0][0])/tim_0500[0]
su_0864 = (tim_0864[0][0])/tim_0864[0]
su_1372 = (tim_1372[0][0])/tim_1372[0]
su_2048 = (tim_2048[0][0])/tim_2048[0]

x = np.arange(1,29)

plt.xlabel("Hilos")
plt.ylabel("SpeedUp")
plt.grid(axis="y")
plt.scatter(x, su_0256, marker="o", s=30, label="N = 256")
plt.scatter(x, su_0500, marker="o", s=30, label="N = 500")
plt.scatter(x, su_0864, marker="o", s=30, label="N = 864")
plt.scatter(x, su_1372, marker="o", s=30, label="N = 1372")
plt.scatter(x, su_2048, marker="o", s=30, label="N = 2048")
plt.legend()
plt.savefig("graf-speedup.png", dpi=600)
plt.show()

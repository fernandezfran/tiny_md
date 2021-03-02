#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""  graficos: comparación de optimizaciones secuenciales, flags, SoA vs AoS, etc
  """
import os
import numpy as np
import matplotlib.pyplot as plt

os.system("cat 01-aos/out* | grep ns/day | awk '{print $4}' > metrica-aos.dat")
metrica = np.loadtxt("metrica-aos.dat", unpack=True)
metrica = np.split(metrica, 20)
aos_nsday, aos_nsday_des = [], []
for i in range(0, len(metrica)):
    aos_nsday.append(np.mean(metrica[i]))
    aos_nsday_des.append(np.std(metrica[i]))
x_aos = np.arange(0, len(metrica))
idx = np.argmax(aos_nsday)
print("AoS max: %f +- %f (flag_nro: %d)" % (aos_nsday[idx], aos_nsday_des[idx], idx+1))

os.system("cat 02-soa/out* | grep ns/day | awk '{print $4}' > metrica-soa.dat")
metrica = np.loadtxt("metrica-soa.dat", unpack=True)
metrica = np.split(metrica, 20)
soa_nsday, soa_nsday_des = [], []
for i in range(0, len(metrica)):
    soa_nsday.append(np.mean(metrica[i]))
    soa_nsday_des.append(np.std(metrica[i]))
x_soa = np.arange(0, len(metrica))
idx = np.argmax(soa_nsday)
print("SoA max: %f +- %f (flag_nro: %d)" % (soa_nsday[idx], soa_nsday_des[idx], idx+1))

os.system("cat 03-naiv-loops/out* | grep ns/day | awk '{print $4}' > metrica-mix.dat")
metrica = np.loadtxt("metrica-mix.dat", unpack=True)
metrica = np.split(metrica, 20)
mix_nsday, mix_nsday_des = [], []
for i in range(0, len(metrica)):
    mix_nsday.append(np.mean(metrica[i]))
    mix_nsday_des.append(np.std(metrica[i]))
x_soa = np.arange(0, len(metrica))
idx = np.argmax(mix_nsday)
print("Naive loops max: %f +- %f (flag_nro: %d)" % (mix_nsday[idx], mix_nsday_des[idx], idx+1))


plt.scatter(x_aos, aos_nsday, marker="o", s=30, color="tab:red", label="AoS")
plt.errorbar(x_aos, aos_nsday, aos_nsday_des, ls='none', color="tab:red")
plt.scatter(x_soa, soa_nsday, marker="s", s=30, color="tab:green", label="SoA")
plt.errorbar(x_soa, soa_nsday, soa_nsday_des, ls='none', color="tab:green")
plt.scatter(x_soa, mix_nsday, marker="^", s=30, color="tab:blue", label="naive loops")
plt.errorbar(x_soa, mix_nsday, mix_nsday_des, ls='none', color="tab:blue")
plt.xticks(x_aos)
plt.xlim(np.min(x_aos)-0.5, np.max(x_aos)+0.5)
plt.xlabel("flags")
plt.ylabel("ns/day")
plt.grid(axis="y")
plt.legend()
plt.savefig("nsday.png", dpi=600)
plt.show()

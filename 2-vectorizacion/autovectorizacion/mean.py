#!/usr/bin/env python3
# -*- coding:utf-8 -*-
"""   cálculo del promedio de ns/day para las tres versiones del código con
  "   flags de autovectorizacion
  """
import os
import numpy as np

os.system("cat 01-aos/statr30.dat | grep ns/day | awk '{print $4}' > ns.dat")
ns = np.loadtxt("ns.dat", unpack=True)
ns_aos = np.mean(ns)
std_aos = np.std(ns)
os.system("rm ns.dat")

os.system("cat 02-soa/statr30.dat | grep ns/day | awk '{print $4}' > ns.dat")
ns = np.loadtxt("ns.dat", unpack=True)
ns_soa = np.mean(ns)
std_soa = np.std(ns)
os.system("rm ns.dat")

os.system("cat 03-naiv-loops/statr30.dat | grep ns/day | awk '{print $4}' > ns.dat")
ns = np.loadtxt("ns.dat", unpack=True)
ns_naiv = np.mean(ns)
std_naiv = np.std(ns)
os.system("rm ns.dat")

print("AoS:         (%f +- %f) [ns/day]" % (ns_aos, std_aos))
print("SoA:         (%f +- %f) [ns/day]" % (ns_soa, std_soa))
print("Naive loops: (%f +- %f) [ns/day]" % (ns_naiv, std_naiv))

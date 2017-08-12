#!/usr/bin/env python
import os,subprocess
possible_versions=["python2.7","python2","python2.6","python"]
pyv=""
i=0
#print("Trying to find a working Python version")
while pyv=="":
    py=possible_versions[i]
    "Testing command %s"%py
    try:
        with open(os.devnull, 'w') as devnull:
            testnumpy=subprocess.check_output("%s -c \'import numpy as np; print \"OK\"\'"%py,shell=True, stderr=devnull)
    except:
        testnumpy=""
    if testnumpy==b"OK\n":
        pyv=py
    i+=1
if pyv=="":
    print("PYTHON2_WITH_NUMPY_NOT_FOUND")
else:
    print(pyv)

import subprocess
import time
import os
import sys

def generate_el_commands(kernel):
    commands = []
    if (kernel == "TC"):
      sym_test = [True]
    else:
      sym_test = [True, False]
    for symmetrize in sym_test:
        for npes in range(2,6):
            for v in [9,14]:
                for e in [8,9,11,22,23]:
                    el = "v" + str(v) + "_e" + str(e) + ".el"
                    oshmem = "oshrun -np " + str(npes) + " ./" + kernel + " -f EdgeLists/" + el
                    openmp = "../gapbs/gapbs/src/" + kernel + " -f EdgeLists/" + el
                    if (symmetrize): 
                        oshmem = oshmem + " -s -n 3 -v"
                        openmp = openmp + " -s -n 3 -v"
                    else:
                        oshmem = oshmem + " -n 3 -v"
                        openmp = openmp + " -n 3 -v"                       
                    commands.append(oshmem)
                    commands.append(openmp)
    return commands

def generate_synthetic_commands(kernel):
  commands = []
  for npes in [2,3,4,6]:
    for gtype in [" -u ", " -g "]:
      for size in [10, 15]:
        oshmem = "oshrun -np " + str(npes) + " ./" + kernel + gtype + str(size) + " -n 3 -v"
        openmp = "../gapbs/gapbs/src/" + kernel + gtype + str(size) + " -n 3 -v"
        commands.append(oshmem)
        commands.append(openmp)
  return commands

def generate_real_commands(kernel):
  commands = []
  graphs = ["twitter.el", "USA-road-d.USA.gr"]
  for g in graphs:
    for npes in range(4, 5, 6):
      commands.append("oshrun -np " + str(npes) + " ./" + kernel + " -f ../gapbs/gapbs/benchmark/graphs/raw/" + g + " -n 3 -v")
      commands.append("../gapbs/gapbs/src/" + kernel + " -f ../gapbs/gapbs/benchmark/graphs/raw/" + g + " -n 3 -v")
  return commands

def prettify(commands):
    i = 0
    temp = []
    for c in commands:
      if i % 2 == 0: temp.append(c) 
      i+=1
    with open("harness_results.txt", "r") as infile: 
        relevant = []
        for l in infile.readlines():
          if (l.find("Verification:") > -1 and l.find("PE") == -1):
            relevant.append(l)
        infile.close()
        with open("results.txt", "w") as outfile:
          i = 0
          for c in commands:
            outfile.write(c+"\n")
            outfile.writelines(relevant[i:i+3])
            i += 3
commands = generate_synthetic_commands(kernel)
commands = commands + generate_real_commands(kernel)
commands = commands + generate_el_commands(kernel)
prettify(commands)

                        

import subprocess
import time
import os

def generate_el_commands(kernel):
    commands = []
    if (kernel == "TC"):
      sym_test = [True]
    else:
      sym_test = [True, False]
    for symmetrize in sym_test:
        for npes in range(2,6):
            for v in [9,14]:
                for e in [8,9,10,11,22,23]:
                    el = "v" + str(v) + "_e" + str(e) + ".el"
                    oshmem = "oshrun -np " + str(npes) + " ./" + kernel + " -f EdgeLists/" + el
                    openmp = "../gapbs/gapbs/src/" + kernel + " -f EdgeLists/" + el
                    if (symmetrize): 
                        oshmem = oshmem + " -s -n 1 -v"
                        openmp = openmp + " -s -n 1 -v"
                    else:
                        oshmem = oshmem + " -n 1 -v"
                        openmp = openmp + " -n 1 -v"                       
                    commands.append(oshmem)
                    commands.append(openmp)
    return commands

def generate_synthetic_commands(kernel):
  commands = []
  for npes in range(2, 6):
    for gtype in [" -u ", " -g "]:
      for size in [10, 15, 25]:
        oshmem = "oshrun -np " + str(npes) + " ./" + kernel + gtype + str(size) + " -n 1 -v"
        openmp = "../gapbs/gapbs/src/" + kernel + gtype + str(size) + " -n 1 -v"
        commands.append(oshmem)
        commands.append(openmp)
  return commands

def generate_real_commands(kernel):
  commands = []
  graphs = ["twitter.el"]
  for g in graphs:
    for npes in range(4, 5):
      commands.append("oshrun -np " + str(npes) + " ./" + kernel + " -f ../gapbs/gapbs/benchmark/graphs/raw/" + g + " -n 1 -v")
      commands.append("../gapbs/gapbs/src/" + kernel + " -f ../gapbs/gapbs/benchmark/graphs/raw/" + g + " -n 1 -v")
  return commands

def execute_tests(kernel):
    #commands = generate_synthetic_commands(kernel)
    #commands = generate_real_commands(kernel)
    commands = generate_el_commands(kernel)
    with open ("harness_results.txt", "a") as out: 
        for c in commands:
            out.write(c+"\n")
            output = subprocess.call(c, stdout=out, shell=True)
            #time.sleep(10)              # I'd use wait, but does subprocess.call spawn the call as a child of this parent?
            print(c + " | Output code: " + str(output))
        out.close()

def prettify():
    with open("harness_results.txt", "r") as infile: 
        lines = infile.readlines()
        relevant = []
        for l in lines:
          #if (l.find("-f") > -1):
          #  relevant.append(l)
          if (l.find("Verification:") > -1 and l.find("PE") == -1):
            relevant.append(l)
        infile.close()
        with open("results.txt", "w") as outfile:
            outfile.writelines(relevant)

try:
  os.remove("harness_results.txt")
  os.remove("results.txt")
except:
  print("Files already missing")
execute_tests("CC")
prettify()

                        

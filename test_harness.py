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
  for npes in [2,3,4,5,6]:
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
    for npes in range(5, 7):
      commands.append("oshrun -np " + str(npes) + " ./" + kernel + " -f ../gapbs/gapbs/benchmark/graphs/raw/" + g + " -n 3 -v")
      commands.append("../gapbs/gapbs/src/" + kernel + " -f ../gapbs/gapbs/benchmark/graphs/raw/" + g + " -n 3 -v")
  return commands

def execute_tests(kernel, variety, commands):
  with open (kernel+"harness_results.txt", "a") as out: 
      for c in commands:
          out.write(c+"\n")
          output = subprocess.call(c, stdout=out, shell=True)
          #time.sleep(10)              # I'd use wait, but does subprocess.call spawn the call as a child of this parent?
          print(c + " | Output code: " + str(output))
      out.close()

def prettify(commands, kernel):
    i = 0
    temp = []
    for c in commands:
      if i % 2 == 0: temp.append(c) 
      i+=1
    with open(kernel+"harness_results.txt", "r") as infile: 
        #j = 0
        #relevant = []
        #lines = infile.readlines()
        #for i in range(len(lines)):
        #  l = lines[i]
        #  if j % 3 == 0:        # Assumes all tests run 3 trials
        #    relevant.append(commands[j/3])
        #  if (l.find("Verification:") > -1 and l.find("PE") == -1):
        #    relevant.append(l)
        #    j += 1
        relevant = []
        for l in infile.readlines():
          if (l.find("Verification:") > -1 and l.find("PE") == -1):
            relevant.append(l)
        infile.close()
        print("Trials:", len(relevant), " Commands:", len(temp))
        with open(kernel+"results.txt", "w") as outfile:
          i = 0
          for c in temp:
            outfile.write(c+"\n")
            outfile.writelines(relevant[i:i+3])
            i += 3

try:
  os.remove("harness_results.txt")
  os.remove("results.txt")
except:
  print("Files already missing")
kernel = sys.argv[1]
variety = sys.argv[2]
if (variety == "synthetic"):
  commands = generate_synthetic_commands(kernel)
elif (variety == "real"):
  commands = generate_real_commands(kernel)
elif (variety == "el"):
  commands = generate_el_commands(kernel)
else:
    commands = generate_el_commands(kernel)
    commands = commands + generate_synthetic_commands(kernel)
    commands = commands + generate_real_commands(kernel)
execute_tests(kernel, variety, commands)
prettify(commands, kernel)

                        

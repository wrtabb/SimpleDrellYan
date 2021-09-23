#!/usr/bin/env python
  
import sys
import subprocess

if (len(sys.argv)!=2) :
    print("ERROR: this script expects exactly 1 argument")
    print("Variables: "+sys.argv[1])
    sys.exit(0)

fileName    = sys.argv[1]

print("Input file:  "+fileName)

scriptSpecs = 'analyzeData.C+("'+ str(fileName)+'")'
print(scriptSpecs)
rootCommand = ['root']
rootCommand.append('-l')
rootCommand.append('-b')
rootCommand.append('-q')
rootCommand.append(scriptSpecs)

(out,err) = subprocess.Popen(rootCommand,stdout=subprocess.PIPE).communicate()
print(out)
print("Any errors?")
print(err)

from timeit import default_timer as timer
from subprocess import call
from sys import argv
from os import remove
from os import path

name = "mandelbrot17"
if len(argv) > 1:
  name = argv[1]
fullPbmFileName = "{0}.pbm".format(name)
if path.exists(fullPbmFileName):
    remove(fullPbmFileName)
start = timer()
call("Release/{0} > {1}".format(name, fullPbmFileName) , shell=True)
end = timer()
print(end - start)

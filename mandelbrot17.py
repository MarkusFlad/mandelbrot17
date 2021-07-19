from timeit import default_timer as timer
from subprocess import call
from sys import argv

name = "mandelbrot17"
if len(argv) > 1:
  name = argv[1]
start = timer()
call("Release/{0} > {0}.pbm".format(name) , shell=True)
end = timer()
print(end - start)

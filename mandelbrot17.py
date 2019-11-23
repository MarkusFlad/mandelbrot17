from timeit import default_timer as timer
from subprocess import call

start = timer()
call("Release/mandelbrot17 > mandelbrot17.pbm", shell=True)
end = timer()
print(end - start)

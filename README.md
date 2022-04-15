# Mesh-smoothing
This program smoothes the 3d model file obj provided by the user using different smoothing techniques. After running the make command a smoothing executable is generated. Following commands can be run in order to use different smoothing techniques.

1. Simple laplacian smoothing

./smoothing (input obj file) (output obj file name) (lambda or stepsize) (dt) (iterations)

e.g.- ./smoothing bunny.obj out.obj 1 1 50


2. Using implicit integration in laplacian

./ smoothing -i (eps- small error threshold) (input obj file) (output obj file name) (lambda or stepsize) (dt) (iterations)

e.g.- ./smoothing -i 10e-5 bunny.obj out.obj 1 1 50


3. Global volume preservation

./smoothing -v  (eps- small error threshold) (input obj file) (output obj file name) (lambda or stepsize) (dt) (iterations)

e.g.- ./smoothing -v 10e-5 bunny.obj out.obj 1 1 50


4. Loop subdivision (It takes a lot of time to execute because it involves subdivision and then smoothing)

./smoothing -s (number of subdivision passes) (input obj file) (output obj file name) (lambda or stepsize) (dt) (iterations)

e.g.- ./smoothing -s 1 bunny.obj out.obj 1 1 50


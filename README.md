# GravSim
A simple 3D gravity simulator, written in c/FreeGLUT.

Controls:
  w/a/s/d/q/e = rotate viewfield.
  z/c = zoom viewfield.
  b/v = increment/decrement number of bodies.
  t = enable/disable tracing.
  k/l = save/load file.
  r = restart simulation.
 
Compile: gcc grav.c -o grav -lm -lglut -lGLU -lGL
Run: ./grav

Issues:
  i) Collisions are not accurate in Euclidian space.
  ii) Gravity simulation is bounded between [-1,1].

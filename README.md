# create_bed
create layer by layer particles stratifying desired volume fraction. 

## idea 
assume desire diameter is 1, we can fisrt generate particles with a non-overlapping distance of 0.8, and using simple collision model to increase the distance between particles (move_par_2D) , as well as openmp for acceleration.

## in the code
create 10 layer of particles, each layer at the same height, and the volume fraction is 50%.



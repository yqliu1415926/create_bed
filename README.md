# create_bed
create layer by layer particles stratifying 50% volume fraction. 

## idea 
assume desire diameter is 1, we can fisrt generate particles with a non-overlapping distance of 0.8, and using simple collision model to increase the distance between particles (move_par_2D) , as well as openmp for acceleration.



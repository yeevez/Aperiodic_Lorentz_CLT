matrices = {[1,1,2;1,2,3;2,3,6],[1,2,1;2,5,4;1,4,6],[1,2,2;2,5,6;2,6,9],[1,2,2;2,5,5;2,5,6],[1,1,1;1,2,3;1,3,6]};
%Originally 10^6
bounces = 10^6;
trials = 1;
step = 1;
sideSize =     0.2655;
lowerDim = 2;
n_existing_matrix_samples = 0;
simResult = simulationHypercube(matrices, bounces, trials, step, sideSize, lowerDim, n_existing_matrix_samples);
while (simResult > 0)
     simResult = simulationHypercube(matrices, bounces, trials, step, simResult*0.99, lowerDim, n_existing_matrix_samples);
end
disp("Simulation ended. Final radius of simulation:")
if simResult< 0
    disp(radius)
else
    disp(simResult*0.99)
end


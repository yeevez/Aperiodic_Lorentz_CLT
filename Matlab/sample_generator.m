%matrices = {[1,1,2;1,2,3;2,3,6],[1,2,1;2,5,4;1,4,6],[1,2,2;2,5,6;2,6,9],[1,2,2;2,5,5;2,5,6],[1,1,1;1,2,3;1,3,6]};
matrices = {[1,2,1,1,2;2,5,4,3,5;1,4,6,5,6;1,3,5,7,9;2,5,6,9,14]};
%Originally 10^6
bounces = 20;
trials = 1;
step = 1;
%radius = 0.2320*0.99;
%radius = .1050;
radius = 0.0850;
outdim = 3;

simResult = simulation(matrices, bounces, trials, step, radius,0,outdim);
while (simResult > 0)
     simResult = simulation(matrices, bounces, trials, step, simResult*0.99);
end
disp("Simulation ended. Final radius of simulation:")
disp(simResult)
if simResult< 0
    disp(radius)
else
    disp(simResult*0.99)
end

load('scatter_samples_matrix1.mat', 'paths')
hold on
%plot3(paths(1,:),paths(2,:),paths(3,:))
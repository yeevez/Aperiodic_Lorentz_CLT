%matrices = {[1,1,2;1,2,3;2,3,6],[1,2,1;2,5,4;1,4,6],[1,2,2;2,5,6;2,6,9],[1,2,2;2,5,5;2,5,6],[1,1,1;1,2,3;1,3,6]};

matrices = {[1,2,1,1,2;2,5,4,3,5;1,4,6,5,6;1,3,5,7,9;2,5,6,9,14], 
    [1,3,2,4,3;3,10,7,13,11;2,7,6,10,12;,4,13,10,19,20;3,11,12,20,34],
    [1,4,2,4,2;4,17,10,18,12;2,10,9,16,13;4,18,16,37,23;2,12,13,23,31],
    [1,1,1,3,3;1,2,3,5,5;1,3,6,9,8;3,5,9,18,16;3,5,8,16,16],
    [1,2,4,4,3;2,5,11,11,9;4,11,26,26,23;4,11,26,27,27;3,9,23,27,39]};
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

%% Function saves the time vector and coordinates of each node at each timestep
%%Input : time - vector storing the timesteps; displ - RxNx2xT matrix storing coordinates of each
%%node at each timepoint; filename - name of ".mat" file where other inputs are to be stored. 
%%Output: stores raw simulations results to "filename.mat"

%%Note: function is needed when saving files inside a parfor. for some
%%reason the regular save function does not work inside a parfor loop

%%Author : Tapan Goel
%%Modified: Feb 21, 2022


function parsave(filename,time,displ)
save(filename,'time','displ');
end

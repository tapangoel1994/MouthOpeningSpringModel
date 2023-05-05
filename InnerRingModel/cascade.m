%% Author: Tapan Goel, Feb 23 2022

%%Script takes the name of the simulation file as input
%%and calculates the average compression of the radial springs in each ring as a
%%function of time and plots them.
%%The objective is to see if there is a cascade of compressions from the
%%outerrings to the inner rings.

%%Input: simulationfilename and folder containing a RxNx2xT matrix displ. displ(r,n,i,t) contains the i-th coordinate
%%of the (r,n)th node at timestep t
%%Input: array time containing the actual time at each timestep in seconds.

%%Output: RxT matrix ringcompression, containing the azimuthally averaged extension (+) or
%%compression (-) of radial springs in each ring at each point in time.
%%ringcompression(r,t) gives the average compression of radial springs
%%connnecting the r-th and r+1-th ring at time t.
%%Output: Plot of the average compression of each ring as a function


function ringcompression = cascade(inputfilefolder,inputfilename)

load([inputfilefolder '\' inputfilename]);
init_radial_dist = squeeze(displ(2:end,end,1,1)-displ(1:end-1,end,1,1));   %%The Nth vertex is along the +ve x-axis based on how the nodes are initialized in the simulations

springvec = squeeze( displ(2:end,:,:,:) - displ(1:end-1,:,:,:) );

springdist = squeeze(sqrt(springvec(:,:,1,:).*springvec(:,:,1,:) + springvec(:,:,2,:).*springvec(:,:,2,:)));

meanspringdist = squeeze(mean(springdist,2));

stderrorspringdist = squeeze( std(springdist,0,2) )./sqrt(size(displ,2));

ringcompression = meanspringdist - init_radial_dist;

fig = figure('visible','off');
for i = 1:1:size(ringcompression,1)
    
    errorbar(time(1:10:end),ringcompression(i,1:10:end),stderrorspringdist(i,1:10:end),'-','LineWidth',1);
    hold on;
end
hold off;
xlabel('Time(s)');
ylabel('Extension (um)');
legend('Ring 1','Ring 2','Ring 3','Ring 4','Ring 5');
title(inputfilename);
figfilename = [inputfilefolder '\' erase(inputfilename,'.mat') '-Cascade.fig'];
pngfilename = [inputfilefolder '\' erase(inputfilename,'.mat') '-Cascade.png'];

saveas(fig,figfilename);
saveas(fig,pngfilename);

close all;

end
 
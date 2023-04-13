%% This script sets the parameters for the simulations defined by the differential equations provided in the accompanying manuscript. 
%%Keeping all othe parameters fixed, this script varies
%%the strength of the active force (f0) and the arrival
%%rate of the active force (1/tau). For each set of parameter values (f0,tau),
%%the script runs the simulation 5 times, for different realizations of the
%%sequence of stochastic forces. For each run, the script generates a .mat
%%file containing the coordinates of each node in the network at each time
%%point.The mat files are stored in a user specificed folder for further
%%processing
%%Nodes (R=5,N=2) and (R=5,N=6) are pulled with constant force "microasp_force" (added constant radial external forcing on the stationary nodes in this ring.)

%% System Geometry parameter
radius_central_ring = 10;%initial radius of innermost ring (in micron) - about 1 cell radius
hypostome_radius = 200; %radius of the hypostome (in micron) - from images in the mouth paper

R = 4;%Number of movable rings
N = 8;%Number of nodes per ring

radial_dist = (hypostome_radius-radius_central_ring)/(R); % the immovable ring sets tha radius of the hypostome


%% System mechanics parameters
k0 = 100; %E*h in gram-(micron/s)^2 or nN/um. E = 5kPa (carter et al); h ~ 20 micron 
l0 = 20; %cross over length in microns
gamma0 = k0/(l0^2); %coefficient of non-linear term. Units of nN/um^3
n = 4; %exponent of non-linear term in energy function
mu = 1000; %viscosity = relaxation_time*E*h (Carter et al). relaxation time ~ 10s. in units of nN-s/micron
microasp_force = 5000; %external force from microaspiration in nano-newtons on nodes (R+1),2 and (R+1),6

%% Simulation parameters
dt = .02; %Time step size in seconds (set as 2/5 of sigma)
T = 25; %Simulation duration in seconds - from Carter et al, d = 4.2 secs on average for ectoderm
sigma = .05; %width of 'delta' function in seconds - (HWHM in fig 2D) from Dupre & Yuste, Current Biology, 2017

%% Changing parameters
forcing_strength = logspace(1,6,21); %10*2.^(0:20); % f0: nN
tau = logspace(-1,1,7);%.05*2.^(0:7);  % average arrival time of contractions in seconds.

%% Initialize vertex positions
for r = 1:1:(R+1)
    for theta = 1:1:N
        r0(r,theta,:) = [(radius_central_ring + radial_dist*(r-1))*cosd(theta*360/N) (radius_central_ring + radial_dist*(r-1))*sind(theta*360/N)];
    end
end

%% Create folder to store .mat files    
foldername = ['IndepedentOuterRing' num2str(N) 'Node' num2str(R) 'RingMicroAsp5000nN'];
mkdir(foldername);    



%% Start parameter sweep

rng(1); %% Seed Random Number Generator

poolobj = parpool(16);
tic
parfor ii = 1:1:length(forcing_strength)%[.001 .005 .01 .05 .1 .5 1 5 10]
    for j = 1:1:length(tau) %[.01 .05 .1 .5 1 5 10]; //average time between contractile forcing
        for k = 1:1:length(gamma0)
            for ensamble = 1:1:5  
                
            %% Generate iid poisson force magnitudes
            poissonforce = poissonforcegenerator(R,N,T,dt,tau(j),sigma,forcing_strength(ii));
            %%Only outermost ring has forces. everything else is zero.
            poissonforce(1:end-1,:,:) = 0;
            
            
            %% Initialize vector of radial and azimuthal spring constants, and natural length of azimuthal springs.
            k_rad = k0*ones(R,1);
            k_azim = k0*ones(R,1);
            gamma_rad = gamma0(k)*ones(R,1);
            gamma_azim = gamma0(k)*ones(R,1);
            b_azim = 2*( radius_central_ring + (0:1:(R-1))*radial_dist)*sind(180/N); % b_azim(r) gives the natural length of the azimuthal springs in the rth ring.    

            %% Run simulations for given set of parameters
            [time, displ]  = Oscillator_2D(r0,dt,T,mu,k_rad,gamma_rad,k_azim,gamma_azim,b_azim,n,poissonforce,radial_dist,microasp_force);
               
            %% Save Raw Data from the simulation

            filename1 = [foldername '/force=' num2str(forcing_strength(ii))...
                            ',tau=' num2str(tau(j)) ',gam=' num2str(gamma0(k)) ',i=' num2str(ensamble) '.mat'];
            parsave(filename1,time,displ);
            end
        end

    end
end
toc
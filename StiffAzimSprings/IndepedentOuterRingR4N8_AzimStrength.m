%% This script sets the parameters for the simulations defined by the differential equations provided in the accompanying manuscript. 
%%We fix l0 = 30 microns. Keeping all othe parameters fixed, this script
%%varies the stiffness of the azimuthal springs (k_azim),
%%the strength of the active force (f0) and the arrival
%%rate of the active force (1/tau). For each set of parameter values (k_azim/k_rad, f0,tau),
%%the script runs the simulation 5 times, for different realizations of the
%%sequence of stochastic forces. For each run, the script generates a .mat
%%file containing the coordinates of each node in the network at each time
%%point.The mat files are stored in a user specificed folder for further
%%processing

%%Note: when we change k_azim, we also change gamma_azim since we want to
%%keep the cross over length scale the same for the radial and azimuthal
%%directions

%%Author: Tapan Goel
%%Modified: April 12th, 2023


%% System Geometry parameters
radius_central_ring = 10;%initial radius of innermost ring (in micron) - about 1 cell radius
hypostome_radius = 200; %radius of the hypostome (in micron) - from images in the mouth paper

R = 4;%Number of movable rings
N = 8;%Number of nodes per ring

radial_dist = (hypostome_radius-radius_central_ring)/(R); % the immovable ring sets tha radius of the hypostome



%% System mechanics parameters
k0 = 100; %E*h in gram-(micron/s)^2 or nN/um. E = 5kPa (carter et al); h ~ 20 micron 

n = 4; %exponent of non-linear term in energy function
mu = 1000; %viscosity = relaxation_time*E*h (Carter et al). relaxation time ~ 10s. in units of nN-s/micron %Note: this is inverse of the mu defined in SI Table 1 in the manuscript.


%% Simulation parameters
dt = .02; %Time step size in seconds (set as 2/5 of sigma)
T = 25; %Simulation duration in seconds - from Carter et al, mouth opening time = 4.2 secs on average for ectoderm
sigma = .05; %width of 'delta' function in seconds - (HWHM in fig 2D) from Dupre & Yuste, Current Biology, 2017

%% Changing parameters;
forcing_strength = logspace(1,6,21); %10*2.^(0:20); % f0: nN
tau = logspace(-1,1,7);%.05*2.^(0:7);  % average arrival time of contractions in seconds.
l0 = 30; %cross over length scale in microns
gamma0 = k0/(l0^2)%1./(1:1:4).^2; %coefficient of non-linear term assuming the cross over length is 1 cell rad, 2cell rad and so on. Units of nN/um^3


%% Initialize vertex positions
for r = 1:1:(R+1)
    for theta = 1:1:N
        r0(r,theta,:) = [(radius_central_ring + radial_dist*(r-1))*cosd(theta*360/N) (radius_central_ring + radial_dist*(r-1))*sind(theta*360/N)];
    end
end

%% Create folder to store .mat files    
foldername = ['IndepedentOuterRing' num2str(N) 'Node' num2str(R) 'Ring_AzimStrength'];
mkdir(foldername);



%% Start parameter sweep

rng(1); %% Seed Random Number Generator

poolobj = parpool(16); %for parallel computing with 16 cores.
tic
parfor ii = 1:1:length(forcing_strength) % change "parfor" to "for" if not using parallel computing.
    for j = 1:1:length(tau) 
        for azim_strength = [.5 2 4]
            for ensamble = 1:1:5    

                % Generate iid poisson force magnitudes
                poissonforce = poissonforcegenerator(R,N,T,dt,tau(j),sigma,forcing_strength(ii));
                poissonforce(1:end-1,:,:) = 0; %Only outermost ring has forces. everything else is zero
                
                
                % Initialize vector of radial and azimuthal spring constants, and natural length of azimuthal springs.
                k_rad = k0*ones(R,1);
                k_azim = azim_strength*k0*ones(R,1);
                gamma_rad = gamma0*ones(R,1);
                gamma_azim = azim_strength*gamma0*ones(R,1);
                b_azim = 2*( radius_central_ring + (0:1:(R-1))*radial_dist)*sind(180/N); % b_azim(r) gives the natural length of the azimuthal springs in the rth ring.    
    
                % Run simulations for given set of parameters
                [time, displ]  = Oscillator_2D(r0,dt,T,mu,k_rad,gamma_rad,k_azim,gamma_azim,b_azim,n,poissonforce,radial_dist);
                   
                % Save time vector (time) and coordinates of all nodes at all times (displ) from the simulation
    
                filename1 = [foldername '/force=' num2str(forcing_strength(ii))...
                                ',tau=' num2str(tau(j)) ',gam=' num2str(gamma0) ',azim_strength=' num2str(azim_strength) ',i=' num2str(ensamble) '.mat'];
                parsave(filename1,time,displ);
            end
        end

    end
end

toc

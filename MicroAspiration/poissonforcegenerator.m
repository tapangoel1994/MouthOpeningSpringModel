%% Function generates the magnitude of indepedent poissonforces at each node at each time point
%% For each node, a poisson process in time is created, with average arrival rate 1/tau to generate a sequence of delta function forces.
%% The delta forces are converted to guassians with time width sigma and force = forcing_strength for each kick and added up to get a continuous time series.

%%Input : Number of movable rings (R), Number of nodes per ring (N), Total
%%sim duration (T), timestepsize (dt), average time between tugs (tau),
%%width of gaussian (sigma), strength of force (forcing strength)

%%Output: RxNxT matrix of force magnitudes (force direction is always along the
%%direction of the radial spring). forcing(r,n,t) gives the magnitude of
%%external force on node (r,n) at time t.


%%Note: The force on each node is independent from the other nodes.

%%Author : Tapan Goel
%%Modified: Feb 21, 2022



function forcing = poissonforcegenerator(R,N,T,dt,tau,sigma,forcing_strength)
   
    lambda=1/tau; % arrival rate 
    event=zeros(R,N,1+floor(T/dt)); % initialize counters to zeros
    forcing = zeros(R,N,1+floor(T/dt));
    temp=rand(size(event)); % generate a random array with size as "event"
    event(temp<lambda*dt)=1; % event occur if R < lambda*delta
    inds=find(event==1); % getting indices of arrivial
    

    
    
    for r = 1:1:R
        for theta = 1:1:N
            p = event(r,theta,:);
            p = p(:);
            inds = find(p==1);
            f2 = zeros(size(p));
            for i = 1:1:length(f2)
            f2(i) = sum(exp(-(i*dt-inds*dt).^2./(2*sigma*sigma)));

            end
    
            forcing(r,theta,:) = forcing_strength.*f2./sqrt(2*pi);
        end
    end

end

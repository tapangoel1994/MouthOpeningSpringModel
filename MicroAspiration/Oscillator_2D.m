%% Uses an RK4 integrator to integrate the network configuration over time 
%% and output the position of each node at each time point


function [time,x]  = Oscillator_2D(r0,dt,T,mu,k_rad,gamma_rad,k_azim,gamma_azim,b_azim,n,poissonforce,radial_dist,microasp_force)

%% Initial conditions:

timesteps = length(0:dt:T);
p = size(r0);
rows = p(1);
N = p(2);
x = zeros(rows,N,2,timesteps);
x(:,:,:,1) = r0;


%% Iterate over timesteps

for i = 2:1:timesteps
    
    time = dt*(i-2);
       
    k1 = func2D(time,x(:,:,:,i-1),mu,k_rad,gamma_rad,k_azim,gamma_azim,b_azim,poissonforce,dt,n,radial_dist,microasp_force);
    k2 = func2D(time + dt/2, x(:,:,:,i-1) + dt.*k1/2,mu,k_rad,gamma_rad,k_azim,gamma_azim,b_azim, poissonforce,dt,n,radial_dist, microasp_force);
    k3 = func2D(time + dt/2, x(:,:,:,i-1) + dt.*k2/2,mu,k_rad,gamma_rad,k_azim,gamma_azim,b_azim, poissonforce,dt,n,radial_dist, microasp_force);
    k4 = func2D(time + dt, x(:,:,:,i-1) + dt.*k3, mu,k_rad,gamma_rad,k_azim,gamma_azim,b_azim, poissonforce,dt,n,radial_dist, microasp_force);
    
    x(:,:,:,i) = x(:,:,:,i-1) + dt.*(k1 + 2*k2 + 2*k3 + k4)./6;
       
end
 
time = 0:dt:T;

end



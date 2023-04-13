%% Calculates the time derivative of the position vectors of each node in the grid given the current configuration of the network
%%derivative is calculated by finite difference between term slightly above and slightly below
%%r0 is a (R+1)xNx2 matrix that contains the positions of all the nodes in the network. The R+1th 
%%row contains the positions of the fixed boundary nodes.
%%k is the spring constant, gamma is coefficient of non-linearity, n is the exponent of the non linear
%%term, mu is the coefficient of mobility.
%%radial_dist is the radial distance between rings
%%radius_central_ring is the radius of the central ring
%%R,N are the number of movable rings and nodes/ring respectively.
%%dl is the spatial step-size that defines the grid. Note that accuracy of derivative calculation
%%depends on choice of dl
%%tau is the time period of periodic signalling at a node.
%%radial_speed is the speed of the neuronal signal that causes contraction
%%forcing_strength is the strength of the delta function force illicted by the neuronal spike
%%sigma is the width of the neuronal spike

%%Author: Tapan Goel
%%Modified: Feb 21st 2022

function drdt = func2D(t0,r0,mu,k_rad,gamma_rad,k_azim,gamma_azim, b_azim, poissonforce,dt,n,radial_dist, microasp_force)

 
 drdt = zeros(size(r0));  %%Initialize the derivative vector. 

 drdt = (1/mu)*TotalForce(r0,k_rad,gamma_rad,k_azim,gamma_azim, b_azim, poissonforce,t0,dt,n,radial_dist,microasp_force);
    

 
 
end
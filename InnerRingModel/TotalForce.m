%% Calculates the total  = spring+poisson force on every node given the current configuration of the network by direct evaluation of force expression;

%% Inputs: 
%  r0: (R+1)xNx2 matrix containing current location of all nodes. 
%  The (R+1)th row contains the positions of the fixed boundary nodes.
%  k_rad(r), gamma_rad(r) are constants for springs connecting the rth and (r+1)th rings
%  k_azim(r), gamma_azim(r) are azimuthal spring constants for the rth ring.
%  b_azim(r) is the natural spring length in the 'r'th ring.
%  poissonforce(r,theta,N) is the magnitude of the poisson force on the
%  (r,theta) node at timestep N
%  t0 is the time at the current timestep
%  dt is the length of the timestep in the simulation
%  n is the order of the non-linear term in the energy function
%  radial_dist is the radial distance between rings
%  radius_central_ring is the radius of the central ring
%


%% Output:
%  force: RxNx2 vector of total = elastic + poisson forces on nodes. force(r,theta,:) gives
%  the total force on the (r,theta)-th node

%% Note: 
%  1. code uses dirvec to get the direction of the spring force. dirvec
%  returns a zero vector if two nodes are on top of each other.
%  2. poisson force is always radial - not along the springs ***Changed***

%%Author : Tapan Goel, Jan 13 2022
%%Modified: Tapan Goel, Feb 15 2022: PoissonForce is along the springs, not radial.

function force  = TotalForce(r0,k_rad,gamma_rad,k_azim,gamma_azim,b_azim,...
                    poissonforce,t0,dt,n,radial_dist)

p = size(r0); % to output number of rings (including boundary) and #vertices per ring
R = p(1)-1;  % r0 has 1 row more than the number of movable rings
N = p(2);
thetacc = circshift(1:1:N,-1); %indices go counterclockwise by 1 unit;
thetac = circshift(1:1:N,1); %indices go clockwise by 1 unit;

force = NaN(R,N,2);  %Initializes all forces to be NaN so that an error pops up if there is some error in the computation

%% Forces on the nodes in the first ring
    for theta = 1:1:N
    
        %%Calculate extensions of all springs connecting the node
        dr2 = norm(squeeze( r0(1,theta,:)- r0(2,theta,:) ) ) - radial_dist;
        dazimcc = norm(squeeze(r0(1,theta,:))-squeeze(r0(1,thetacc(theta),:))) - b_azim(1);
        dazimc = norm(squeeze(r0(1,theta,:))-squeeze(r0(1,thetac(theta),:))) - b_azim(1);
        
        %%Calculate total force on the node; NB: dirvec returns zero vector
        %%if nodes are on top of each other
        f = [0;0];
        f = f -( k_rad(1)*dr2 + gamma_rad(1)*(dr2^(n-1)) )* dirvec( squeeze(r0(1,theta,:))-squeeze(r0(2,theta,:)) );
        f = f -(k_azim(1)*dazimcc + gamma_azim(1)*(dazimcc^(n-1)))* dirvec(squeeze(r0(1,theta,:))-squeeze(r0(1,thetacc(theta),:)) );
        f = f -(k_azim(1)*dazimc + gamma_azim(1)*(dazimc^(n-1)))* dirvec(squeeze(r0(1,theta,:))-squeeze(r0(1,thetac(theta),:)));

        force(1,theta,:) = f + poissonforce(1,theta,1+floor(t0/dt)).*dirvec(squeeze(r0(2,theta,:))-squeeze(r0(1,theta,:)));
    end

%% Forces on the nodes in the remaining rings    
for r = 2:1:R
    for theta = 1:1:N
        
       
        dr_rlow = norm(squeeze(r0(r,theta,:))-squeeze(r0(r-1,theta,:))) - radial_dist;
        dr_rhigh = norm(squeeze(r0(r,theta,:))-squeeze(r0(r+1,theta,:))) - radial_dist;
        dazimcc = norm(squeeze(r0(r,theta,:))-squeeze(r0(r,thetacc(theta),:))) - b_azim(r);
        dazimc = norm(squeeze(r0(r,theta,:))-squeeze(r0(r,thetac(theta),:))) - b_azim(r);
        
        f = [0;0];
        
        f = -(k_rad(r-1)*dr_rlow + gamma_rad(r-1)*(dr_rlow^(n-1)) )* dirvec( (squeeze(r0(r,theta,:))-squeeze(r0(r-1,theta,:))) );
        f = f -(k_rad(r)*dr_rhigh + gamma_rad(r)* (dr_rhigh^(n-1)))* dirvec( (squeeze(r0(r,theta,:))-squeeze(r0(r+1,theta,:))) );
        
        f = f -(k_azim(r)*dazimcc + gamma_azim(r)*(dazimcc^(n-1)))* dirvec(squeeze(r0(r,theta,:))-squeeze(r0(r,thetacc(theta),:)) );
        f = f -(k_azim(r)*dazimc + gamma_azim(r)*(dazimc^(n-1)))* dirvec(squeeze(r0(r,theta,:))-squeeze(r0(r,thetac(theta),:)));

        force(r,theta,:) = f + poissonforce(r,theta,1+floor(t0/dt)).*dirvec(squeeze(r0(r+1,theta,:))-squeeze(r0(r,theta,:)));    
    end
end



end
 


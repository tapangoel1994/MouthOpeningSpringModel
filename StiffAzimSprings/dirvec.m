%% Function calculates direction vector from a given vector:
%%Input : v, an nx1 vector
%%Output: direction, either v normalized or [0;0] itself if norm(v) == 0

%%Author : Tapan Goel, Jan 13 2022


function direction = dirvec(v)
 
if(norm(v) == 0)
    error('Nodes are on top of each other');
    direction = [0;0];
else
    direction = v/norm(v);
end

function [fx,fy] = compute_fx_fy(rho,rho_obs,grid_x,grid_y,q_T_x,q_T_y)

%% Compute forcing term in x and y directions
% Input:
% rho: density of the swarm
% rho_obs: density of the traces
% grid_x: meshgrid of the x-coordinate
% grid_y: meshgrid of the y-coordinate
% q_T_x: discretized kernel on the x-coordinate
% q_T_y: discretized kernel on the y-coordinate
% Output:
% fx: forcing term along the x-coordinate
% fy: forcing term along the y-coordinate

% Size of the grid
N_xC = size(grid_x,1);
N_yC = size(grid_x,2);

delta_x = grid_x(2,1)-grid_x(1,1);
delta_y = grid_y(1,2)-grid_y(1,1);

%% Computation of the forcing
fx = zeros(size(rho)); % Pre-allocation
fy = zeros(size(rho)); % Pre-allocation
for i = 1:N_xC
    for j = 1:N_yC
       fx(i,j) = rho(i,j)*delta_x*delta_y*sum(squeeze(q_T_x(i,j,:,:)).*(rho+rho_obs),"all"); % This is a convolution
       fy(i,j) = rho(i,j)*delta_x*delta_y*sum(squeeze(q_T_y(i,j,:,:)).*(rho+rho_obs),"all"); % This is a convolution
    end
end
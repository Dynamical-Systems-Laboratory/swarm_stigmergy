function rho_obs = find_rho_obs(f,mesh,L,rho_obs_mpi)
%% Compute the density of the obstacles
% Input: 
% f: Objective forcing term
% mesh: Mesh on the circle
% L: Length-scale of the repulsive kernel
% rho_obs_mpi: Value of rho_obs at -pi (it does not affect the final
% result)
% Output:
% rho_obs: Convolution over the circle at each point of the mesh

% We use a smaller version of f and mesh, as it is periodic and the value at -pi
% is already equal to that at pi
f_d = f(1:end-1);
mesh_d = mesh(1:end-1);

% We extend them from periodicity to be able to compute centered
% finite differences
f_ext = [f_d(end); f_d; f_d(1)];
dx = mesh(2)- mesh(1);

% Compute centered finite differences
df = (f_ext(3:end)-f_ext(1:end-2))/(2*dx);

% Compute rho_obs from deconvolution
rho_obs_d = 0.5*(df-df(1)-1/L^2*cumtrapz(mesh_d,f_d))+rho_obs_mpi;
% Extend rho_obs to become periodic again
rho_obs = [rho_obs_d; rho_obs_d(1)];
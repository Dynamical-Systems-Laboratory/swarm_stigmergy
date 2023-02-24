function count = discretize_fun(f,mesh,N_scale)
%% Discretize the position of the robots
% Inputs:
% f: Density function to be discretized
% mesh: Mesh on the circle
% N_scale: Number of discrete units per unit area of f (that is, assuming that f is unity)
% Output: 
% count: Counts at each of the midpoint elements

dx = mesh(2)-mesh(1); % Element size

area = trapz(mesh,f); % Area of f
f_sc = f/area; % Scaled f to have unit area

fC_sc = (f_sc(1:end-1)+f_sc(2:end))/2; % f at the midpoint elements

N_tot = round(N_scale*area); % Approximation of the number of elements to be deployed

count = round(N_tot*fC_sc*dx); % Computation of the counts at each of the midpoint elements
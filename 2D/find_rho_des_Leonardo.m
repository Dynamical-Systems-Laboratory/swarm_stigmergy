function [rho_des,x_tot,y_tot] = find_rho_des_Leonardo(R_0, A, N_theta, k_x, k_y, N, dt, t_span_fr3, k_y_amp)
%% This function creates the desired density for the 2D case of Leonardo's lion
% Input:
% R_0: Radius of the circle
% A: Amplitude of the circle
% N_theta: Number of angles along the circle
% k_x: Parameter of von Mises distribution
% k_y: Parameter of von Mises distribution
% N: Size of the grid (equal dimension of each side)
% dt: Time interval in the finite volume algorithm
% t_span_fr3: Decompose the time vector in 3 parts
% k_y_amp: Max amplitude of von Mises distribution
% Output:
% rho_des: Desired density over time (X_i, Y_j, t_k)
% x_tot: Parametrized x coordinates of the centroid of the trajectory (flower shape)
% y_tot: Parametrized y coordinates of the centroid of the trajectory (flower shape)

%% Preliminary steps
theta_tot = linspace(0,2*pi,N_theta)'; % Meshing angular coordinate

x_f = @(theta) (R_0+A*cos(6*theta)).*cos(theta); % Function of the wavy circle in x-direction
y_f = @(theta) (R_0+A*cos(6*theta)).*sin(theta); % Function of the wavy circle in y-direction

x_tot = x_f(theta_tot); % Compute centroid trajectory
y_tot = y_f(theta_tot); % Compute centroid trajectory


%% Meshing
mesh = linspace(-pi,pi,N); % Create mesh

[grid_x,grid_y] = meshgrid(mesh,mesh); % Create meshgrid
grid_x = grid_x';
grid_y = grid_y';

%% Time change
t_vec = [0:dt:t_span_fr3*3]'; % Time vector
N_t = length(t_vec); % Length of the time span

t_ind_p = find(t_vec==t_span_fr3); % Find index of the time vector for the first transition
t_ind_r = 2*t_ind_p; % Find index of the time vector for the second transition

% Modulate k_y over time to allow the swarm to "open" and "close"
k_y_vec = [ones(t_ind_p,1)*k_y; ones(t_ind_p,1)*k_y+k_y_amp*[[1:t_ind_p/2]'; t_ind_p/2; [t_ind_p/2:-1:1]']/(t_ind_p/2); ones(t_ind_p,1)*k_y];

omega = pi/t_span_fr3; % Speed

% Angle of the centroid of the swarm (mu for von Mises)
theta = zeros(size(t_vec)); %Preallocation
theta_0 = -pi/2; % Initial bias
theta(1:t_ind_p) = theta_0+omega*t_vec(1:t_ind_p); % First segment
theta(t_ind_p+1:t_ind_r) = theta(t_ind_p); % Second segment (stop)
theta(t_ind_r+1:end) = theta(t_ind_p)+omega*(t_vec(t_ind_r+1:end)-t_vec(t_ind_r)); % Last segment

%% Compute desired density profile
rho_des = zeros(N,N,N_t); % Preallocation
for i = 1:length(t_vec)
    rho_des(:,:,i) = exp(k_x*cos(grid_x-x_f(theta(i)))).*exp(k_y_vec(i)*cos(grid_y-y_f(theta(i))))/(4*pi^2*besseli(0,k_x)*besseli(0,k_y_vec(i)));
end
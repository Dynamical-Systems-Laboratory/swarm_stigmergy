clear all; close all; clc;

%% Parameters
N = 5000; % Number of robots
M = 1; % Total mass
n_tot = 250; % Number of mesh points to be used on the circle (rule of thumb: approx N/20 is enough)
L = 1; % Interaction scale for repulsive potential
k = 2; % Parameter of von Mises distribution
mu = -1; % Parameter of von Mises distribution
dt = 1.e-1; % Time interval in the simulation (in this case it does not matter as we use ode45 from Matlab)
t_span = 100; % Time span of the simulation
sigma = 2; % Parameter of the kernel for density estimation
scale_factor = 1; % Scale factor in the continuification (needed if densities have different integrals)

m = M/N; % Mass of each individual robot

%% Repulsive kernel
q = @(x) 1/(exp(2*pi/L)-1)*sign(x).*(exp((2*pi-abs(x))/L)-exp(abs(x)/L));

%% Meshing and distributions
mesh = linspace(-pi,pi,n_tot)'; % Mesh on the circle
dx = mesh(2)-mesh(1); % Element size
meshC = mesh(1:end-1)+dx/2; % Mesh of the element midpoints
t_vec = [0:dt:t_span]';

rho_obj = M*exp(k*cos(mesh-mu))./(2*pi*besseli(0,k)); % Objective density
rho_obj = rho_obj/trapz(mesh,rho_obj); % Scale objective density to integrate to 1
rho_obj_t = zeros(size(rho_obj)); % Time derivative of objective density

%% Preliminary computations
f_obj = compute_f_obj_rep(rho_obj,rho_obj_t,L,mesh); % Objective forcing term f
rho_obs_mpi = 0; % Set the value of rho_obs at -pi (it does not affect the final result)
rho_obs = find_rho_obs(f_obj,mesh,L,rho_obs_mpi); % Compute the density of the obstacles
rho_obs = rho_obs-min(rho_obs); % Translate the density of the obstacles so that it is pointwise not negative
f = compute_conv_rep(rho_obs,L,mesh); % Compute the forcing term from the distribution of the obstacles

%% Discretization of the obstacles
counts_obs = discretize_fun(rho_obs,mesh,N); % Compute the counts of the number of obstacles on the mesh
N_obs = sum(counts_obs); % Number of obstacles

% Find the position of the obstacles
pos_obs = zeros(N_obs,1); % Pre-allocation
i = 1;
j = 1;
while i<N_obs
    pos_obs(i:i+counts_obs(j)-1) = mesh(j);
    i = i+counts_obs(j);
    j = j+1;
end

%% Computations
x0 = linspace(-pi,pi-eps,N)'; % Initial position (uniform) of the robots

[t_vec,x_sol] = ode45(@(t,x) discrete_eqs(x, pos_obs, q, m,t), t_vec, x0); %Simulation of the discrete equations
x_sol = wrapToPi(x_sol); % Wrap to [-pi,pi]

%% KL divergence
KL_div = zeros(size(t_vec));
for i = 1:length(t_vec)-1
    density = continuize_fun(x_sol(i,:)',mesh,sigma,scale_factor); % Continuify the density of the obstacles
    KL_div(i) = trapz(mesh,density.*log(density./rho_obj)); % Computation of the KL divergence
end    

%% Plots
set(groot, 'defaultAxesFontName', 'Arial')
set(groot, 'defaultTextFontName', 'Arial')

% Final density
figure
density = continuize_fun(x_sol(end,:)',mesh,sigma,scale_factor);
plot(mesh,density)
hold on
plot(mesh, rho_obj)

% KL divergence plot
figure
plot(t_vec,KL_div)
xlabel('Time','FontName', 'Arial')
ylabel('KL Divergence','FontName', 'Arial')
set(gca,'FontSize',14)

%% Export
exp_KL = [t_vec, KL_div];

save('all_data_discrete.mat')
writematrix(exp_KL,'Static_KL.txt','Delimiter','tab');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Function for evolution of each unit in the swarm
function x_dot = discrete_eqs(x, x_obs, q, m, t)
    % Function for evolution of each unit in the swarm
    % Input: 
    % x: ode45 state (position of the robots)
    % x_obs: Position of the obstacles
    % q: Repulsive kernel
    % m: Mass of each robot
    % t: Time
    % Output:
    % x_dot: Time derivative of the state (velocity of the robots)
    x = wrapToPi(x); % Wrap to [-pi,pi]
    N = length(x); % Number of robots
    N_obs = length(x_obs); % Number of obstacles
    x_dot = m*(sum(q(repmat(x,[1,N])-repmat(x',[N,1])),2)+sum(q(repmat(x,[1,N_obs])-repmat(x_obs',[N,1])),2)); % Time derivative
    t % Print time to see evolution
end
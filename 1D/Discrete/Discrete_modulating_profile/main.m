clear all; close all; clc;

%% Parameters
N = 5000; % Number of robots
M = 1; % Total mass
n_tot = 250; % Number of mesh points to be used on the circle (rule of thumb: approx N/20 is enough)
L = 1; % Interaction scale for repulsive potential
k = 2; % Parameter of von Mises distribution
mu = -1; % Parameter of von Mises distribution
v = 0.05; % Speed of translation
Omega = 0.1; % Oscillations' natural frequency
A = 0.5; % Amplitude of oscillations
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

t_ind = find(t_vec==0); % You can change the time the objective density starts changing

% Modulation of k
k_mod = 2+A*sin(Omega*(t_vec(t_ind:end)'-t_vec(t_ind))-pi);

% Construction of the objective density
rho_obj = M*exp([k_mod k_mod(end)].*cos(mesh-mu+v*[t_vec;t_vec(end)+dt]'))./(2*pi*besseli(0,[k_mod k_mod(end)]));
% Time derivative of the objective density
rho_obj_t = (rho_obj(:,2:end)-rho_obj(:,1:end-1))/dt;

%% Preliminary computations
% Preallocation
f_obj = zeros(length(mesh),length(t_vec));
rho_obs = zeros(length(mesh),length(t_vec));
counts_obs = zeros(length(mesh)-1,length(t_vec));
N_obs = zeros(length(t_vec),1);

rho_obs_mpi = 0; % Set the value of rho_obs at -pi (it does not affect the final result)
pos_obs = cell(length(t_vec),1); % Position of the obstacles over time
for i = 1:length(t_vec)
    f_obj(:,i) = compute_f_obj_rep(rho_obj(:,i),rho_obj_t(:,i),L,mesh); % Objective forcing term f
    rho_obs(:,i) = find_rho_obs(f_obj(:,i),mesh,L,rho_obs_mpi); % Compute the density of the obstacles
    rho_obs(:,i) = rho_obs(:,i)-min(rho_obs(:,i)); % Translate the density of the obstacles so that it is pointwise not negative
    counts_obs(:,i) = discretize_fun(rho_obs(:,i),mesh,N); % Compute the counts of the number of obstacles on the mesh
    N_obs(i) = sum(counts_obs(:,i)); % Number of obstacles
    pos_obs{i} = find_pos(counts_obs(:,i),N_obs(i),mesh); % Find the position of the obstacles
end

%% Computations
x0 = linspace(-pi,pi-eps,N)'; % Initial position (uniform) of the robots

[t_vec,x_sol] = ode45(@(t,x) discrete_eqs(x, pos_obs, q, m,t,t_vec), t_vec, x0); %Simulation of the discrete equations
x_sol = wrapToPi(x_sol); % Wrap to [-pi,pi]

%% KL divergence
KL_div = zeros(size(t_vec));
for i = 1:length(t_vec)-1
    density = continuize_fun(x_sol(i,:)',mesh,sigma,scale_factor); % Continuify the density of the obstacles
    KL_div(i) = trapz(mesh,density.*log(density./rho_obj(:,i))); % Computation of the KL divergence
end    

%% Plots
set(groot, 'defaultAxesFontName', 'Arial')
set(groot, 'defaultTextFontName', 'Arial')

% Final density
figure
density = continuize_fun(x_sol(end,:)',mesh,2, 1);
plot(mesh,density)
hold on
plot(mesh, rho_obj(:,end))

% KL divergence plot 
figure
plot(t_vec,KL_div)
xlabel('Time','FontName', 'Arial')
ylabel('KL Divergence','FontName', 'Arial')
set(gca,'FontSize',14)

%% Export
exp_KL = [t_vec, KL_div];

writematrix(exp_KL,'Dynamic_modulating_KL.txt','Delimiter','tab');
save('all_data_discrete.mat')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Function for evolution of each unit in the swarm
function x_dot = discrete_eqs(x, x_obs, q, m, t,t_vec)
    % Function for evolution of each unit in the swarm
    % Input:
    % x: ode45 state (position of the robots)
    % x_obs: Position of the obstacles
    % q: Repulsive kernel
    % m: Mass of each robot
    % t: Time
    % t_vec: Time of vectors for the obstacles
    % Output:
    % x_dot: Time derivative of the state (velocity of the robots)
    closest = interp1(t_vec,t_vec,t,'nearest','extrap'); % Find closest time in the time of vectors
    x_obs = x_obs{t_vec==closest}; % Set the position of the obstacles equal to the closest time
    x = wrapToPi(x); % Wrap to [-pi,pi]
    N = length(x); % Number of robots
    N_obs = length(x_obs); % Number of obstacles
    x_dot = m*(sum(q(repmat(x,[1,N])-repmat(x',[N,1])),2)+sum(q(repmat(x,[1,N_obs])-repmat(x_obs',[N,1])),2)); % Time derivative
    t % Print time to see evolution
end

clear all; close all; clc;

%% Parameters
N = 501; % Number of points
L = 1; % Interaction scale for repulsive potential
k = 2; % Parameter of von Mises distribution
mu = -1; % Parameter of von Mises distribution
dt = 1.e-1; % Time interval in the finite volume algorithm
t_span = 100; % Time span of the simulation

%% Meshing and distributions
mesh = linspace(-pi,pi,N)'; % Mesh on the circle
t_vec = [0:dt:t_span]'; % Vector of time instances

rho_in = 1/(2*pi)*ones(size(mesh)); % Initial density (uniform)
rho_obj = exp(k*cos(mesh-mu))./(2*pi*besseli(0,k)); % Objective density
% rho_obj = rho_obj/trapz(mesh,rho_obj); %You can add this line to ensure the distribution has unitary integral
rho_obj_t = zeros(size(rho_obj)); % Time derivative of objective density

%% Preliminary computations
f_obj = compute_f_obj_rep(rho_obj,rho_obj_t,L,mesh); % Objective forcing term f
rho_obs_mpi = 0; % Set the value of rho_obs at -pi (it does not affect the final result)
rho_obs = find_rho_obs(f_obj,mesh,L,rho_obs_mpi); % Compute the density of the obstacles
f = compute_conv_rep(rho_obs,L,mesh); % Compute the forcing term from the distribution of the obstacles

%% Computations of the density evolution
rho = zeros(N,length(t_vec)); % Preallocation
rho(:,1) = rho_in; % Initial density
for i = 1:length(t_vec)-1
   rho(:,i+1) = lax_friedrichs(rho(:,i), dt, mesh, L, f); % Finite volume evolution
end

%% Error and KL divergence
err = abs(rho-repmat(rho_obj,[1,length(t_vec)])); % Error

KL_div = zeros(size(t_vec)); % Preallocation
for i = 1:length(t_vec)
    KL_div(i) = trapz(mesh,rho(:,i).*log(rho(:,i)./rho_obj)); % Computation of the KL divergence
end    

%% Plots
% Set fonts
set(groot, 'defaultAxesFontName', 'Arial')
set(groot, 'defaultTextFontName', 'Arial')

% Final density
figure
plot(mesh,rho(:,end))
hold on
plot(mesh,rho_obj,'--','LineWidth',1.5)
axis([-pi pi 0 0.6])
xticks([-pi -pi/2 0 pi/2 pi])
xticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'})
xlabel('Position','FontName', 'Arial')
ylabel('Density of the swarm','FontName', 'Arial')
set(gca,'FontSize',14)

% Obstacle density
figure
plot(mesh,rho_obs)
axis([-pi pi -0.5 0.1])
xticks([-pi -pi/2 0 pi/2 pi])
xticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'})
xlabel('Position','FontName', 'Arial')
ylabel('Density of obstacles','FontName', 'Arial')
set(gca,'FontSize',22)

% Error plot
[grid_x,grid_t] = meshgrid(mesh,t_vec); % Grid for the plot
figure
s = surf(grid_t,grid_x,err');
s.EdgeColor = 'none';
hold on
view(2)
xlabel('Time','FontName', 'Arial')
ylabel('Position','FontName', 'Arial')
colorbar
axis([0 t_span -pi pi])
caxis([0 0.35])
yticks([-pi -pi/2 0 pi/2 pi])
yticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'})
set(gca,'FontSize',14)

% KL divergence plot
figure
plot(t_vec,KL_div)
xlabel('Time','FontName', 'Arial')
ylabel('KL Divergence','FontName', 'Arial')
set(gca,'FontSize',22)

% Export video
v = VideoWriter('S1.avi');

open(v);
figure
for i = 1:5:size(rho,2) % 5x speed
    plot(mesh,rho(:,i))
    hold on
    plot(mesh,rho_obj,'--')
    set(gca,'FontSize',14)
    xlabel('Position','FontName', 'Arial')
    ylabel('Density of the swarm','FontName', 'Arial')
    axis([-pi pi 0 0.6])
    xticks([-pi -pi/2 0 pi/2 pi])
    xticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'})
    hold off
    pause(0.01)
    frame = getframe(gcf);
    writeVideo(v,frame);
end
close(v);

%% Export data as txt files
exp_profiles = [mesh, rho(:,end),rho_obj];
exp_profile_obs = [mesh, rho_obs];
exp_KL_div = [t_vec,KL_div];

% Downsampling for the colorplot
downsample_x = downsample(mesh,5);
ind_x = downsample([1:length(mesh)]',5);
downsample_t = downsample(t_vec,5);
ind_t = downsample([1:length(t_vec)]',5);
N_x_down = length(downsample_x);
N_t_down = length(downsample_t);

grid_x_down = grid_x(ind_t, ind_x);
grid_t_down = grid_t(ind_t, ind_x);
err_down = err(ind_x, ind_t);

exp_colorplot = [reshape(grid_t_down,[N_x_down*N_t_down,1]) reshape(grid_x_down,[N_x_down*N_t_down,1]) reshape(err_down',[N_x_down*N_t_down,1])];

% Export
writematrix(exp_profiles,'Static_profiles.txt','Delimiter','tab');
writematrix(exp_profile_obs,'Static_profile_obs.txt','Delimiter','tab');
writematrix(exp_KL_div,'Static_KL_div.txt','Delimiter','tab');
writematrix(exp_colorplot,'Static_colorplot.txt','Delimiter','tab');

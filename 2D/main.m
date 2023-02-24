clear all; close all; clc;

%% Load the discretized repulsive kernel over the grid
load('qTs.mat'); % You need to recompute qTs if you change the parameters

%% Parameters
N_xC = 30; % Number of cells along the x-axis
N_yC = 30; % Number of cells along the y-axis
L_x = 2*pi; % Length along the x-axis
L_y = 2*pi; % Length along the y-axis
dt = 1.e-1; % Time interval in the finite volume algorithm
t_span = 300; % Time span of the simulation
R_0 = 2.5; % Radius of the circle
A = 0.3; % Amplitude of the circle
N_theta = 500; % Number of angles along the circle
k_x = 2; % Parameter of von Mises distribution
k_y = 2; % Parameter of von Mises distribution
k_y_amp = 3; % Max amplitude of von Mises distribution

%% Meshing
N_bar_x = N_xC+1; % Number of mesh points along the x-axis
N_bar_y = N_yC+1; % Number of mesh points along the y-axis
mesh_x = linspace(-L_x/2,L_x/2,N_bar_x)'; % Mesh along the x-axis
mesh_y = linspace(-L_y/2,L_y/2,N_bar_y)'; % Mesh along the y-axis
mesh_xC = mesh_x(1:end-1)+(mesh_x(2)-mesh_x(1))/2; % Mesh of the midpoint elements along the x-axis
mesh_yC = mesh_y(1:end-1)+(mesh_y(2)-mesh_y(1))/2; % Mesh of the midpoint elements along the y-axis
delta_x = mesh_x(2)-mesh_x(1); % Element size along the x-axis
delta_y = mesh_y(2)-mesh_y(1); % Element size along the y-axis
[grid_x,grid_y] = meshgrid(mesh_xC,mesh_yC); % Mesh grid of the midpoints of the cells
[grid_x_ext,grid_y_ext] = meshgrid(mesh_x,mesh_y); % Mesh grid of the mesh
grid_x = grid_x';
grid_y = grid_y';

t_vec = [0:dt:t_span]'; % Time vector used in the finite volume algorithm 
N_t = length(t_vec); % Length of the time vector
t_span_fr3 = t_span/3; % Decompose the time vector in 3 parts

%% Distribution of the density
% Desired density for Leonardo's robot
[rho_des,x_fl,y_fl] = find_rho_des_Leonardo(R_0, A, N_theta, k_x, k_y, N_xC, dt, t_span_fr3, k_y_amp);
% We neglect the time derivative of the desired density (quasi-static problem)
% One can put the time derivative, but it requires increasing the
% resolution in both space and time
rho_des_t = zeros(size(rho_des(:,:,1:end-1))); 

% Initial density (uniform)
rho_in = ones(size(grid_x))/(L_x*L_y);

%% Computation of the density
rho = zeros(N_xC,N_yC,length(t_vec)); % Preallocation
rho(:,:,1) = rho_in; % Initial density
rho_obs0 = -rho_des(:,:,1); % Initial try for rho_obs
for i = 1:length(t_vec)-1
   rho_obs = fsolve(@(x) find_rho_obs(x,rho_des(:,:,i),rho_des_t(:,:,i),grid_x,grid_y,q_T_x,q_T_y,dt),rho_obs0); % Compute the needed density of the obstacles
   rho_obs0 = rho_obs; % Updated initial try for rho_obs
   rho(:,:,i+1) = lax_friedrichs_2D(rho(:,:,i), rho_obs, dt, grid_x, grid_y,q_T_x,q_T_y); % Compute the updated density
   i % To track progress in the loop
end

%% KL divergence
KL_div = zeros(size(t_vec));
for i = 1:length(t_vec)
    KL_div(i) = compute_2D_integral(rho(:,:,i).*log(rho(:,:,i)./rho_des(:,:,i)), delta_x, delta_y);
end    

%% Plots
% Set font
set(groot, 'defaultAxesFontName', 'Arial')
set(groot, 'defaultTextFontName', 'Arial')

% KL Divergence
figure
plot(t_vec,KL_div)
axis([0 t_span 0 2])
xlabel('Time','FontName', 'Arial')
ylabel('KL Divergence','FontName', 'Arial')
set(gca,'FontSize',14)

% Initial - 3D
figure
surf(grid_x,grid_y,rho(:,:,100))
axis([-pi pi -pi pi -0.01 0.5])
zticks([0 0.25 0.5])
xlabel('x','FontName', 'Arial')
ylabel('y','FontName', 'Arial')
zlabel('Density of the swarm','FontName', 'Arial')
set(gca,'FontSize',14)

% Initial - top view
figure
surf(grid_x,grid_y,rho(:,:,100))
view(2)
pbaspect([1 1 1])
hold on
plot3(x_fl,y_fl,ones(length(x_fl)),'r','Linewidth',3)
axis([-pi pi -pi pi -1 1.5])
xticks([-pi -pi/2 0 pi/2 pi])
xticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'})
yticks([-pi -pi/2 0 pi/2 pi])
yticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'})
hold off
xlabel('x','FontName', 'Arial')
ylabel('y','FontName', 'Arial')
set(gca,'FontSize',22)

% Quarter - 3D
figure
surf(grid_x,grid_y,rho(:,:,500))
axis([-pi pi -pi pi -0.01 0.5])
zticks([0 0.25 0.5])
xlabel('x','FontName', 'Arial')
ylabel('y','FontName', 'Arial')
zlabel('Density of the swarm','FontName', 'Arial')
set(gca,'FontSize',14)

% Quarter - top view
figure
surf(grid_x,grid_y,rho(:,:,500))
view(2)
pbaspect([1 1 1])
hold on
plot3(x_fl,y_fl,ones(length(x_fl)),'r','Linewidth',3)
axis([-pi pi -pi pi -1 1.5])
xticks([-pi -pi/2 0 pi/2 pi])
xticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'})
yticks([-pi -pi/2 0 pi/2 pi])
yticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'})
hold off
xlabel('x','FontName', 'Arial')
ylabel('y','FontName', 'Arial')
set(gca,'FontSize',22)

% Middle - 3D
figure
surf(grid_x,grid_y,rho(:,:,1500))
axis([-pi pi -pi pi -0.01 0.5])
zticks([0 0.25 0.5])
xlabel('x','FontName', 'Arial')
ylabel('y','FontName', 'Arial')
zlabel('Density of the swarm','FontName', 'Arial')
set(gca,'FontSize',14)

% Middle - top view
figure
surf(grid_x,grid_y,rho(:,:,1500))
view(2)
pbaspect([1 1 1])
hold on
plot3(x_fl,y_fl,ones(length(x_fl)),'r','Linewidth',3)
axis([-pi pi -pi pi -1 1.5])
xticks([-pi -pi/2 0 pi/2 pi])
xticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'})
yticks([-pi -pi/2 0 pi/2 pi])
yticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'})
hold off
xlabel('x','FontName', 'Arial')
ylabel('y','FontName', 'Arial')
set(gca,'FontSize',22)

% Third quarted - 3D
figure
surf(grid_x,grid_y,rho(:,:,2500))
axis([-pi pi -pi pi -0.01 0.5])
zticks([0 0.25 0.5])
xlabel('x','FontName', 'Arial')
ylabel('y','FontName', 'Arial')
zlabel('Density of the swarm','FontName', 'Arial')
set(gca,'FontSize',14)

% Third quarted - top view
figure
surf(grid_x,grid_y,rho(:,:,2500))
view(2)
pbaspect([1 1 1])
hold on
plot3(x_fl,y_fl,ones(length(x_fl)),'r','Linewidth',3)
axis([-pi pi -pi pi -1 1.5])
xticks([-pi -pi/2 0 pi/2 pi])
xticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'})
yticks([-pi -pi/2 0 pi/2 pi])
yticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'})
hold off
xlabel('x','FontName', 'Arial')
ylabel('y','FontName', 'Arial')
set(gca,'FontSize',22)

% Export 3D video
v = VideoWriter('S4.avi');
open(v);

figure
for i = 1:10:length(t_vec) % 10x speed
    surf(grid_x,grid_y,rho(:,:,i))
    axis([-pi pi -pi pi -0.1 0.5])
    xlabel('x','FontName', 'Arial')
    ylabel('y','FontName', 'Arial')
    zlabel('Density of the swarm','FontName', 'Arial')
    xticks([-pi -pi/2 0 pi/2 pi])
    xticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'})
    yticks([-pi -pi/2 0 pi/2 pi])
    yticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'})
    pause(0.01)
    frame = getframe(gcf);
    writeVideo(v,frame);
end
close(v);

% Export top view video
v = VideoWriter('video_top.avi');
open(v);

figure
for i = 1:10:length(t_vec) % 10x speed
    surf(grid_x,grid_y,rho(:,:,i))
    view(2)
    pbaspect([1 1 1])
    hold on
    plot3(x_fl,y_fl,ones(length(x_fl)),'r','Linewidth',3)
    hold off
    xlabel('x','FontName', 'Arial')
    ylabel('y','FontName', 'Arial')
    xticks([-pi -pi/2 0 pi/2 pi])
    xticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'})
    yticks([-pi -pi/2 0 pi/2 pi])
    yticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'})
    pause(0.01)
    frame = getframe(gcf);
    writeVideo(v,frame);
end
close(v);

%% Export
exp_initial = [reshape(grid_y,[N_xC*N_yC,1]) reshape(grid_x,[N_xC*N_yC,1]) reshape(rho(:,:,100)',[N_xC*N_yC,1])];
exp_quarter = [reshape(grid_y,[N_xC*N_yC,1]) reshape(grid_x,[N_xC*N_yC,1]) reshape(rho(:,:,500)',[N_xC*N_yC,1])];
exp_mid = [reshape(grid_y,[N_xC*N_yC,1]) reshape(grid_x,[N_xC*N_yC,1]) reshape(rho(:,:,1500)',[N_xC*N_yC,1])];
exp_third_quart = [reshape(grid_y,[N_xC*N_yC,1]) reshape(grid_x,[N_xC*N_yC,1]) reshape(rho(:,:,2500)',[N_xC*N_yC,1])];
exp_KL = [t_vec real(KL_div)];
exp_flower = [x_fl, y_fl, ones(size(x_fl))];

writematrix(exp_initial,'Swarm_initial.txt','Delimiter','tab');
writematrix(exp_quarter,'Swarm_quarter.txt','Delimiter','tab');
writematrix(exp_mid,'Swarm_mid.txt','Delimiter','tab');
writematrix(exp_third_quart,'Swarm_third_quarter.txt','Delimiter','tab');
writematrix(exp_KL,'KL_div.txt','Delimiter','tab');
writematrix(exp_flower,'flower.txt','Delimiter','tab');

save('All_data.mat')

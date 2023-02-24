clear all; close all; clc;

%% This script computes the discretized version of the kernel once for all, so that one does not need to recompute it every time the desired density is changed
%%  You need to recompute qTs if you change the parameters of the mesh

%% Parameters
N_xC = 30; % Number of cells along the x-axis
N_yC = 30; % Number of cells along the y-axis
L_x = 2*pi; % Length along the x-axis
L_y = 2*pi; % Length along the y-axis
N_sum_q = 50; % Number of summands to be used to approximate the kernel


%% Meshing
N_bar_x = N_xC+1; % Number of mesh points along the x-axis
N_bar_y = N_yC+1; % Number of mesh points along the y-axis
mesh_x = linspace(-L_x/2,L_x/2,N_bar_x)'; % Mesh along the x-axis
mesh_y = linspace(-L_y/2,L_y/2,N_bar_y)'; % Mesh along the y-axis
mesh_xC = mesh_x(1:end-1)+(mesh_x(2)-mesh_x(1))/2; % Mesh of the midpoint elements along the x-axis
mesh_yC = mesh_y(1:end-1)+(mesh_y(2)-mesh_y(1))/2; % Mesh of the midpoint elements along the y-axis

[grid_x,grid_y] = meshgrid(mesh_xC,mesh_yC); % Mesh grid of the midpoints of the cells
grid_x = grid_x';
grid_y = grid_y';

%% Compute q_T
[q_T_x,q_T_y] = compute_q_T(grid_x,grid_y,N_sum_q,L_x,L_y);

%% Export
save('qTs.mat','q_T_x','q_T_y')
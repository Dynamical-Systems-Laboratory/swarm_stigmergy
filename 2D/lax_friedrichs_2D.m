function rho_np1 = lax_friedrichs_2D(rho_n, rho_obs, dt, grid_x, grid_y,q_T_x,q_T_y)

%% Lax-Friedrichs in 2D on a rectangular grid
% Input:
% rho_n: density at time n
% rho_obs: density of the traces
% dt: Time interval in the finite volume algorithm
% grid_x: meshgrid of the x-coordinate
% grid_y: meshgrid of the y-coordinate
% q_T_x: discretized kernel on the x-coordinate
% q_T_y: discretized kernel on the y-coordinate
% Output:
% rho_np1: density at time n+1


N_xC = size(grid_x,1);
N_yC = size(grid_x,2);

dx = grid_x(2,1)-grid_x(1,1);
dy = grid_y(1,2)-grid_y(1,1);

% Compute forcing term on the x and y
[fx,fy] = compute_fx_fy(rho_n,rho_obs,grid_x,grid_y,q_T_x,q_T_y);

% Compute entire density
rho_tot = rho_n+rho_obs;

%% x-axis fluxes
% Extend for periodicity
fxm_ext = [fx(N_xC,:); fx];
fxp_ext = [fx; fx(1,:)];

rhom_ext = [rho_tot(N_xC,:); rho_tot];
rhop_ext = [rho_tot; rho_tot(1,:)];

% Fluxes for each cell in the x-direction
Fx = 0.5*(fxm_ext+fxp_ext)-dx/(4*dt)*(rhop_ext-rhom_ext);

%% y-axis fluxes
% Extend for periodicity
fym_ext = [fy(:,N_yC) fy];
fyp_ext = [fy fy(:,1)];

rhom_ext = [rho_tot(:,N_yC) rho_tot];
rhop_ext = [rho_tot rho_tot(:,1)];

% Fluxes for each cell in the y-direction
Fy = 0.5*(fym_ext+fyp_ext)-dy/(4*dt)*(rhop_ext-rhom_ext);

%% Finite volumes
Fx_L = Fx(1:N_xC,:); % Fluxes left cell
Fx_R = Fx(2:N_xC+1,:); % Fluxes right cell

Fy_B = Fy(:,1:N_yC); % Fluxes bottom cell
Fy_T = Fy(:,2:N_yC+1); % Fluxes top cell

% Lax-Friedrichs
rho_np1 = rho_n - dt/dx*(Fx_R-Fx_L) - dt/dy*(Fy_T-Fy_B);
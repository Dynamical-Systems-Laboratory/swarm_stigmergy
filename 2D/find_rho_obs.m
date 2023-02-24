function res = find_rho_obs(rho_obs,rho_des,rho_des_t,grid_x,grid_y,q_T_x,q_T_y,dt)

%% Compute the density of the traces (through fsolve)
% Input:
% rho_obs: tentative density of the traces
% rho_des: desired density of the robots
% rho_des_t: time derivative of the desired density of the robots
% grid_x: meshgrid of the x-coordinate
% grid_y: meshgrid of the y-coordinate
% q_T_x: discretized kernel on the x-coordinate
% q_T_y: discretized kernel on the y-coordinate
% dt: Time interval in the finite volume algorithm
% Output:
% res: Residual of the equations for fsolve


dx = grid_x(2,1)-grid_x(1,1);
dy = grid_y(1,2)-grid_y(1,1);

% Compute forcing term on the x and y
[fx,fy] = compute_fx_fy(rho_des,rho_obs,grid_x,grid_y,q_T_x,q_T_y);

% Compute entire density
rho_tot = rho_des+rho_obs;

%% x-axis fluxes
% Extend for periodicity
fxm_ext = [fx(end,:); fx];
fxp_ext = [fx; fx(1,:)];

rhom_ext = [rho_tot(end,:); rho_tot];
rhop_ext = [rho_tot; rho_tot(1,:)];

% Fluxes for each cell in the x-direction
Fx = 0.5*(fxm_ext+fxp_ext)-dx/(4*dt)*(rhop_ext-rhom_ext);

%% y-axis fluxes
% Extend for periodicity
fym_ext = [fy(:,end) fy];
fyp_ext = [fy fy(:,1)];

rhom_ext = [rho_tot(:,end) rho_tot];
rhop_ext = [rho_tot rho_tot(:,1)];

% Fluxes for each cell in the y-direction
Fy = 0.5*(fym_ext+fyp_ext)-dy/(4*dt)*(rhop_ext-rhom_ext);

%% Finite volumes (Lax Friedrichs)
Fx_L = Fx(1:end-1,:); % Fluxes left cell
Fx_R = Fx(2:end,:); % Fluxes right cell

Fy_B = Fy(:,1:end-1); % Fluxes bottom cell
Fy_T = Fy(:,2:end); % Fluxes top cell

res = rho_des_t*dt+dt/dx*(Fx_R-Fx_L)+dt/dy*(Fy_T-Fy_B); % Residual from finite volumes
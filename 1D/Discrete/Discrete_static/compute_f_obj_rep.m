function f = compute_f_obj_rep(rho_obj,rho_obj_t,L,mesh)

%% Compute the objective forcing term for repulsive kernel
% Input: 
% rho_obj: Objective density function
% rho_obj_t: Time derivative of the objective density function
% L: Length-scale of the repulsive kernel
% mesh: Mesh on the circle
% Output:
% f: Objective forcing term

v_rho = compute_conv_rep(rho_obj,L,mesh); % Objective velocity

f0 = -v_rho-1./rho_obj.*cumtrapz(mesh,rho_obj_t); % First part of the objective forcing term

A = -trapz(mesh,f0)/trapz(mesh,1./rho_obj); % Computation of the "constant" in the forcing term to have zero integral

f = f0+A./rho_obj; % Objective forcing term

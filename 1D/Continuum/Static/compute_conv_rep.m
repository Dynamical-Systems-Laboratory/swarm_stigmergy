function conv = compute_conv_rep(f,L,mesh)
%% Compute the convolution with the repulsive kernel
% Input: 
% f: Function to compute the convolution of
% L: Length-scale of the repulsive kernel
% mesh: Mesh on the circle
% Output:
% conv: Convolution over the circle at each point of the mesh

% Separated here in two for convenience
conv_up = exp((2*pi-mesh)/L).*cumtrapz(mesh,exp(mesh/L).*f)-exp((mesh)/L).*cumtrapz(mesh,exp(-mesh/L).*f);
conv_down = -exp((2*pi+mesh)/L).*(trapz(mesh,exp(-mesh/L).*f)-cumtrapz(mesh,exp(-mesh/L).*f))+exp((-mesh)/L).*(trapz(mesh,exp(mesh/L).*f)-cumtrapz(mesh,exp(mesh/L).*f));

conv = 1/(exp(2*pi/L)-1)*(conv_up+conv_down);
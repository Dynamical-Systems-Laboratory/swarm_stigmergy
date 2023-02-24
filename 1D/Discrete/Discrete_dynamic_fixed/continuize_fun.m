function counts = continuize_fun(f,mesh,sigma, scale_factor)
%% Continuify positions to a density
% Input:
% f: Positions to continuify
% mesh: Mesh on the circle
% sigma: Parameter of the kernel for density estimation
% scale_factor: Integral of the continuification
% Output:
% counts: Continuified density

%% Transform in degrees
deg_mesh = (mesh+pi)*180/pi;
deg_f = (f+pi)*180/pi;

%% Compute the PDF from kernel approximation
[pdf, ~] = circ_ksdensity(deg_f, deg_mesh, 'msn',sigma);
pdf(end) = pdf(1); % Make the pdf periodic

counts = pdf/trapz(mesh,pdf)*scale_factor; % Continuified density with integral different from 1
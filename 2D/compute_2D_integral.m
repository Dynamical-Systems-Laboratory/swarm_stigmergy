function integ = compute_2D_integral(f, dx, dy)

%% Computes a 2D integral with trapezoidal rule
% Input:
% f: function to integrate in 2D on an uniform grid
% dx: mesh step along x-axis
% dy: mesh step along y-axis

N_x = size(f,1);
N_y = size(f,2);

%% Integration
integ = 0;
for i = 1:N_x-1
    for j = 1:N_y-1
        integ = integ + (f(i,j)+f(i+1,j)+f(i,j+1)+f(i+1,j+1))/4*dx*dy;
    end
end
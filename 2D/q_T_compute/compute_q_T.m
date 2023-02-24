function [q_T_x,q_T_y] = compute_q_T(grid_x,grid_y,N_sum_q,L_x,L_y)
%% This function computes the discretized version of the kernel
% We generate 4D discretized kernel, which provide the kernel in the x and y directions
% between points (1,2) and (3,4) (where 1,2,3,4 indicate the dimensions of the 4D matrix)
% Input: 
% grid_x: Mesh grid (x coordinate) of the midpoints of the cells
% grid_y: Mesh grid (y coordinate) of the midpoints of the cells
% N_sum_q: Number of summands (symmetric positive/negative) to be used to approximate the kernel 
% L_x: Length along the x-axis
% L_y: Length along the y-axis
% Output:
% q_T_x: Discretized kernel (4D) in the x-direction
% q_T_y: Discretized kernel (4D) in the y-direction

N_xC = size(grid_x,1); % Number of cells along the x-axis
N_yC = size(grid_x,2); % Number of cells along the y-axis

q_T_x = zeros(N_xC,N_yC,N_xC,N_yC); % Preallocation
q_T_y = zeros(N_xC,N_yC,N_xC,N_yC); % Preallocation

%% Computation
for i = 1:N_xC
    x_i = grid_x(i,1);
    for j = 1:N_yC
        y_j = grid_y(1,j);
        for k = 1:N_xC
            x_k = grid_x(k,1);
            for l = 1:N_yC
                y_l = grid_y(1,l);
                for m = -N_sum_q:N_sum_q
                    for n = -N_sum_q:N_sum_q
                        if i==k && j==l && m==0 && n==0
                            % Avoid this term as it is singular
                        else
                            theta = atan2(y_j-y_l+L_y*n,x_i-x_k+L_x*m); % Angle between the points
                            q_T_x(i,j,k,l) = q_T_x(i,j,k,l)+exp(-sqrt((x_i-x_k+L_x*m)^2+(y_j-y_l+L_y*n)^2))*cos(theta); % Discretized kernel along the x-axis
                            q_T_y(i,j,k,l) = q_T_y(i,j,k,l)+exp(-sqrt((x_i-x_k+L_x*m)^2+(y_j-y_l+L_y*n)^2))*sin(theta); % Discretized kernel along the y-axis
                        end
                    end
                end
            end
        end
    end
    i % To track advancement
end
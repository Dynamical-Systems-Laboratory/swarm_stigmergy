function pos_obs = find_pos(counts_obs,N_obs,mesh)
%% Determine position of the obstacles based on the counts in the mesh
% Input:
% counts_obs: Counts of the obstacles in each segment of the mesh
% N_obs: total number of obstacles
% mesh: Mesh on the circle
% Output:
% pos_obs: Position of the obstacles

pos_obs = zeros(N_obs,1); % Pre-allocation
i = 1;
j = 1;
while i<N_obs
    pos_obs(i:i+counts_obs(j)-1) = mesh(j);
    i = i+counts_obs(j);
    j = j+1;
end
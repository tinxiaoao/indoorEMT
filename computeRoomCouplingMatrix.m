function [sigma_cpl, sigma_diag] = computeRoomCouplingMatrix(edgeTbl, frequency, epsilon_params, height_wall)
% Compute inter-room coupling cross-sections σ_ij (i≠j) for all frequencies.
% Inputs
%   edgeTbl       : output from findRoomInterfaces (RoomI,RoomJ,connType,Width_m,Length_m)
%   frequency     : column vector (Hz)
%   epsilon_params: material Cole-Cole parameter matrix (N×5)
%   height_wall   : wall height (m)
% Outputs
%   sigma_cpl     : N×N×F array (i≠j entries = σ_ij, i=j = 0)
%   sigma_diag    : N×F   matrix, each row i = Σ_{k≠i} σ_ik

rooms = unique([edgeTbl.RoomI; edgeTbl.RoomJ]);
N = max(rooms);
F = numel(frequency);

sigma_cpl  = zeros(N, N, F);
sigma_diag = zeros(N,     F);

for idx = 1:height(edgeTbl)
    i = edgeTbl.RoomI(idx);
    j = edgeTbl.RoomJ(idx);
    conn = char(edgeTbl.connType(idx));
    width_m = edgeTbl.Width_m(idx);

    if strcmp(conn,'door')
        area = (edgeTbl.Length_m(idx) * height_wall) / 2; % 门的面积只考虑几何面积一半
        SUM_T = ones(1,F); % 门视为开口，无材料损耗
    elseif strcmp(conn,'window')
        area = edgeTbl.Length_m(idx) * height_wall;
        layer_vec = [width_m/2, 0.012, width_m/2];
        res = multi_layer_model(frequency, 4, layer_vec, epsilon_params);
        SUM_T = res.SUM_Solid_Angle_mean_T.'; 
    else % wall
        area = edgeTbl.Length_m(idx) * height_wall;
        layer_vec = [0.005, width_m, 0.005];
        res = multi_layer_model(frequency, 2, layer_vec, epsilon_params);
        SUM_T = res.SUM_Solid_Angle_mean_T.';
    end

    sigma_vec = 0.5 * area .* SUM_T;    % ½·A·SUM_T (1×F)

    sigma_cpl(i,j,:) = sigma_cpl(i,j,:) + reshape(sigma_vec,1,1,F);
    sigma_cpl(j,i,:) = sigma_cpl(j,i,:) + reshape(sigma_vec,1,1,F); % assume reciprocity

    sigma_diag(i,:)  = sigma_diag(i,:)  + sigma_vec;
    sigma_diag(j,:)  = sigma_diag(j,:)  + sigma_vec;
end

end

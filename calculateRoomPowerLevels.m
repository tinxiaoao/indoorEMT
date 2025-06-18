function E_V_per_m = calculateRoomPowerLevels(M, P_in_dBm)
% calculateRoomPowerLevels - 计算房间稳态平均电场强度 (V/m)
%
% 输入：
%   M          - 耦合截面传播矩阵 (NxNxF), 单位 m², N为房间数量，F为频点数量
%   P_in_dBm   - 每个房间的输入功率 (Nx1)，单位 dBm
%
% 输出：
%   E_V_per_m  - 各房间电场强度 (NxF)，单位 V/m
%
% 示例：
%   M = rand(2, 2, 51);                  % 示例耦合截面传播矩阵 (m²)
%   P_in_dBm = [0; -Inf];                % 房间1输入0 dBm，房间2无输入
%   E_V_per_m = calculateRoomPowerLevels(M, P_in_dBm);

Z0 = 377;  % 自由空间波阻抗 (Ω)

% 转换输入功率 dBm -> W
P_in_W = 10.^((P_in_dBm - 30)/10);
P_in_W(isinf(P_in_W)) = 0;

F = size(M, 3);
N = size(M, 1);

S_out_W_m2 = zeros(N, F);
E_V_per_m = zeros(N, F);

for idx = 1:F
    M_f = M(:, :, idx);

    % 稳定求解房间功率密度
    if rcond(M_f) < 1e-12
        warning('频点 %d 传播矩阵条件数异常，使用伪逆求解', idx);
        S_out_W_m2(:, idx) = pinv(M_f) * P_in_W;
    else
        S_out_W_m2(:, idx) = M_f \ P_in_W;
    end

    % 避免负值或零值异常
    S_out_W_m2(S_out_W_m2(:, idx) <= 1e-12, idx) = 1e-12;

    % 场强计算 (V/m)
    E_V_per_m(:, idx) = sqrt(S_out_W_m2(:, idx) * Z0);
end

end

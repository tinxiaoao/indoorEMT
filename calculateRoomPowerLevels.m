function S_out_dBm = calculateRoomPowerLevels(M, S_in_dBm)
% calculateRoomPowerLevels - 计算房间稳态平均功率水平
%
% 输入：
%   M        - 传播矩阵 (NxN)
%   S_in_dBm - 源输入功率向量 (Nx1)，单位：dBm
%
% 输出：
%   S_out_dBm - 各房间的稳态平均功率水平 (Nx1)，单位：dBm
%
% 示例：
%   M = [1.2, -0.2; -0.1, 1.1];  % 示例传播矩阵
%   S_in_dBm = [0; -Inf];        % 房间1输入功率0 dBm，房间2无源输入(-Inf)
%   S_out_dBm = calculateRoomPowerLevels(M, S_in_dBm);

% 将输入功率从dBm转换为线性功率(mW)
S_in_mW = 10.^(S_in_dBm / 10);

% 处理输入-Inf (对应功率0 mW)
S_in_mW(isinf(S_in_mW)) = 0;

% 求解稳态方程：M * S_out = S_in (mW)
S_out_mW = M \ S_in_mW;

% 将线性功率(mW)转换回dBm
S_out_dBm = 10 * log10(S_out_mW);

end

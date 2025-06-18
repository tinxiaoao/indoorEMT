function csvFile = saveRoomPowerLong(roomIDs, S_out_dBm, frequency, outDir, baseID)
% saveRoomPowerLong  把 (房间 × 频率) 功率矩阵展平成长表并写 CSV
%
% csvFile = saveRoomPowerLong(roomIDs, S_out_dBm, frequency, outDir, baseID)
%
% 输入:
%   roomIDs     : 1×R string / char 向量或 R×1 cellstr
%   S_out_dBm   : R×N double  —— 功率矩阵 (dBm)
%   frequency   : N×1 double  —— 频率向量 (Hz)
%   outDir      : CSV 保存目录 (string / char)，默认 'E:\code\IndoorEMT'
%   baseID      : 楼层/平面图编号，用于生成文件名，可为空 ''
%
% 输出:
%   csvFile     : 实际写入的 CSV 完整路径 (char)

% -------------------------------------------------------------------------
% 参数检查与默认值
% -------------------------------------------------------------------------
if nargin < 4 || isempty(outDir)
    outDir = "E:\code\IndoorEMT";
end
if nargin < 5
    baseID = "";
end

% 确保目录存在
if ~isfolder(outDir)
    mkdir(outDir);
end

% 统一 RoomID 为 string 列向量
roomIDs = string(roomIDs(:));           % R×1
[R, N]  = size(S_out_dBm);              % R=房间数, N=频点

% 频率转 GHz 列向量
freqGHz = frequency(:) / 1e9;           % N×1

% -------------------------------------------------------------------------
% 组装长表 (vectorized, no loop)
% -------------------------------------------------------------------------
RoomCol      = repmat(roomIDs , N, 1);              % (R*N)×1
FreqCol_GHz  = repmat(freqGHz , R, 1);              % (R*N)×1
PowerCol_dBm = reshape(S_out_dBm.', R*N, 1);        % 列优先展开, (R*N)×1

resultTable_long = table(RoomCol, FreqCol_GHz, PowerCol_dBm, ...
                         'VariableNames', ...
                         {'RoomID','Frequency_GHz','SteadyStatePower_dBm'});

% -------------------------------------------------------------------------
% 写 CSV (UTF-8)
% -------------------------------------------------------------------------
timestamp = datestr(now, 'yyyymmdd_HHMMSS');
if strlength(baseID) > 0
    fname = sprintf('roomPower_%s_%s.csv', baseID, timestamp);
else
    fname = sprintf('roomPower_%s.csv', timestamp);
end
csvFile = fullfile(outDir, fname);
writetable(resultTable_long, csvFile, 'Encoding','UTF-8');

fprintf('Room-power CSV written: %s\n', csvFile);
end

%% batch_sigma_loss.m
% -------------------------------------------------------------
% 批量计算 σ_loss (房间功率损耗截面积, 单位 m²)
% 遍历 roomself_excel/<ID>.xlsx 与 topology_excel/<ID>.xlsx 配对
% —— 仅打印结果，不写文件 ——
% 修复: 若 RoomID 列非 cellstr，自动转换为 string 向量，
%       避免 "此类型的变量不支持使用花括号进行索引" 错误。
% -------------------------------------------------------------
% 中文列头：
%   房间表:  房间ID | 类型 | 面积 | 周长 | 外墙 length | 外墙 width
%   拓扑表: 房间1   | 房间2 | 连接类型 | 数量 | length | width
% -------------------------------------------------------------

close all; clear; clc;

%% 0. 配置 ---------------------------------------------------------------
rootRoom  = 'E:/code/CubiCasa5k/output/roomself_excel';
rootTopo  = 'E:/code/CubiCasa5k/output/topology_excel';
h_real    = 2.2;   % 层高 (m)
door_real = 0.9;   % 标准门宽 (m)

roomMap = {'外墙 length','outerWall_len_px';
           '外墙 width','outerWall_wid_px';
           '房间ID','RoomID'};

topoMap = {'连接类型','Type';
           '房间1','RoomA';
           '房间2','RoomB';
           'length','len_px';
           'width','wid_px'};

%% 1. 遍历房间表 ---------------------------------------------------------
roomFiles = dir(fullfile(rootRoom,'*.xlsx'));

for f = 1:numel(roomFiles)
    fnameRoom = fullfile(rootRoom, roomFiles(f).name);
    fnameBase = erase(roomFiles(f).name,'.xlsx');
    fnameTopo = fullfile(rootTopo, [fnameBase '.xlsx']);

    fprintf('\n>> Processing "%s"\n', fnameBase);

    try
        %% 2. 读取房间表 --------------------------------------------------
        optsRoom = detectImportOptions(fnameRoom,'VariableNamingRule','preserve');
        for k = 1:size(roomMap,1)
            idx = strcmp(optsRoom.VariableNames, roomMap{k,1});
            if any(idx), optsRoom.VariableNames{idx} = roomMap{k,2}; end
        end
        tblRoom = readtable(fnameRoom, optsRoom);

        %% 3. 读取拓扑表 --------------------------------------------------
        optsTopo = detectImportOptions(fnameTopo,'VariableNamingRule','preserve');
        for k = 1:size(topoMap,1)
            idx = strcmp(optsTopo.VariableNames, topoMap{k,1});
            if any(idx), optsTopo.VariableNames{idx} = topoMap{k,2}; end
        end
        tblTopo = readtable(fnameTopo, optsTopo);

        %% 4. 估算像素→米比例 r (门宽法，若找不到门则 r=0.01) -------------
        doorFlag = strcmp(tblTopo.Type,'door');
        if any(doorFlag)
            door_w_px = tblTopo.wid_px(doorFlag);
            r = door_real / mean(door_w_px,'omitnan');
        else
            r = 0.01;  % 默认比例
        end
        fprintf('   r = %.4f m/px\n', r);

        %% 5. 统一 RoomID/RoomA/RoomB 为字符串 --------------------------
        tblRoom.RoomID = string(tblRoom.RoomID);
        tblTopo.RoomA  = string(tblTopo.RoomA);
        tblTopo.RoomB  = string(tblTopo.RoomB);

        %% 6. 计算 σ_loss 并打印 -----------------------------------------
        for i = 1:height(tblRoom)
            roomID = tblRoom.RoomID(i);
            Aw_px  = tblRoom.outerWall_len_px(i) * tblRoom.outerWall_wid_px(i);
            maskA  = tblTopo.RoomA==roomID & tblTopo.Type=="wall";
            maskB  = tblTopo.RoomB==roomID & tblTopo.Type=="wall";
            Ai_px  = sum(tblTopo.len_px(maskA|maskB) .* tblTopo.wid_px(maskA|maskB));
            sigma_loss = (Aw_px + Ai_px) * (h_real/r) * r^2; % m²
            fprintf('      Room %-10s  σ_loss = %.3f m²\n', roomID, sigma_loss);
        end

    catch ME
        fprintf(2,'   ERROR: %s\n', ME.message);
    end
end

fprintf('\nAll files processed.\n');

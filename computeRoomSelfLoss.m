function sigma_loss_tbl = computeRoomSelfLoss(wideTbl, frequency, px2m, height_wall, epsilon_params)
    % 1. 初始化
    roomIDs = wideTbl.RoomID;
    N = height(wideTbl);
    numFreq = length(frequency);
    sigma_data = zeros(N, numFreq);
    
    % 2. 自动识别内墙段列
    varNames = wideTbl.Properties.VariableNames;
    segLenMask = ~cellfun('isempty', regexp(varNames, '^Seg\d+Len_px$'));
    segWidMask = ~cellfun('isempty', regexp(varNames, '^Seg\d+Wid_px$'));
    segLenCols = varNames(segLenMask);
    segWidCols = varNames(segWidMask);
    
    % 3. 遍历每个房间计算损耗截面积
    for i = 1:N
        % 提取该房间行的数据
        row = wideTbl(i,:);
        % 计算 S_ext, S_int, S_cf
        S_ext = row.ExtWallLength_px * height_wall * px2m;
        % 内墙长度求和（跳过 NaN）
        segLens = row{:, segLenCols};
        total_int_length = nansum(segLens);
        S_int = total_int_length * height_wall * px2m;
        % 顶棚和地板面积
        S_cf = row.Area_px * 2 * (px2m^2);
        
        % 计算平均内墙厚度（若无隔墙则给0）
        segWids = row{:, segWidCols};
        if all(isnan(segWids))
            avg_int_width = 0;
        else
            avg_int_width = nanmean(segWids) * px2m;
        end
        
        % 4. 获取各构件损耗系数 R1, R2, R3
        R1 = multi_layer_model(frequency, epsilon_params, [15 19 15], row.ExtWallWidth_px * px2m);
        R2 = multi_layer_model(frequency, epsilon_params, [15 16 15], avg_int_width);
        R3 = multi_layer_model(frequency, epsilon_params, [15 19 15], 0.1);  % 0.1 m 固定厚度
        
        % 5. 计算各频率下损耗截面积 sigma_vec
        sigma_vec = 0.5 * (S_ext * R1 + S_int * R2 + S_cf * R3);
        sigma_data(i, :) = sigma_vec;
    end
    
    % 6. 构造输出表格
    freq_GHz = frequency / 1e9;
    colNames = arrayfun(@(x) sprintf('SigmaLoss_f%.1fGHz', x), freq_GHz, 'UniformOutput', false);
    colNames = strrep(colNames, '.', '_');
    sigma_loss_tbl = array2table(sigma_data, 'VariableNames', colNames);
    sigma_loss_tbl.RoomID = roomIDs;
    sigma_loss_tbl = movevars(sigma_loss_tbl, 'RoomID', 'Before', 1);
end

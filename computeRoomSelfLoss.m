function [sigma_loss, roomIDs] = computeRoomSelfLoss(wideTbl, frequency, px2m, height_wall, epsilon_params)
% Compute room self-loss (sigma_loss) based on multi_layer_model

roomIDs = wideTbl.RoomID;
N = height(wideTbl);
numFreq = length(frequency);
sigma_loss = zeros(N, numFreq);

% Identify inner wall segment columns
varNames = wideTbl.Properties.VariableNames;
segLenCols = regexp(varNames, '^Seg\d+Len_px$', 'match');
segWidCols = regexp(varNames, '^Seg\d+Wid_px$', 'match');
segLenCols = [segLenCols{:}];
segWidCols = [segWidCols{:}];

for i = 1:N
    % Extract room data
    row = wideTbl(i,:);

    % External wall area and thickness
    S_ext = row.ExtWallLength_px * height_wall * px2m;
    ext_thickness = row.ExtWallWidth_px * px2m * 0.1;  % 0.1是系数修正，因为有些平面图外墙故意画的很粗，导致远超实际

    % Internal wall area (per segment)
    sigma_int = zeros(1, numFreq);
    for k = 1:numel(segLenCols)
        if isnan(row.(segLenCols{k})) || isnan(row.(segWidCols{k}))
            continue;
        end
        seg_length = row.(segLenCols{k}) * px2m;
        seg_thickness = row.(segWidCols{k}) * px2m;
        seg_area = seg_length * height_wall;
        res_int = multi_layer_model(frequency, 2, seg_thickness, epsilon_params);
        sigma_int = sigma_int + seg_area .* res_int.SUM_Solid_Angle_mean_R' / 2;
    end

    % Ceiling and floor area
    S_cf = row.Area_px * 2 * px2m^2;

    % Call multi_layer_model for each component
    res_ext = multi_layer_model(frequency, 1, ext_thickness, epsilon_params);
    res_cf  = multi_layer_model(frequency, 3, 0.1, epsilon_params);
    sigma_ext = S_ext * res_ext.SUM_Solid_Angle_mean_R' / 2;
    sigma_cf  = S_cf  * res_cf.SUM_Solid_Angle_mean_R' / 2;

    % Sum sigma_loss for each frequency
    sigma_loss(i, :) = sigma_ext + sigma_int + sigma_cf;
end

end

function edgeTbl = roomTransArea(topoTbl)
%ROOMTRANSAREA  Collect wall / door / window connections and compute areas.
%   edgeTbl = ROOMTRANSAREA(topoTbl)  (topology table with Chinese headers)
%
%   Return columns:
%       RoomA   RoomB   ConnType    len_px  wid_px  area_px
%   ConnType = 'wall' | 'door' | 'window'
%   All numeric outputs are pixel units.
%
% -----------------------------------------------------------------------
% Robust against cell array column names (uses strcmp).
% -----------------------------------------------------------------------

% ---------- Column name mapping ----------------------------------------
map = {
    '连接类型','ConnType';
    '房间1','RoomA';
    '房间2','RoomB';
    'length','len_px';
    'width','wid_px'};

for k = 1:size(map,1)
    colCN = map{k,1};
    colEN = map{k,2};
    idx   = strcmp(topoTbl.Properties.VariableNames, colCN);
    if any(idx)
        topoTbl.Properties.VariableNames{idx} = colEN;
    end
end

% ---------- keep only wall / door / window rows ------------------------
keepRows = ismember(topoTbl.ConnType, {'wall','door','window','墙','门','窗'});
subTbl   = topoTbl(keepRows, {'RoomA','RoomB','ConnType','len_px','wid_px'});

% convert Chinese conn type to English uniform --------------------------
subTbl.ConnType = replace(subTbl.ConnType, {'墙','门','窗'}, {'wall','door','window'});

% ---------- compute pixel area ----------------------------------------
subTbl.area_px = subTbl.len_px .* subTbl.wid_px;

% ---------- cast IDs to string ----------------------------------------
subTbl.RoomA   = string(subTbl.RoomA);
subTbl.RoomB   = string(subTbl.RoomB);

edgeTbl = subTbl;
end

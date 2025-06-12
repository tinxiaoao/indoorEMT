function edgeTbl = findRoomInterfaces(topoTbl, px2m, height_wall)
% FINDROOMINTERFACES Extract wall/door/window interfaces between room pairs
% from a topology table.
%
% Inputs
% -------
% topoTbl : table imported from topology_excel. Must contain columns
% RoomA, RoomB, ConnectionType, Length, Width (pixel units).
% px2m : scalar. Pixel-to-metre scale factor.
% height_wall : wall height in metres.
%
% Output
% ------
% edgeTbl : table, each row = one interface segment.
% Variables:
% RoomI, RoomJ uint16 / double (ascending order)
% connType categorical {'wall','door','window'}
% Length_m, Width_m double (converted to m)
% Area_m2 double (full geometric area)
%
% Notes
% -----
% - open/opening types are ignored.
% - Duplicate segments with the same rooms, type, and thickness are merged.
% - For doors, the geometric area is kept complete. A factor of 1/2 is applied later
% during coupling cross-section calculation due to σ_ij = Area×T/2 with T=1.

arguments
    topoTbl table
    px2m    (1,1) double {mustBePositive}
    height_wall (1,1) double {mustBePositive}
end

%% Strict column names (case-sensitive)
requiredCols = {'RoomA','RoomB','ConnectionType','Length','Width'};
missing = setdiff(requiredCols, topoTbl.Properties.VariableNames);
if ~isempty(missing)
    error('Topology table lacks required columns: %s', strjoin(missing,', '));
end

cA = 'RoomA';
cB = 'RoomB';
cT = 'ConnectionType';
cL = 'Length';
cW = 'Width';

%% Filter rows: only wall/door/window
validMask = ismember(lower(string(topoTbl.(cT))), {'wall','door','window'});
T = topoTbl(validMask, :);

%% Build edge rows
nSeg = height(T);
RoomI   = zeros(nSeg,1,'uint16');
RoomJ   = zeros(nSeg,1,'uint16');
connStr = strings(nSeg,1);
Len_m   = zeros(nSeg,1);
Wid_m   = zeros(nSeg,1);
Area_m2 = zeros(nSeg,1);

for r = 1:nSeg
    a = T.(cA)(r);  b = T.(cB)(r);
    if a == b, continue; end % Skip self-loops
    [RoomI(r), RoomJ(r)] = deal(min(a,b), max(a,b));
    connStr(r) = lower(string(T.(cT)(r)));
    Len_m(r)   = T.(cL)(r) * px2m;
    Wid_m(r)   = T.(cW)(r) * px2m;
    Area_m2(r) = Len_m(r) * height_wall; % Complete geometric area for doors
end

%% Assemble output table, drop zero rows (self loops)
keep = RoomI ~= 0;
edgeTbl = table(RoomI(keep), RoomJ(keep), categorical(connStr(keep)), ...
                Len_m(keep), Wid_m(keep), Area_m2(keep), ...
                'VariableNames', {'RoomI','RoomJ','connType','Length_m','Width_m','Area_m2'});

%% Merge duplicate segments with same rooms, type, thickness
[gID,~,gidx] = findgroups(edgeTbl.RoomI, edgeTbl.RoomJ, edgeTbl.connType, edgeTbl.Width_m);
if numel(gID) < height(edgeTbl)
    aggLen  = splitapply(@sum, edgeTbl.Length_m, gidx);
    aggArea = splitapply(@sum, edgeTbl.Area_m2, gidx);
    edgeTbl = unique(edgeTbl(:,{'RoomI','RoomJ','connType','Width_m'}));
    edgeTbl.Length_m = aggLen;
    edgeTbl.Area_m2  = aggArea;
end

end

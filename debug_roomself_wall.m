% Debug and test the overall electromagnetic propagation system

close all; clear; clc;

%% 参数设置
BASE_ID  = "1";
rootRoom = "E:/code/CubiCasa5k/output/roomself_excel";
rootTopo = "E:/code/CubiCasa5k/output/topology_excel";

roomFile = fullfile(rootRoom, BASE_ID + ".xlsx");
topoFile = fullfile(rootTopo, BASE_ID + ".xlsx");

%% 读取数据
roomTbl = readtable(roomFile, 'VariableNamingRule', 'preserve');
topoTbl = readtable(topoFile, 'VariableNamingRule', 'preserve');

epsilon_params = load('epsilon_paper_tabel.txt');

%% 转换系数 (示例: 已知实际面积和对应的像素面积)
real_area = 55.5;   % 示例：20平方米
pixel_area = 364620; % 示例：对应4000像素面积
px2m = sqrt(real_area / pixel_area);

%% 墙高设置
height_wall = 2.2; % 墙高2.2米

%% 频率设置
frequency = (1e9:1e8:6e9)'; % 1GHz 到 6GHz，步进100MHz

%% 计算房间墙段信息
wideTbl = roomWallSegments(roomTbl, topoTbl);

%% 计算房间自身损耗
[sigma_loss, roomIDs] = computeRoomSelfLoss(wideTbl, frequency, px2m, height_wall, epsilon_params);

%% 计算房间接口（墙、门、窗）
edgeTbl = findRoomInterfaces(topoTbl, px2m, height_wall);

%% 计算房间耦合矩阵
[sigma_cpl, sigma_diag] = computeRoomCouplingMatrix(edgeTbl, frequency, epsilon_params, height_wall);

%% 构建完整传播矩阵
M = assemblePropagationMatrix(sigma_loss, sigma_diag, sigma_cpl);

%% 定义输入源矩阵（示例：第一个房间0dBm，其他房间无源-Inf dBm）
S_in_dBm = -Inf(size(roomIDs));
S_in_dBm(1) = 0;

%% 计算各房间稳态功率水平
S_out_dBm = calculateRoomPowerLevels(M, S_in_dBm);

%% 显示结果
resultTable = table(roomIDs, S_out_dBm, 'VariableNames', {'RoomID', 'SteadyStatePower_dBm'});
disp(resultTable);

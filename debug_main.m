%% run_roomwall_demo.m  —  Example usage for a single floor plan
%% debug_main.m  —  run full workflow on a single floor‑plan (ID = 1)
% -----------------------------------------------------------------------
% 修改 BASE_ID 以切换到其他编号的平面图；脚本将：
%   1) 读取房间信息表
%   2) 读取拓扑信息表
%   3) 调用 allroomwall 计算每个房间的墙体统计
%   4) 调用 roomTransArea 生成房间对的连接信息
%   5) 显示结果表
% -----------------------------------------------------------------------

close all; clear; clc;

%% 0. 参数：文件编号 ------------------------------------------------------
BASE_ID = "1";   % 改成 '2'、'3' 等可切换测试编号
rootRoom = "E:/code/CubiCasa5k/output/roomself_excel";
rootTopo = "E:/code/CubiCasa5k/output/topology_excel";

roomFile = fullfile(rootRoom, BASE_ID + ".xlsx");
topoFile = fullfile(rootTopo, BASE_ID + ".xlsx");

fprintf("Reading files:\n  roomFile = %s\n  topoFile = %s\n\n", roomFile, topoFile);

%% 1. 读取房间表 ---------------------------------------------------------
optsR    = detectImportOptions(roomFile,'VariableNamingRule','preserve');
roomTbl  = readtable(roomFile, optsR);

%% 2. 读取拓扑表 ---------------------------------------------------------
optsT    = detectImportOptions(topoFile,'VariableNamingRule','preserve');
topoTbl  = readtable(topoFile, optsT);

%% 3. 各房间墙体统计-------------------------------------------
resultTbl = roomWallSegments(roomTbl, topoTbl);

fprintf("\n=== Room Wall Statistics (pixels) ===\n");
disp(resultTbl);

% %% 4. 房间对连接信息 (墙/门/窗) -----------------------------------------
% edgeTbl = roomTransArea(topoTbl);
% 
% fprintf("\n=== Room‑to‑Room Connections (pixels) ===\n");
% disp(edgeTbl);


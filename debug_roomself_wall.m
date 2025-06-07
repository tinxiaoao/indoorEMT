close all; clear; clc;

BASE_ID  = "1";
rootRoom = "E:/code/CubiCasa5k/output/roomself_excel";
rootTopo = "E:/code/CubiCasa5k/output/topology_excel";

roomFile = fullfile(rootRoom, BASE_ID + ".xlsx");
topoFile = fullfile(rootTopo, BASE_ID + ".xlsx");

optsR   = detectImportOptions(roomFile, 'VariableNamingRule','preserve');
roomTbl = readtable(roomFile, optsR);

optsT   = detectImportOptions(topoFile, 'VariableNamingRule','preserve');
topoTbl = readtable(topoFile, optsT);

wideTbl = roomWallSegments(roomTbl, topoTbl);
disp(wideTbl);

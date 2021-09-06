close all; clearvars; clc;
% A test script for isStructureComplete

%%
struct_a.a = 15;
struct_a.b = 2;

if isStructureComplete(struct_a)
    disp('Pass')
else
    disp('Fail')
end

%%
struct_a.struct_b = struct_a;

if isStructureComplete(struct_a)
    disp('Pass')
else
    disp('Fail')
end

%%
struct_a.struct_b.c = [];

if ~isStructureComplete(struct_a)
    disp('Pass')
else
    disp('Fail')
end

%%
struct_a.struct_b.c = 56;
struct_a.c = [];

if ~isStructureComplete(struct_a)
    disp('Pass')
else
    disp('Fail')
end
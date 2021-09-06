clearvars; close all; clc;
% Test script for testing faster algorithms in addValueToStruct and
% isCompleteStructure functions

struct_a = struct('a',[],'b',[],'c',[]);
struct_b = struct('b',[],'c',[],'d',[]);
struct_c = struct('c',[],'d',[],'e',[]);
struct_d = struct('a',struct_a,'b',struct_b,'c',struct_c);


struct_d = addValuesToStruct(struct_d,{'a','b','c','d','e'},{1,2,3,4,5});

struct_d = addValuesToStruct(struct_d,struct('a',5,'b',4,'d',2,'e',1));

struct_d = addValuesToStruct(struct_d,'a',2,'b',1,'d',5,'e',4);
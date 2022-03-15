%% mergeStruct Test

StructA = struct('A',struct('a',1,'b',2,'c',3),'B',struct('c',5,'d',6));
StructB = struct('A',struct('a',55,'c',[]),'B',struct('b',4,'c',5,'d',6),'C',struct('a',3,'e',5,'n',[]));

mergeStructs(StructA,StructB)
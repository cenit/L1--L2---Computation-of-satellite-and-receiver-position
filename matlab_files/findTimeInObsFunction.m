function [ n,satamount ] = findTimeInObsFunction( lovFile,mine )
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here
lov033b = lovFile;
IndexC = strfind(lov033b, 'END');
Index = find(not(cellfun('isempty', IndexC)));
n = Index + 1;
%%
satamount = sscanf(lov033b{n,8}, '%f',[1 Inf]);
sam = str2double(lov033b{n,2});
samp=[str2double(lov033b{n,3}),str2double(lov033b{n,4}),str2double(lov033b{n,5}),str2double(lov033b{n,6})];
myTiming = isequal(mine,samp);
%%
for i = 1:100
    n = n + satamount*2+1;
    satamount = sscanf(lov033b{n,8}, '%f',[1 Inf]);
    sam = str2double(lov033b{n,2});
    samp=[str2double(lov033b(n,3)),str2double(lov033b{n,4}),str2double(lov033b{n,5}),str2double(lov033b{n,6})];
    myTiming = isequal(mine,samp);
    if myTiming
        return
    end
end
[strsplit(lov033b{n,8},'G'),strsplit(lov033b{20,9},'G'),strsplit(lov033b{20,10},'G'),strsplit(lov033b{20,10},'G')]
end


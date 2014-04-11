function [ xvalue,yvalue,zvalue ] = samplefunction( lov033b )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
IndexC = strfind(lov033b, 'APPROX');
Index = find(not(cellfun('isempty', IndexC)));
%%
xvalue = str2double(lov033b{Index-length(lov033b)});
yvalue = str2double(lov033b{Index-2*length(lov033b)});
zvalue = str2double(lov033b{Index-3*length(lov033b)});

end


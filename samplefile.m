lov033b = importObserverFileAsString('0lov033b.04o', 1, 5629);
IndexC = strfind(lov033b, 'APPROX');
Index = find(not(cellfun('isempty', IndexC)));
%%
xvalue = str2double(lov033b{Index-length(lov033b)});
yvalue = str2double(lov033b{Index-2*length(lov033b)});
zvalue = str2double(lov033b{Index-3*length(lov033b)});
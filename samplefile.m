clear all
lov033b = importObserverFileAsString('0lov033b.04o', 1, 5629);
IndexC = strfind(lov033b, 'END');
Index = find(not(cellfun('isempty', IndexC)));
lov033b{Index+1,2}
%%
xyc1 = sscanf(lov033b{Index+1,8}, '%f',[1 Inf])
mine=[2,1,14,0];
sam = str2double(lov033b{Index+1,2});
[str2double(lov033b{Index+1,6})]
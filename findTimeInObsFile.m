
lov033b = importObserverFileAsString('0lov033b.04o', 1, 5629);
IndexC = strfind(lov033b, 'END');
Index = find(not(cellfun('isempty', IndexC)));
n = Index + 1;
lov033b{n,2}
mine=[2,1,14,0];
%%
satamount = sscanf(lov033b{n,8}, '%f',[1 Inf])
sam = str2double(lov033b{n,2});
samp=[str2double(lov033b{n,3}),str2double(lov033b{n,4}),str2double(lov033b{n,5}),str2double(lov033b{n,6})]
myTiming = isequal(mine,samp);
%%
for i = 1:100
    n = n + satamount*2+1;
    satamount = sscanf(lov033b{n,8}, '%f',[1 Inf]);
    sam = str2double(lov033b{n,2});
    samp=[str2double(lov033b(n,3)),str2double(lov033b{n,4}),str2double(lov033b{n,5}),str2double(lov033b{n,6})];
    myTiming = isequal(mine,samp);
    if myTiming
        disp('yay');
        disp(sprintf('%g satellites,(%g days,%g hour,%g minutes,%g seconds)',satamount,samp(1),samp(2),samp(3),samp(4)));
        return
    end
end

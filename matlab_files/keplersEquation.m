function Ek = keplersEquation( mk,ec )
E0(1) = mk;
n=7;
for i = 1:n
    E0(i+1) = mk + ec*sin(E0(i));
    eps = abs(E0(i)-E0(i+1));
end
Ek = E0(end);
end


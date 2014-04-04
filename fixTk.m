function tk = fixTk( tk )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
if tk > 302400
    tk = tk - 604800;
elseif tk < -302400
    tk = tk + 604800;
else
    tk = tk;
end


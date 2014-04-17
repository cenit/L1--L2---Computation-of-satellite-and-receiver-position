% Function used to fix tk value
function tk = fixTk( tk )
if tk > 302400
    tk = tk - 604800;
elseif tk < -302400
    tk = tk + 604800;
else
    tk = tk;
end


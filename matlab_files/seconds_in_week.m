% Function that converts days,hours,min,seconds to total seconds
function total_seconds = seconds_in_week( days, hours, min, sec )
total_seconds = (days-1)*24*3600 + hours*3600 + min*60 + sec;
end


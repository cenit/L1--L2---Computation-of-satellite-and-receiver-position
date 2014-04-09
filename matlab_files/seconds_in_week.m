function total_seconds = seconds_in_week( days, hours, min, sec )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
total_seconds = (days-1)*24*3600 + hours*3600 + min*60 + sec;
end


function [ApproxPos, ObsTypes] = ReadObsHeader(fid)
%ReadObsHeader   read selected header values from RINEX observation file
%
%   Usage:
%      [ApproxPos, ObsTypes] = ReadObsHeader(fid)
%   Input:
%      fid is file identifier for input file
%   Output:
%      ApproxPos is approximate Cartesian position vector of 
%      the receiver in meters
%      ObsTypes is array of observation type strings, i.e. C1, P1, D2
%
%   See also ReadObsRecord

while ~feof(fid)
	line = fgets(fid);
	tag = line(61:end);
	if strncmpi(tag, 'APPROX POSITION XYZ', 19)
		ApproxPos = parsef(line, {3,'F14.4'});
	elseif strncmpi(tag, '# / TYPES OF OBSERV', 19)
		[num, ObsTypes] = parsef(line, {'I6',9,{4,'X','A2'}});
		ObsTypes = ObsTypes(1:num);
	elseif (strncmpi(tag, 'END OF HEADER', 13))
		break
	end
end
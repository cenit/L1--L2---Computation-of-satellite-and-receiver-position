function [Epoch, Data] = ReadObsRecord(fid, numobs, PRN)
%ReadObsRecord   read one epoch of data values from RINEX observation file
%
%   Usage:
%      [Epoch, Data] = ReadObsRecord(fid, numobs, [PRN])
%   Input:
%      fid is file identifier for input file (determined by using fopen)
%      numobs is number of observation types, i.e. if the observations
%      types are L1 L2 C1, then  numobs = 3
%      PRN is (optional) list of selected satellites.  Only use it 
%        if you only want one satellite.  Otherwise, all are returned.
%   Output:
%      Epoch is struct containing year,mon,day,hour,min,sec of epoch
%      this means to access this information, you use
%         Epoch.year
%         Epoch.mon
%         Epoch.day etc
%   
%      Data is struct array of PRN, Val, LLI, SNR for each satellite
%         if you want the C1 data for PRN 12, you need to know
%         where the C1 data are stored and where the PRN 12 data are stored.
%
%         If the observable types are stored in ObsTypes (read in ReadHeader)     
%         C1index = strmatch('C1', ObsTypes) gives that index
%
%         The PRN number will be stored in Data.PRN
%         If you are looking for PRN 12
%         satindex = find([Data.PRN] == 12);
%         
%         Thus, Data(satindex).Val(C1index) gives you C1 data for PRN 12
%
%   See also ReadObsHeader

if nargin < 3
	PRN = [];
end

% parse data record header line from file
[date,secs,flag,numsats,sats] = parsef(fgets(fid), ...
		{5,{'X','I2'},'F11.7',2,'X','I1','I3',12,{'X','I2'}});

% put epoch vals into struct
Epoch = struct('year',date(1),'mon',date(2),'day',date(3),'hour',date(4),'min',date(5),'sec',secs);

% create empty Data struct
Data = struct('PRN',{},'Val',{},'LLI',{},'SNR',{});

% read the data record items, one for each satellite
for satnum = 1:numsats
	record = '';
	for linenum = 0:(numobs-1)/5
		[line, term] = fgets(fid);
		record = [record, line(1:end-length(term)), blanks(80-length(line)+length(term))];
	end
	if isempty(PRN) | find(PRN == sats(satnum))
		data = parsef(record, {numobs,{'F14.3','I1','I1'}});
		Data(end+1) = struct('PRN',sats(satnum),'Val',data{1},'LLI',data{2},'SNR',data{3});
	end
end
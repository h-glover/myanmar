function adcp = adcp_ascii(filnam)
% ADCP_ASCII read classic WinRiver TXT ASCII files
% input: filename (string)
% output: ADCP structure
% file structure is based on:
% "WinRiver II User's Guide" (Teledyne RD Instruments, Feb. 2007) pp. 64-65
%
% no error checking is done, so this could fail miserably on malformed
% files
%
% D. Nowacki nowacki@uw.edu Feb. 2012 
%
% Revisions:
% 28 May 2014 - more comments

fid = fopen(filnam);

% these three lines are at the beginning of each file.
adcp.comment1 = fgetl(fid);
adcp.comment2 = fgetl(fid);

format = '%f %f %f %u8 %u8 %u8 %u8';
segarray = textscan(fid, format, 1);

adcp.binsize = cell2mat(segarray(:,1));
adcp.blank = cell2mat(segarray(:,2));
adcp.firstbin = cell2mat(segarray(:,3));
adcp.numcells = cell2mat(segarray(:,4));
adcp.numpings = cell2mat(segarray(:,5));
adcp.tpe = cell2mat(segarray(:,6));
adcp.watermode = cell2mat(segarray(:,7));


% now read blocks throughout file
% note that, to conserve space, ensembles are indexed sequentially, not by
% ensemble number.

n = 1;

while(~feof(fid))

format = '%f %f %f %f %f %f %f %f %f %f %f %f %f';
segarray = textscan(fid, format, 1);
if isempty(cell2mat(segarray(:,8)))
    break
end

% ROW 1
% R Field Description
% 1 1     ENSEMBLE TIME -Year (at start of ensemble) 
%   2                   - Month 
%   3                   - Day 
%   4                   - Hour 
%   5                   - Minute 
%   6                   - Second 
%   7                   - Hundredths of seconds 
%   8     ENSEMBLE NUMBER (or SEGMENT NUMBER for processed or averaged raw  data)   
%   9     NUMBER OF ENSEMBLES IN SEGMENT (if averaging ON or processing data) 
%   10    PITCH – Average for this ensemble (degrees) 
%   11    ROLL – Average for this ensemble (degrees)  
%   12    CORRECTED HEADING - Average ADCP heading (corrected for one cycle error) + heading offset + magnetic variation 
%   13    ADCP TEMPERATURE - Average for this ensemble (°C) 

adcp.year(n) = cell2mat(segarray(:,1));
adcp.month(n) = cell2mat(segarray(:,2));
adcp.day(n) = cell2mat(segarray(:,3));
adcp.hour(n) = cell2mat(segarray(:,4));
adcp.minute(n) = cell2mat(segarray(:,5));
adcp.second(n) = cell2mat(segarray(:,6));
adcp.hunsec(n) = cell2mat(segarray(:,7));
adcp.ensnum(n) = cell2mat(segarray(:,8));
adcp.numens(n) = cell2mat(segarray(:,9));
adcp.pitch(n) = cell2mat(segarray(:,10));
adcp.roll(n) = cell2mat(segarray(:,11));
adcp.heading(n) = cell2mat(segarray(:,12));
adcp.temp(n) = cell2mat(segarray(:,13));

format = '%f %f %f %f %f %f %f %f %f %f %f %f';
segarray = textscan(fid, format, 1);

% ROW 2
% 2	1   BOTTOM-TRACK VELOCITY - East(+)/West(-); average for this ensemble (cm/s or ft/s) 
%   2   Reference = BTM       - North(+)/South(-) 
%   3                         - Vertical (up[+]/down[-]) 
%   4                         - Error 
% 2 1   BOTTOM-TRACK VELOCITY – GPS (GGA or VTG) Velocity (calculated from GGA String) Reference = GGA  East(+)/West (-1) 
%   2   Reference = VTG - GPS (GGA or VTG) North(+)/South(-) Velocity 
%   3                         - BT  (up[+]/down[-]) Velocity
%   4                         - BT Error

adcp.btvel_x(n) = cell2mat(segarray(:,1));
adcp.btvel_y(n) = cell2mat(segarray(:,2));
adcp.btvel_z(n) = cell2mat(segarray(:,3));
adcp.btvel_err(n) = cell2mat(segarray(:,4));
adcp.depthsounder(n) = cell2mat(segarray(:,5));
adcp.ggaalt(n) = cell2mat(segarray(:,6));
adcp.ggadeltaalt(n) = cell2mat(segarray(:,7));
adcp.ggahdop(n) = cell2mat(segarray(:,8));
adcp.depth(1,n) = cell2mat(segarray(:,9));
adcp.depth(2,n) = cell2mat(segarray(:,10));
adcp.depth(3,n) = cell2mat(segarray(:,11));
adcp.depth(4,n) = cell2mat(segarray(:,12));

format = '%f %f %f %f %f';
segarray = textscan(fid, format, 1);

% ROW 3
adcp.elapdist(n) = cell2mat(segarray(:,1));
adcp.elaptime(n) = cell2mat(segarray(:,2));
adcp.distnorth(n) = cell2mat(segarray(:,3));
adcp.disteast(n) = cell2mat(segarray(:,4));
adcp.distmadegood(n) = cell2mat(segarray(:,5));

format = '%f %f %f %f %f';
segarray = textscan(fid, format, 1);

% ROW 4
adcp.lat(n) = cell2mat(segarray(:,1));
adcp.lon(n) = cell2mat(segarray(:,2));
% dont care about other values on this line for now

format = '%f %f %f %f %f %f %f %f %f';
segarray = textscan(fid, format, 1);

% ROW 5
adcp.qmid(n) = cell2mat(segarray(:,1));
adcp.qtop(n) = cell2mat(segarray(:,2));
adcp.qbot(n) = cell2mat(segarray(:,3));
% dont care about other values on this line for now

% ROW 6
% 6 1	NUMBER OF BINS TO FOLLOW 
%   2   MEASUREMENT UNIT – cm or ft
%   3   VELOCITY REFERENCE – BT, GGA, VTG, or NONE for current velocity data rows 7-26 fields 2-7 
%   4   INTENSITY UNITS - dB or counts
%   5   INTENSITY SCALE FACTOR – in dB/count 
%   6   SOUND ABSORPTION FACTOR – in dB/m

format = '%u8 %s %s %s %f %f';
segarray = textscan(fid, format, 1);

adcp.binstofollow(n) = cell2mat(segarray(:,1));
adcp.measunit(n) = segarray(:,2);
adcp.velref(n) = segarray(:,3);
adcp.intensunits(n) = segarray(:,4);
adcp.intensscale(n) = cell2mat(segarray(:,5));
adcp.soundabsorp(n) = cell2mat(segarray(:,6));

% ROW 7-26

% 7-26  1    DEPTH – Corresponds to depth of data for present bin (depth cell); includes ADCP depth and blanking value; in m or ft. 
%       2   VELOCITY MAGNITUDE 
%       3   VELOCITY DIRECTION
%       4   EAST VELOCITY COMPONENT – East(+)/West(-) 
%       5   NORTH VELOCITY COMPONENT - North(+)/South(-) 
%       6   VERTICAL VELOCITY COMPONENT - Up(+)/Down(-) 
%       7   ERROR VELOCITY 
%       8   BACKSCATTER     – Beam 1
%       9                   - Beam 2 
%       10                  - Beam 3 
%       11                  - Beam 4
%       12  PERCENT-GOOD 
%       13  DISCHARGE

format = '%f %f %f %f %f %f %f %f %f %f %f %f %f';
segarray = textscan(fid, format, adcp.numcells);

adcp.z(:,n) = cell2mat(segarray(:,1));
adcp.spd(:,n) = cell2mat(segarray(:,2));
adcp.dir(:,n) = cell2mat(segarray(:,3));
adcp.east(:,n) = cell2mat(segarray(:,4));
adcp.north(:,n) = cell2mat(segarray(:,5));
adcp.up(:,n) = cell2mat(segarray(:,6));
adcp.err(:,n) = cell2mat(segarray(:,7));
adcp.echo(:,n,1) = cell2mat(segarray(:,8));
adcp.echo(:,n,2) = cell2mat(segarray(:,9));
adcp.echo(:,n,3) = cell2mat(segarray(:,10));
adcp.echo(:,n,4) = cell2mat(segarray(:,11));
adcp.percgood(:,n) = cell2mat(segarray(:,12));
adcp.q(:,n) = cell2mat(segarray(:,13));

n = n+1;
end 


% do brief QAQC on the data
adcp.spd(adcp.spd == -32768) = nan;
adcp.east(adcp.east == -32768) = nan;
adcp.north(adcp.north == -32768) = nan;
adcp.up(adcp.up == -32768) = nan;
adcp.err(adcp.err == -32768) = nan;
adcp.lat(adcp.lat == 30000) = nan;
adcp.lon(adcp.lon == 30000) = nan;

adcp.echo(adcp.echo == 255) = nan;

fclose(fid);

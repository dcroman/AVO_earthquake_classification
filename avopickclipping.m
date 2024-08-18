function [Tout]=avopickclipping(mode,ignore,cliplength,FIL1,FIL2,FIH1,FIH2,Fc1,Fc2,pha_file, mseed_file) %
% Reads in AVO .pha files (hypoinverse format) and corresponding .mseed files (1-minute, multiplexed)
% Converts phase file summary line into variables
% Identifies first three P-wave picks
% Reads in corresponding .mseed file and extracts clip from first three
% picked waveforms for analysis
% Creates stacked spectrum from these
% Identifies Frequency Index, top five peaks, and skewness/kurtosis (full
% and band-limited)
% If mode 2 or 3 selected, called eval_result to plot interactive manual
% checking figure
%
%   Last updated Apr 16, 2022 1:19pm by Diana Roman


%   Open the phase file and extract station and p-pick info

warning('off', 'all')

fid = fopen(pha_file);
C=textscan(fid,'%s','delimiter','\n', 'HeaderLines',1); C=char(C{:});
fclose(fid);

fid = fopen(pha_file);
D=textscan(fid,'%s','delimiter','\n'); D=char(D{:});
fclose(fid);

stn=strings(size(C,1)-1,1);
cmp=strings(size(C,1)-1,1);
net=strings(size(C,1)-1,1);
loc=strings(size(C,1)-1,1);
picks=strings(size(C,1)-1,5);

for n=1:size(C,1)-1
if C(n,15) == 'P' && C(n,12) == 'Z' 
stn(n,1)=C(n,1:4);
stn=strtrim(stn);
end
end

for n=1:size(C,1)-1
if C(n,15) == 'P' && C(n,12) == 'Z'
net(n,1)=C(n,6:7);
end
end

for n=1:size(C,1)-1
if C(n,15) == 'P' && C(n,12) == 'Z' 
cmp(n,1)=C(n,10:12);
end
end

for n=1:size(C,1)-1
if C(n,15) == 'P' && C(n,12) == 'Z'
loc(n,1)=C(n,112:113);
end
end

for n=1:size(C,1)-1
if C(n,15) == 'P' && C(n,12) == 'Z'
picks(n,1)=C(n,18:21);
end
end

for n=1:size(C,1)-1
if C(n,15) == 'P' && C(n,12) == 'Z'
picks(n,2)=day(datetime(C(n,18:25), 'InputFormat', 'yyyyMMdd'), 'dayofyear');
end
end

for n=1:size(C,1)-1
if C(n,15) == 'P' && C(n,12) == 'Z'
picks(n,3)=C(n,26:27);
end
end

for n=1:size(C,1)-1
if C(n,15) == 'P' && C(n,12) == 'Z' 
picks(n,4)=C(n,28:29);
end
end

for n=1:size(C,1)-1
if C(n,15) == 'P' && C(n,12) == 'Z'
picks(n,5)=str2num(C(n,31:34))/100;
end
end

oyear=str2num(D(1,1:4));
omonth=str2num(D(1,5:6));
oday=str2num(D(1,7:8));
ohour=str2num(D(1,9:10));
ominute=str2num(D(1,11:12));
osecond=str2num(D(1,13:16))/100;
origintime=datetime(oyear,omonth,oday,ohour,ominute,osecond);
evlat=str2num(D(1,17:18))+(str2num(D(1,20:23))/100/60);
evlon=str2num(D(1,24:26))+(str2num(D(1,28:31))/100/60);
if D(1,27) == ' ' 
    evlon=evlon*-1;
end
evdep=str2num(D(1,32:36))/100;
mag=str2num(D(1,148:150))/100;
azgap=str2num(D(1,43:45));
dist=str2num(D(1,46:48))/100;
rms=str2num(D(1,49:52))/100;
herr=str2num(D(1,86:89))/100;
derr=str2num(D(1,90:93))/100;
evid=str2num(D(1,137:146)); 

picks=str2double(picks);
rows = any(isnan(picks),2);
picks(rows,:) = [];

%Fix magnitudes using MLtrue table
%load('MLtrue.mat')
%MLtrueindex=find(ML{:,:}==evid);
%mag=ML{MLtrueindex,2};


for n=1:size(picks,1)
% Here's the start of a solution to Error 6
    if picks(n,5)>60
        picks(n,5)=picks(n,5)-60;
        picks(n,4)=picks(n,4)+1;
    end

Pick(n,:)=horzcat(num2str(picks(n,1)), '/', num2str(picks(n,2),'%03.f'), '/', num2str(picks(n,3),'%02.f'), '/', num2str(picks(n,4),'%02.f'), '/', num2str(picks(n,5), '%05.2f'));
pick_datetime = datetime(Pick,'InputFormat','yyyy/DDD/HH/mm/ss.SS');
end

%   Get list of picked components with which to query mseed ChannelFullName
 picklist=strings(size(picks,1),1);

for n=1:1:size(stn,1)
 if C(n,15) == 'P'  && C(n,12) == 'Z'
    picklist(n,1)=strcat(net(n),':',stn(n),':',loc(n),':',cmp(n));
    
end
end

% Convert to char and clean up picklist
picklist(strcmp(picklist(:, 1), ':::'), :) = [];
picklist=strrep(picklist,' ','');
picklist=strrep(picklist,'--','');
picklist=picklist(~strcmp(picklist(:,1),""),:);
station=stn(~strcmp(stn(:,1),""),:); 
picklist=rmmissing(picklist);

% Add datetime to picklist and station and sort by first arrival (so that only first
% three arrivals will be used to create stacked periodogram
picklist_sorted=strings(size(picks,1),2);
station_sorted=strings(size(station,1),2);
for n=1:1:size(picklist,1)
   picklist_sorted(n,1)=picklist(n,1);
   picklist_sorted(n,2)= pick_datetime(n,1);
end 
for n=1:1:size(station,1)
   station_sorted(n,1)=station(n,1);
   station_sorted(n,2)= pick_datetime(n,1);
end 
picklist_sorted=sortrows(picklist_sorted, 2);
station_sorted=sortrows(station_sorted, 2);

% Remove any stations on the 'ignore' list from picklist_sorted and
% station_sorted
for n=1:size(picklist_sorted,1)
    TF = contains(ignore,char(station_sorted(n)));
    if TF==1
         picklist_sorted(n,:) = [];
         station_sorted(n,:) = [];
    else
    end
end

% Read in mseed file
[X,I]=rdmseed(mseed_file);

 for n=1:3
    match = reshape(strcmp({I.ChannelFullName}, picklist_sorted(n)), size(I));
    k = find(match, 1, 'first');
    stationlab = station_sorted(n);
    stationlab=strrep(stationlab,' ','');
    try  
    q=I(k).XBlockIndex;  
    picked_waves.(stationlab) = cat(1,X(q).d);    
        catch ME
    fprintf('file %s :', pha_file)
    fprintf('station %s ', stationlab)
    fprintf('not processed: %s\n', ME.message);
    %continue
    end
end

for n=1:numel(fieldnames(picked_waves))  
    try
    match = reshape(strcmp({X.ChannelFullName}, picklist_sorted(n)), size(X));
    h(n,1) = find(match, 1, 'first');
    try
        start_times(n,1:5)=X(h(n)).RecordStartTime;  %this has to be sequence 1
    sample_rate(n,1)=X(h(n)).SampleRate;
    orig_sample_rate(n,1)=sample_rate(n,1);
         catch ME
    fprintf('file %s :', pha_file)
    fprintf('station %s ', stationlab)
    fprintf('not processed: %s\n', ME.message);
    continue
    end
    end
end 


%Generate differences in seconds between start and pick times

for n=1:size(start_times,1)
Start(n,:)=horzcat(num2str(start_times(n,1)), '/', num2str(start_times(n,2)), '/', num2str(start_times(n,3)), '/', num2str(start_times(n,4)), '/', num2str(start_times(n,5), '%05.2f'));
start_datetime = datetime(Start,'InputFormat','yyyy/DDD/HH/mm/ss.SSS');
end

% Downsample 100 Hz channels to 50 Hz and adjust sample rates - kludges
% added to account for noninteger sample rates
for n=1:numel(fieldnames(picked_waves))
    if sample_rate(n,1)-100<1 && sample_rate(n,1)-100>=0
            stationlab = strcat(station_sorted(n));
            stationlab=strrep(stationlab,' ','');
            picked_waves.(stationlab) = downsample(picked_waves.(stationlab),2);
            sample_rate(n) = round(sample_rate(n)/2);
    end
end

%Bandpass filter the picked waveforms
for n=1:numel(fieldnames(picked_waves))
    stationlab = strcat(station_sorted(n));
    stationlab=strrep(stationlab,' ','');
    demeaned.(stationlab)=detrend(picked_waves.(stationlab),0);  
    detrended.(stationlab)=detrend(demeaned.(stationlab),1);
    tapered.(stationlab)=detrended.(stationlab).*hamming(length(detrended.(stationlab)));
    %Design second-order Bu bandpass filter at 0.5-20 Hz
    Fs = 50;  % Sampling Frequency
    N   = 2;    % Order
    % Construct causal 2nd order butterworth filter
    Fc1n=Fc1/(Fs/2);
    Fc2n=Fc2/(Fs/2);
    [B,A]=butter(N,[Fc1n Fc2n]);   
    % Bandpass the clipped wave
    filtered.(stationlab)=filter(B,A,tapered.(stationlab));
end


%Clip out wave from the pick time to specified number of seconds after pick
%Sort pick_datetime now

% Error 3 happens in here
%Problem if pick is late in the waveform - not enough to clip
%Need to generate equal-length clips for frequency indexing. 
%Error is "event not processed: Index exceeds the number of array elements

pick_datetime_sorted=sortrows(pick_datetime, 1);
for n=1:size(start_datetime,1)
    Duration(n)=pick_datetime_sorted(n)-start_datetime(n);           
    if Duration(n)<0; Duration(n) = 0; end
    startclip(n) = seconds(Duration(n))*sample_rate(n);
    if Duration(n)==0; startclip(n)=startclip(n)+1; end;
    endclip(n) = startclip(n)+(sample_rate(n)*cliplength);
    stationlab = strcat(station_sorted(n));
    stationlab=strrep(stationlab,' ','');

    clipped_waves.(stationlab) = filtered.(stationlab)(startclip(n):endclip(n));
end

    
%Now compute spectral parameters for each wave in clipped_waves following JP's GS_STACK code

for n=1:numel(fieldnames(clipped_waves))
    stationlab = strcat(station_sorted(n));
    stationlab=strrep(stationlab,' ','');
    N=length(clipped_waves.(stationlab));
    FR=1/N*0.01;
    K=(1/N)*fft(clipped_waves.(stationlab));
    K(1)=[];
    POWER=2*((abs(K(1:floor(N/2))).^2)/FR);
    M=max(POWER);
    power.(stationlab)=sqrt(POWER/M);
end

% Stack all of the power spectra - need to zero pad to length of longest
% but should not have to do this if lengths are equal (why aren't they?)
fns=fieldnames(power);
first_fns = char(fns(1));
stacklength = max(structfun(@numel,power));
stack = zeros(stacklength,1);

for n = 1:numel(fns)
    if sample_rate(n) == 50  %Hardcoded to toss out 20 Hz sampled channels
    stationlab = strcat(station_sorted(n));
    stationlab = strrep(stationlab,' ','');
    stack=stack + power.(stationlab);
    end
end

Nstack=length(clipped_waves.(stationlab));
nyquist=0.5*sample_rate(1,1);
f=(1:Nstack/2)/(Nstack/2)*nyquist;
  

FIL1s=FIL1*5;
FIL2s=FIL2*5;
FIH1s=FIH1*5;
FIH2s=FIH2*5;

FIBandStack=stack(FIL1s:FIH2s);
fBandStack=f(FIL1s:FIH2s);


%Determine top five peak frequencies
stack_copy=FIBandStack;
firstxIndex = find(stack_copy == max(stack_copy), 1, 'first');
firstpeakfreq = fBandStack(firstxIndex);
stack_copy(firstxIndex) = -Inf;
%
secondxIndex = find(stack_copy == max(stack_copy), 1, 'first');
secondpeakfreq = fBandStack(secondxIndex);
stack_copy(secondxIndex) = -Inf;
%
thirdxIndex = find(stack_copy == max(stack_copy), 1, 'first');
thirdpeakfreq = fBandStack(thirdxIndex);
stack_copy(thirdxIndex) = -Inf;
%
fourthxIndex = find(stack_copy == max(stack_copy), 1, 'first');
fourthpeakfreq = fBandStack(fourthxIndex);
stack_copy(fourthxIndex) = -Inf;
%
fifthxIndex = find(stack_copy == max(stack_copy), 1, 'first');
fifthpeakfreq = fBandStack(fifthxIndex);
stack_copy(fifthxIndex) = -Inf;

%Determine Frequency Index a la Eq. 1 of Buurman and West (2010) using specified bands
%Note Rodgers et al. 2015 uses base-2 log to create meaningful ratios.
% Error 7 is here - stacklength does not equal 125 (it is 63)

if stacklength == cliplength*25
A_upper1=stack(FIH1s:FIH2s);
A_lower1=stack(FIL1s:FIL2s);
FI15v610 = log10(mean(A_upper1)/mean(A_lower1));
else     fprintf('stacklength %i :', stacklength)
end

firstpick=station_sorted(1);
secondpick=station_sorted(2);
thirdpick=station_sorted(3);

%Skewness and kurtosis
%Positive skewness is a dominance of lower freqs
skew=skewness(stack);
skewBL=skewness(FIBandStack);
%Positive kurtosis is a more peaked spectrum
kurt=kurtosis(stack);
kurtBL=kurtosis(FIBandStack); 

if mode==1
    manualcheck=0;
else 
    manualcheck=eval_result(mode,B,A,cliplength,FIL1,FIL2,FIH1,FIH2,picked_waves,firstpick,secondpick,thirdpick,startclip,endclip,sample_rate,f,stack,firstpeakfreq,secondpeakfreq,thirdpeakfreq,fourthpeakfreq,fifthpeakfreq,evid,origintime,evdep,kurt,kurtBL,skew,skewBL,mag,FI15v610,power,picklist_sorted);       
    %manualcheck = menu('Manual Classification', 'HF', 'LF', 'Uncertain');
    %close
end

%Create a table line to add to cumulative table
Tout=table(origintime,evlat,evlon,evdep,mag,azgap,dist,rms,herr,derr,evid,firstpeakfreq,secondpeakfreq,thirdpeakfreq,fourthpeakfreq,fifthpeakfreq,FI15v610,kurt,skew, kurtBL, skewBL, firstpick, secondpick, thirdpick, manualcheck);
clear picked_waves clipped_waves Start start_datetime fns start_times filtered Duration power
end

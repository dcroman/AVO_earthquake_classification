%AutoclassMS 
%   Code by Diana Roman to parse AVO pha files for P-picks and pull corresponding waveforms from mseed.
%   Combined with algorithm from JP's GS_STACK example to produce a stacked spectrum for all P-picked Z channels
%   Extracts frequency index (FI), top five peaks, and kurtosis/skewness
%   from stacked spectrum
%
%   Last updated Apr 16, 2022 1:19pm by Diana Roman




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%SET THE PARAMETERS FOR ANALYSIS HERE: 

% Set bands for Frequency Index (default is 1-5 and 6-10) - note they must
% be equal in width
    FIL1=1;
    FIL2=5;
    FIH1=6;
    FIH2=10;

% Set length of data after P-wave pick to analyze (default is 5 seconds)
    cliplength=5;

% Set bandpass filter to apply before analysis (default is 0.5-20 Hz)
    Fc1 = 0.5;  % First Cutoff Frequency
    Fc2 = 20;   % Second Cutoff Frequency
    
% Stations to ignore (needs the trailing whitespace to work)
   ignore = 'AKXX ';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% DO NOT EDIT BELOW THIS LINE

%Prompt user to set mode: 
%Mode 1: Analysis (does not make any plots)
%
%Mode 2: Calibration (use on subset to manually select clear LP/VT and plot
%Frequency Index by event class)
%
%Mode 3: Blind Calibration (use on subset to manually select clear LP/VT and plot
%Frequency Index by event class) - does not show autoclassification info in figures

mode=1;

prompt = "Set Mode\n1 - Analysis Mode\n2 - Calibration Mode (show autoclassification)\n3 - Calibration Mode (blind)\n";
mode=input(prompt);

% Create empty table with headers for output
phafiles = dir('*.pha');
seedfiles = dir('*.mseed');
headers = {'origintime' 'evlat' 'evlon' 'evdep' 'mag' 'azgap' 'dist' 'rms' 'herr' 'derr' 'evid' 'firstpeakfreq' 'secondpeakfreq' 'thirdpeakfreq' 'fourthpeakfreq' 'fifthpeakfreq' 'FI15v610' 'kurt' 'skew' 'kurtBL' 'skewBL' 'firstpick' 'secondpick' 'thirdpick' 'manualcheck'};
data = cell(0,25);
Table = cell2table(data);
Table.Properties.VariableNames = headers;
clear data headers

% Loop over all .pha-.mseed pairs in directory to produce output Table
for n=1:length(phafiles)
    try
    pha_file=phafiles(n).name;
    mseed_file=seedfiles(n).name;
    [Tout]=avopickclipping(mode,ignore,cliplength,FIL1,FIL2,FIH1,FIH2,Fc1,Fc2,pha_file,mseed_file); 
  Table=[Table; Tout];
    catch ME
    fprintf('file %s :', pha_file)
    fprintf('event not processed: %s\n', ME.message);
    continue
    end
end

clear mode ignore Tout cliplength ME mseed_file n pha_file prompt Fc1 Fc2 FIH1 FIH2 FIL1 FIL2 phafiles seedfiles


%  Vestigal code - for future version
%if mode==4;
%make a plot like Fig x of Buurman and West
%FI=Table{:,18};
%Manclass=Table{:,27};
%plot(Manclass,FI,'+')
%xlabel('Manual Class: 1-HF, 2-LF, 3-uncertain' )
%ylabel('Frequency Index')
%xlim([0 4])
%end

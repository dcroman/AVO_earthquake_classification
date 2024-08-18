function manualcheck=eval_result(mode,B,A,cliplength,FIL1,FIL2,FIH1,FIH2,picked_waves,firstpick,secondpick,thirdpick,startclip,endclip,sample_rate,f,stack,firstpeakfreq,secondpeakfreq,thirdpeakfreq,fourthpeakfreq,fifthpeakfreq,evid,origintime,evdep,kurt,kurtBL,skew,skewBL,mag,FI15v610,power,picklist_sorted)
%
%   Last updated Apr 16, 2022 1:19pm by Diana Roman

%Filter first waveform using specified band (just reapplies butterworth
%filter from avopickclipping
tmp=picked_waves.(firstpick);
demean_tmp=detrend(tmp,0);
detrend_tmp=detrend(demean_tmp,1);
taper_tmp=detrend_tmp.*hamming(length(detrend_tmp));
filt_tmp_first=filter(B,A,taper_tmp);

%Cut a small chunk of filt_temp - enough for spectrogram (startclip-4s (=200 samples) to startclip+10 = 500 samples)
bitstart1=startclip(1)-200;
bitend1=startclip(1)+500;
bitfirst=filt_tmp_first(bitstart1:bitend1);

%First Spectrogram
window = 256;
nfft = 256;
overlap = 128;
[sfirst,fofirst,tfirst] = spectrogram(bitfirst,window,overlap,nfft,sample_rate(1));
sfirst = (abs(sfirst)).^2;

%Filter second waveform
tmp=picked_waves.(secondpick);
demean_tmp=detrend(tmp,0);
detrend_tmp=detrend(demean_tmp,1);
taper_tmp=detrend_tmp.*hamming(length(detrend_tmp));
%filt_tmp_second=filter(Hd,taper_tmp);
filt_tmp_second=filter(B,A,taper_tmp);

%Cut a small chunk of filt_temp - enough for spectrogram (startclip-4s (=200 samples) to startclip+10 = 500 samples)
bitstart2=startclip(2)-200;
bitend2=startclip(2)+500;
bitsecond=filt_tmp_second(bitstart2:bitend2);

%Second Spectrogram
window = 256;
nfft = 256;
overlap = 128;
[ssecond,fosecond,tsecond] = spectrogram(bitsecond,window,overlap,nfft,sample_rate(2));
ssecond = (abs(ssecond)).^2;

%Filter third waveform
tmp=picked_waves.(thirdpick);
demean_tmp=detrend(tmp,0);
detrend_tmp=detrend(demean_tmp,1);
taper_tmp=detrend_tmp.*hamming(length(detrend_tmp));
%filt_tmp_third=filter(Hd,taper_tmp);
filt_tmp_third=filter(B,A,taper_tmp);


%Cut a small chunk of filt_temp - enough for spectrogram (startclip-4s (=200 samples) to startclip+10 = 500 samples)
bitstart3=startclip(3)-200;
bitend3=startclip(3)+500;
bitthird=filt_tmp_third(bitstart3:bitend3);

%Third Spectrogram
window = 256;
nfft = 256;
overlap = 128;
[sthird,fothird,tthird] = spectrogram(bitthird,window,overlap,nfft,sample_rate(3));
sthird = (abs(sthird)).^2;
    

%%Make figure
mcfig=figure('units','normalized','outerposition',[0 0 1 1]);

if mode~=1
%stacked spectrum
subplot(6,3,[1 4]);
plot(f,stack);
hold('on');
ylabel('Normalized Count Spectrum')
xlabel('Hz')
xlim([0 25]);
ylim([0 5])
title('Stacked Spectrum')
rectangle('Position',[FIL1 0 (FIL2-FIL1) max(f)], 'FaceColor', [1 0 0 0.2]);
rectangle('Position',[FIH1 0 (FIH2-FIH1) max(f)], 'FaceColor', [0 0 1 0.2]);
else
end

if mode==2
meanoffive=mean([firstpeakfreq secondpeakfreq thirdpeakfreq fourthpeakfreq fifthpeakfreq]);
labelx= [firstpeakfreq secondpeakfreq thirdpeakfreq fourthpeakfreq fifthpeakfreq meanoffive];
labely= [2 2 2 2 2 2.5];
plot(labelx,labely,'+', 'MarkerSize', 10)
labels = {'1','2','3','4','5','avg'};
text(labelx,labely,labels,'VerticalAlignment','top','HorizontalAlignment','center')
box('on');
grid('on')
hold('off');
else
end

if mode==2
ax=subplot(6,3,[7 10]);
%textbox
%dim = [0.6 0.7 0.15 0.15];
%dimclass = [0.6 0.65 0.15 0.15];
str = {append('CUSPID: ', num2str(evid)), append('Time: ', datestr(origintime)), append('Depth: ',num2str(evdep), ' km BSL'), append('Kurt-Full: ', num2str(kurt)), append('Kurt-BL: ', num2str(kurtBL)), append('Skew-Full: ', num2str(skew)), append('Skew-BL: ', num2str(skewBL)), append('ML: ', num2str(mag)), append('FI: ', num2str(FI15v610))};
text(0,0.5, str,'FontSize',16);
set (ax, 'visible', 'off')
if FI15v610<0 autoclass= "Low Frequency";
    text(0,0, autoclass,'FontSize',16, 'Color',[1 1 1],'BackgroundColor',[1 0 0],'FontWeight','bold', 'Position',[0.695412844036697 0.771084337349398 0]);
    else autoclass="High Frequency" ;
            text(0,0,autoclass,'FontSize',16, 'Color',[1 1 1],'BackgroundColor',[0 0 1],'FontWeight','bold','Position',[0.695412844036697 0.771084337349398 0]);
end


annotation(mcfig,'rectangle',...
    [0.0956051602814699 0.416543574593796 0.275778733385457 0.539881831610045],...
    'LineWidth',2);
annotation(mcfig,'textbox',...
    [0.0956051602814699 0.95642540620384 0.134261923377639 0.0228951255539144],...
    'String',{'Automatic Classification:'},...
    'LineStyle','none',...
    'FontWeight','bold',...
    'FontSize',14,...
    'FitBoxToText','off');
else
end

subplot(6,3,[2 5]);
%Plot first waveform PSD
hold('on');
PSDfirst=power.(firstpick);
plot(f,PSDfirst)
ylabel('Normalized Count Spectrum')
xlabel('Hz')
ylim([0 5])
title(picklist_sorted(1))
rectangle('Position',[FIL1 0 (FIL2-FIL1) max(f)], 'FaceColor', [1 0 0 0.2]);
rectangle('Position',[FIH1 0 (FIH2-FIH1) max(f)], 'FaceColor', [0 0 1 0.2]);
box('on');
hold('off'); 

subplot(6,3,3);
%Plot first filtered waveform
hold('on');
set(gca,'XTick',...
    [0 50 100 150 200 250 300 350 400 450 500 550 600 650 700],'XTickLabel',...
    {'0','1','2','3','4','5','6','7','8','9','10','11','12','13','14'})
xlim([100 550])
plot(bitfirst)
xlabel('Time, s')
ylabel('Amplitude')
xline(200, 'LineWidth', 2);
xline(200+(50*cliplength), 'LineWidth', 2);
title(picklist_sorted(1))
box('on');
hold('off');

subplot(6,3,6);
%Plot first Spectrogram - logarithmic
hold('on');
contourf(tfirst,fofirst,log10(sfirst));
xlabel('Time, s')
ylabel('Magnitude, dB')
xlim([2 11])
ylim([0 25]);
title(picklist_sorted(1))
box('on');
hold('off');

subplot(6,3,[8 11]);
%Plot second waveform PSD
hold('on');
PSDsecond=power.(secondpick);
plot(f,PSDsecond)
ylabel('Normalized Count Spectrum')
xlabel('Hz')
title({'';picklist_sorted(2)})
ylim([0 5])
rectangle('Position',[FIL1 0 (FIL2-FIL1) max(f)], 'FaceColor', [1 0 0 0.2]);
rectangle('Position',[FIH1 0 (FIH2-FIH1) max(f)], 'FaceColor', [0 0 1 0.2]);
box('on');
hold('off');

subplot(6,3,9);
%Plot second filtered waveform
hold('on');
plot(bitsecond);
set(gca,'XTick',...
    [0 50 100 150 200 250 300 350 400 450 500 550 600 650 700],'XTickLabel',...
    {'0','1','2','3','4','5','6','7','8','9','10','11','12','13','14'})
xlim([100 550])
xlabel('Time, s')
ylabel('Amplitude')
title({'';picklist_sorted(2)})
xline(200, 'LineWidth', 2);
xline(200+(50*cliplength), 'LineWidth', 2);
box('on');
hold('off');

subplot(6,3,12);
%Plot second Spectrogram - logarithmic
hold('on');
contourf(tsecond,fosecond,log10(ssecond));
xlabel('Time, s')
ylabel('Magnitude, dB')
title({'';picklist_sorted(2)})
xlim([2 11])
ylim([0 25]);
box('on');
hold('off');

subplot(6,3,[14 17]);
%Plot third waveform PSD
hold('on');
PSDthird=power.(thirdpick);
plot(f,PSDthird)
ylabel('Normalized Count Spectrum')
xlabel('Hz')
title({'';picklist_sorted(3)})
ylim([0 5])
rectangle('Position',[FIL1 0 (FIL2-FIL1) max(f)], 'FaceColor', [1 0 0 0.2]);
rectangle('Position',[FIH1 0 (FIH2-FIH1) max(f)], 'FaceColor', [0 0 1 0.2]);
box('on');
hold('off');

subplot(6,3,15);
%Plot Third filtered waveform
hold('on');
plot(bitthird);
set(gca,'XTick',...
    [0 50 100 150 200 250 300 350 400 450 500 550 600 650 700],'XTickLabel',...
    {'0','1','2','3','4','5','6','7','8','9','10','11','12','13','14'})
xlim([100 550])
xlabel('Time, s')
ylabel('Amplitude')
title(picklist_sorted(3))
xline(200, 'LineWidth', 2);
xline(200+(50*cliplength), 'LineWidth', 2);
box('on');
hold('off');

subplot(6,3,18);
%Plot third Spectrogram - logarithmic
hold('on');
contourf(tthird,fothird,log10(sthird));
xlabel('Time, s')
ylabel('Magnitude, dB')
title(picklist_sorted(3))
xlim([2 11])
ylim([0 25]);
box('on');
hold('off');


%Add buttons to set manualcheck variable (1=HF,    2=LF,    3=uncertain)
hfb=uicontrol('style','push','units','pix','position',[450 400 300 40],'fontsize',14,'string','High-Frequency','callback',{@hfbutton_call, mcfig});
lfb=uicontrol('style','push','units','pix','position',[450 350 300 40],'fontsize',14,'string','Low-Frequency','callback',{@lfbutton_call,mcfig});
ufb=uicontrol('style','push','units','pix','position',[450 300 300 40],'fontsize',14,'string','Uncertain','callback',{@uncbutton_call,mcfig});
annotation(mcfig,'rectangle',...
    [0.0953870211102424 0.194977843426883 0.275605942142299 0.152141802067947],...
    'LineWidth',2);
annotation(mcfig,'textbox',...
    [0.0959960906958561 0.350812407680945 0.119406567630962 0.0206794682422451],...
    'String',{'Manual Classification:'},...
    'LineStyle','none',...
    'FontWeight','bold',...
    'FontSize',14,...
    'FitBoxToText','off');
uiwait(mcfig)
%disp('waiting for click...');
delete(mcfig)


    function hfbutton_call(hfb,EventData,mcfig)
        manualcheck=1;
        uiresume(mcfig)
    end
    
    function lfbutton_call(lfb,EventData,mcfig)
        manualcheck=2;
        uiresume(mcfig)
    end

    function uncbutton_call(ufb,EventData,mcfig)
        manualcheck=3;
        uiresume(mcfig)
    end

%Vestigal code - for future version
%Plots a PDF that doesn't look very nice - disabled for now
%set(gcf, 'PaperUnits', 'inches');
%x_width=8.5; y_width=11;
%set(gcf, 'PaperPosition', [0 0 x_width y_width]); 
%filename = [num2str(evid) '.pdf'];
%saveas(gcf,filename)

end

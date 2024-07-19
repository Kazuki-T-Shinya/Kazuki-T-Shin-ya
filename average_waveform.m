%%%% waveform average

clear
clc
close all force

set(0,'defaultAxesFontSize',20);
set(0,'defaultAxesFontName','Arial');
set(0,'defaultTextFontSize',20);
set(0,'defaultTextFontName','Arial');

% [Fname_freq,Pname_freq]=uigetfile('*.mat','Results','C:\Users\Takafumi\Desktop\ユビ');
frequency=[2;4;8;16;32;50;60];

[Fname,Pname]=uigetfile('*.mat','Results','C:\Users\');
mkdir(Pname,strtok(Fname(1:4)))

[pathstr,name,ext] = fileparts(Fname);

DATA=load([Pname,Fname]);

waveform=DATA.amplifier_data';
t_amplifier=DATA.t_amplifier';
sample_rate=DATA.sample_rate;
trigger_data=DATA.board_adc_data;

% DATA=load([Pname,Fname]);
% % token1=strtok(Fname(2:4);
% token_Ch1 =['V',strtok(Fname(1:4)),'_Ch1'];
% token_Ch2 =['V',strtok(Fname(1:4)),'_Ch2'];

% Ch1=getfield(DATA, token_Ch1);
% Ch2=getfield(DATA, token_Ch2);
% 
% trig_time=0:Ch1.interval:Ch1.length*Ch1.interval;
% fs=round(1/Ch1.interval);
% Resampling_values=resample(Ch1.values,round(1/Ch2.interval),fs);
% 
% % plot(DATA.V60dB_Ch1.values);
% % hold on
% % plot(Resampling_values,'r');
% % hold off
% 
% waveform_time=0:Ch2.interval:Ch2.length*Ch2.interval;
% fs1=round(1/Ch2.interval);
% 
% if length(trig_time(1:end))==length(Ch1.values)
%     trig=[trig_time(1:end)',Ch1.values];
%     
% else
%     trig=[trig_time(1:end-1)',Ch1.values];
%     
% end
% 
% if length(waveform_time(1:end))==length(Resampling_values)
%     trig_2=[waveform_time(1:end)',Resampling_values(1:end)];
%     
% else
%     trig_2=[waveform_time(1:end-(length(waveform_time)-length(Resampling_values)))',Resampling_values(1:end)];
%     
% end
% 
% 
% if length(waveform_time(1:end))==length(Ch2.values)
%     waveform=[waveform_time(1:end)',Ch2.values];
%     
% else
%     waveform=[waveform_time(1:end-1)',Ch2.values];
%     
% end



% n=10;
% LPF_freq=100;
% HPF_freq=0.1;
% LPF_freq=LPF_freq/((1/Ch2.interval)/2);
% HPF_freq=HPF_freq/((1/Ch2.interval)/2);
% LPF_freq=fir1(n,[HPF_freq,LPF_freq]);
% LPF_trig=filtfilt(LPF_freq,1,trig_2(1:end,2));

% plot(LPF_trig);
trigger_threshold=1.5;

upper_threshold=find(trigger_data<trigger_threshold);
trigger_dulation=diff(upper_threshold);
trigger_location=find(trigger_dulation>sample_rate*0.1);
trigger_location_2=[1;trigger_location'+1];
trigger_location=[1;trigger_location'];
upper_threshold=upper_threshold';

for i=1:length(trigger_location)-1
   trigger_dulation_2(i,1)=upper_threshold(trigger_location(i+1,1),1)-upper_threshold(trigger_location(i,1)+1,1); 
end

trigger_dulation_2=[trigger_dulation_2;upper_threshold(length(upper_threshold),1)-upper_threshold(trigger_location(length(trigger_location),1)+1,1)];
trigger_dulation_2=trigger_dulation_2/sample_rate*1000;

figure;
plot(t_amplifier,trigger_data,'b');
hold on
plot(t_amplifier(upper_threshold(trigger_location_2(1:end)),1),trigger_threshold,'+');
% plot(t_amplifier(upper_threshold_deviant(trigger_location_2_deviant(1:end),1),1),trigger_threshold_deviant,'+');
% plot(trig(upper_threshold_deviant(trigger_location_2_deviant(1:end),1),1),trigger_threshold_deviant,'+');
hold off

% figure;
% plot(trig(1:end,1),DATA.V60dB_Ch1.values,'b');
% hold on
% plot(trig(upper_threshold(trigger_location_2(1:end),1),1),-0.03,'+');
% hold off

% figure;
% plot(trig_2(1:end,1),trig_2(1:end,2),'b');
% hold on
% plot(trig_2(upper_threshold(trigger_location_2(1:end),1),1),-0.03,'+');
% hold off

%%
trigger_threshold=(max(trigger_data)+min(trigger_data))/2;

part_trig=ones(93000,length(trigger_location_2));
part_waveform=ones(93000,length(trigger_location_2));

n=100;
LPF_freq=3000;
HPF_freq=300;
LPF_freq=LPF_freq/(sample_rate/2);
HPF_freq=HPF_freq/(sample_rate/2);
LPF_freq=fir1(n,[HPF_freq,LPF_freq]);
waveform1=[waveform(1:end,1),filtfilt(LPF_freq,1,waveform(1:end,1))];

for j=1:length(trigger_location_2)
    
    if j<length(frequency)
        
        part_trig(1:upper_threshold(trigger_location_2(j+1,1))-upper_threshold(trigger_location_2(j,1))+1,j)=trigger_data(upper_threshold(trigger_location_2(j,1),1):upper_threshold(trigger_location_2(j+1,1),1));
        part_waveform(1:upper_threshold(trigger_location_2(j+1,1))-upper_threshold(trigger_location_2(j,1))+1,j)=waveform1(upper_threshold(trigger_location_2(j,1)):upper_threshold(trigger_location_2(j+1,1)),2);
    
    else
        
        part_trig(1:length(trigger_data)-upper_threshold(trigger_location_2(j,1))+1,j)=trigger_data(upper_threshold(trigger_location_2(j,1),1):end);
        part_waveform(1:length(waveform1)-upper_threshold(trigger_location_2(j,1))+1,j)=waveform1(upper_threshold(trigger_location_2(j,1),1):end,2); 
    
    end
end
%%
% for j=1:27

%     parttrig=[trig_2(upper_threshold(trigger_location_2(j,1),1):upper_threshold(trigger_location_2(j+1,1),1),1),trig_2(upper_threshold(trigger_location_2(j,1),1):upper_threshold(trigger_location_2(j+1,1),1),2)];
%     partwaveform=[waveform(upper_threshold(trigger_location_2(j,1),1):upper_threshold(trigger_location_2(j+1,1),1),1),waveform(upper_threshold(trigger_location_2(j,1),1):upper_threshold(trigger_location_2(j+1,1),1),2)];
%     
ave_waveform_data_include=zeros(round(sample_rate*0.05),length(frequency));

for j=1:length(frequency)
%     j=1;
    upper_threshold_trig=find(part_trig(1:end,j)<trigger_threshold);
    trigger_dulation_trig=diff(upper_threshold_trig);
    trigger_location_trig=find(trigger_dulation_trig>sample_rate*0.001);
    trigger_location_trig_2=[1;trigger_location_trig+1];
    trigger_location_trig=[1;trigger_location_trig];

    time=0:(1/sample_rate):894*(1/sample_rate);
    time=(time-0.01)*1000;
    x=time(1:893);
    % partwaveform=[waveform(upper_threshold(trigger_location_2(1,1),1):upper_threshold(trigger_location_2(2,1),1),1),waveform(upper_threshold(trigger_location_2(1,1),1):upper_threshold(trigger_location_2(2,1),1),2)];

    time_trig=0:(1/sample_rate):length(part_trig(1:end,1))/sample_rate;
    time_trig=time_trig(1:end-1)';
    figure;
    plot(time_trig,part_trig(1:end,1));
%     plot(part_trig(1:end,1),parttrig(1:end,2),'b');
    hold on
    plot(time_trig(upper_threshold_trig(trigger_location_trig_2(1:end,1),1),1),trigger_threshold,'+');
    hold off

    waveform_data=zeros(893,length(trigger_location_trig_2)-1);

    for i=1:length(trigger_location_trig_2)-2
        partwaveform_2=part_waveform(upper_threshold_trig(trigger_location_trig_2(i,1),1):upper_threshold_trig(trigger_location_trig_2(i+1,1),1),j);
        waveform_data(1:length(partwaveform_2),i)=partwaveform_2(find(partwaveform_2~=0));

    end

    ave_waveform_data=mean(waveform_data,2);
    std_waveform_data=std(waveform_data,0,2);
    
    ave_waveform_data=ave_waveform_data(1:893,1);
    std_waveform_data=std_waveform_data(1:893,1);
    
    Base_ave_waveform_data=max(ave_waveform_data(1:round(sample_rate*0.01)));
    
    ave_waveform_data_include(1:893,j)=ave_waveform_data;
%     figure;
    
%     plot(time(1:893),ave_waveform_data(1:893));
%     H1 = plot(x, ave_waveform_data(1:893), 'Color', 'k', 'LineWidth', 1);
%     hold on;
%     H2 = plot(x, [ave_waveform_data(1:893) - 0.5*std_waveform_data(1:893), ave_waveform_data(1:893) + 0.5*std_waveform_data(1:893)], 'Color', 'r','LineWidth', 2);
%     hold off
%     axis([-10,40,-0.0175,-0.0075]);
    FIG=figure;
%     H1 = plot(time(1:893), ave_waveform_data(1:893), 'Color', 'k', 'LineWidth', 1);
%     H2 = shadedErrorBar(time(1:893), waveform_data, {@mean, @(time(1:893)) 0.5*std(time(1:893))  }, '-m', 0);
    H(2) = shadedErrorBar(x, waveform_data(1:893,1:end)', {@mean, @(x) std(x)/sqrt(100)  }, {'-b', 'LineWidth', 2}, 0);
    axis([-3,15,-5,5]);
%     TI=num2str(frequency(1,j));
    TI=num2str(frequency(j,1));
    re_line1=refline([0,Base_ave_waveform_data+0.25]);
    title(TI)
    xlabel('Time [ms]')
    ylabel('Amplitude [μV]')
    saveas(gcf,[Pname,strtok(Fname(1:4)),'\',TI,'.jpg'])
    hold off
end
% file_name=[strtok(Fname(1:4)),'\total.mat'];
save([Pname,strtok(Fname(1:4)),'\',strtok(Fname(1:4)),'total.mat'], 'ave_waveform_data_include');
% close all force
msgbox('Finish!!!')
% end


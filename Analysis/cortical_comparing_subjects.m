
nSubj = 44;
subj_ID = ["P003","P006", "P007","P009","P012","P013","P016",...
    "P017","P019","P020","P024","P027","P028","P029","P030",...
    "P036","P037","P039","P040","P044","P046","P048","P049",...
    "P051","P053","P054","P056","P062","P064","P065","P066",...
    "P067","P068","P076","P077","P078","P079","P083","P085",...
    "P086","P090","P092","P094","P095"];


fs = 500;
nChan = 14;
nTrials = 480;


for i = 1:nSubj
    % [avg_hi_trial_var(i,:,:), avg_lo_trial_var(i,:,:), itpc_hi(i,:,:), itpc_lo(i,:,:), ~, ~,hz,t,phase_hi(i,:,:), phase_lo(i,:,:),phase_shift_hi(i,:)] = cortical_analysis_simplified(fs,nChan,subj_ID(i),nTrials,chan_used);
    if subj_ID(i) == 'P030'
        chan_used = [1 2 4 5 6 8 9 22 24 25 26 28 30 31]; %14 channels
    else
        chan_used = [1 2 4 5 6 8 9 23 25 26 27 29 31 32]; %14 channels
    end
    [avg_hi_trial_var(i,:,:),avg_lo_trial_var(i,:,:),t,fftmag_lo(i,:),fftmag_hi(i,:),itpc_hi(i,:,:),itpc_lo(i,:,:),phase_hi(i,:,:),phase_lo(i,:,:),avg_phase_both_cond(i,:),~]= single_subject_cortical_analysis(fs,nChan,subj_ID(i),nTrials,chan_used);
end


tones = 0:0.25:1.5;
tones = tones*1000;
lo_tones = tones(1,[1,3,5]);
hi_tones = tones(1,[2,4,6]);


% for each subject (1st dimension), average the cortical, itpc, and phase for each
% electrode (2nd dimension) across the average epoch/length (3rd dim)
for j = 1:nSubj
    %     attn_avg_hi_chan(j,:) = mean(squeeze(avg_hi_chan_var(j,:,:))); %channel variation
    %     attn_avg_lo_chan(j,:) = mean(squeeze(avg_lo_chan_var(j,:,:))); %channel variation
    attn_avg_hi_trial(j,:) = mean(squeeze(avg_hi_trial_var(j,:,:)),2); %trial variation
    attn_avg_lo_trial(j,:) = mean(squeeze(avg_lo_trial_var(j,:,:)),2); %trial variation
    itpc_avg_hi(j,:) = mean(squeeze(itpc_hi(j,:,:)));
    itpc_avg_lo(j,:) = mean(squeeze(itpc_lo(j,:,:)));
    phase_avg_hi(j,:) = mean(squeeze(phase_hi(j,:,:)));
    phase_avg_lo(j,:) = mean(squeeze(phase_lo(j,:,:)));
end

fft_mag_hi_2Hz = fftmag_hi(:,5);
fft_mag_lo_2Hz = fftmag_lo(:,5);
fft_mag_difference = fft_mag_hi_2Hz - fft_mag_lo_2Hz;
 
% figure
% scatter(behav_high_avg_dp,fft_mag_difference)
% ylabel('Magnitude difference at 2Hz')
% xlabel('High D-prime')

% [correlation,P] = corrcoef(behav_high_aavg_dp,fft_mag_difference);

%% Plotting cortical data -- SEM across trials
%%%%%% Plotting cortical responses %%%%%%
addpath errorbar_files/
f1 = figure;
figure(f1)
s_attn = shadedErrorBar(t,mean(attn_avg_hi_trial,1),std(attn_avg_hi_trial)/sqrt(nSubj),'lineProps','r');
hold on
s_attn.mainLine.LineWidth = 1.5;
s_dist = shadedErrorBar(t,mean(attn_avg_lo_trial,1),std(attn_avg_lo_trial)/sqrt(nSubj),'lineProps','b');
s_dist.mainLine.LineWidth = 1.5;
xline(hi_tones(2),'r--', 'High Onset','LabelVerticalAlignment','bottom', 'LineWidth',1.5,'FontWeight','bold');
xline(hi_tones(1,[1 3]),'r--','LineWidth',1.5,'FontWeight','bold');

xline(lo_tones(2),'b--','Low Onset','LabelVerticalAlignment','bottom','LineWidth',1.5,'FontWeight','bold');
xline(lo_tones(1,[1 3]),'b--','LineWidth',1.5,'FontWeight','bold');

xline(tones(1,7),'-',{'End of', 'Sequence'},'LabelVerticalAlignment','bottom','FontWeight','bold','LineWidth',1.5)
legend({'Attn Hi', 'Attn Lo'},'FontSize',14)
set(gca,'FontSize',14,'FontWeight','bold');
hold off
ylabel('Amplitude (microV)','FontSize',14,'FontWeight','bold')
xlabel('Time (ms)','FontSize',14,'FontWeight','bold')
title(['Average Epoch across all Trials (n=' num2str(nSubj) ')'])
xlim([-100 2000])

% exportgraphics(f2,'Cortical_avg.png','Resolution',1000)

%%%%%% Plotting phase shift %%%%%%
for i = 1:nSubj
    % phase_shift(i,:) = phase_avg_hi(i,5) - phase_avg_lo(i,5);
    phase_shift(i,:) = wrapToPi(phase_avg_hi(i,5) - phase_avg_lo(i,5));
end

% phase_shift = wrapToPi(phase_shift);

figure()
polarplot(abs(phase_shift),itpc_avg_hi(:,5),'*b','MarkerSize',13,'MarkerFaceColor','auto','LineWidth',1)
ax = gca;
ax.GridAlpha =0.5;
ax.FontSize = 12;
ax.FontWeight = 'bold';
title(['Average Neural Phase Shift at 2Hz (n=' num2str(nSubj) ')'])

% itpc_both_cond = [itpc_avg_hi(:,5) itpc_avg_lo(:,5)];
% avg_itpc_both_cond = mean(itpc_both_cond,2);




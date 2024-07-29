function  [avg_hi_trialvar,avg_lo_trialvar,t,fftmag_lo,fftmag_hi,itpc_hi,itpc_lo,phase_hi,phase_lo,avg_phase_both_cond,phase_shift]= single_subject_cortical_analysis(fs,nChan,subject_id,nTrials,chan_used)
% % %%%%%%%%%%%% Uncomment below if running single subject
% clear all
% close all

% subject_id = 'P094';
% nChan = 14;
% if subject_id == 'P030'
%     chan_used = [1 2 4 5 6 8 9 22 24 25 26 28 30 31]; %14 channels
% else
%     chan_used = [1 2 4 5 6 8 9 23 25 26 27 29 31 32]; %14 channels
% end
% 
% fs = 500;
% nTrials = 480;
% tones = 0:0.25:1.5;
% tones = tones*1000;
% lo_tones = tones(1,[1,3,5]);
% hi_tones = tones(1,[2,4,6]);

%%%%%%%%%% Run below 
% subject_id = 'P017' only has 477 lo channels --> RE DO 

if subject_id == 'P030'
    attn_Hi = load(['/Users/victoriafigarola/Documents/1_CMU/Barb_Lab/Projects/SASA/Data/SASA_' num2str(subject_id) '/AttnHi100_redo.mat']);
    attn_Lo = load(['/Users/victoriafigarola/Documents/1_CMU/Barb_Lab/Projects/SASA/Data/SASA_' num2str(subject_id) '/AttnLo100_redo.mat']);
else
    attn_Hi = load(['/Users/victoriafigarola/Documents/1_CMU/Barb_Lab/Projects/SASA/Data/SASA_' num2str(subject_id) '/AttnHi100.mat']);
    attn_Lo = load(['/Users/victoriafigarola/Documents/1_CMU/Barb_Lab/Projects/SASA/Data/SASA_' num2str(subject_id) '/AttnLo100.mat']);
end

t = attn_Hi.EEG.times;

%%%%%%%%%%% Compiling Channels Data
for j = 1:nChan %variability across channels 
    ch_hi(j,:,:) = squeeze(attn_Hi.EEG.data(chan_used(1,j),:,:));
    avg_hi_chanvar(j,:) = mean(squeeze(ch_hi(j,:,:)),2); %(14 channels x 1050 time points)

    ch_lo(j,:,:) = squeeze(attn_Lo.EEG.data(chan_used(1,j),:,:));
    avg_lo_chanvar(j,:) = mean(squeeze(ch_lo(j,:,:)),2); %(14 channels x 1050 time points)
end
if subject_id == "P017"
    ch_lo(:,:,478:480) = 0;
end

for k = 1:nTrials %variability across trials 
    if subject_id == 'P017'
        avg_hi_trialvar(:,k) = detrend(mean(squeeze(ch_hi(:,:,k))),2); %(1050 time points x 480 trials)
        for k = 1:477
            avg_lo_trialvar(:,k) = detrend(mean(squeeze(ch_lo(:,:,k))) ,2); %(1050 time points x 477 trials)
            avg_lo_trialvar(:,478:480) = 0;
            %   SEM_low_trial = std(avg_lo_trialvar,[],2)./sqrt(477);
        end
    else
        avg_hi_trialvar(:,k) = detrend(mean(squeeze(ch_hi(:,:,k))) ,2); %(1050 time points x 480 trials)
        avg_lo_trialvar(:,k) = detrend(mean(squeeze(ch_lo(:,:,k))) ,2); %(1050 time points x 480 trials)
    end

end

SEM_high_trial = std(avg_hi_trialvar,[],2)./sqrt(nTrials); 
SEM_high_trial = SEM_high_trial';
SEM_low_trial = std(avg_lo_trialvar,[],2)./sqrt(nTrials);
SEM_low_trial = SEM_low_trial';

%% Plots 

%%%%%%%%%% SEM across trials
% addpath errorbar_files/
% % f1 = figure;
% % figure(f1)
% figure()
% s_attn=shadedErrorBar(t,mean(avg_hi_trialvar,2),SEM_high_trial,'lineProps','r');
% hold on
% s_attn.mainLine.LineWidth = 1.5;
% s_dist=shadedErrorBar(t,mean(avg_lo_trialvar,2),SEM_low_trial,'lineProps','b');
% s_dist.mainLine.LineWidth = 1.5;
% xline(hi_tones(2),'r--', 'High Onset','LabelVerticalAlignment','bottom', 'LineWidth',1.5,'FontWeight','bold');
% xline(hi_tones(1,[1 3]),'r--','LineWidth',1.5,'FontWeight','bold');
% xline(lo_tones(2),'b--','Low Onset','LabelVerticalAlignment','bottom','LineWidth',1.5,'FontWeight','bold');
% xline(lo_tones(1,[1 3]),'b--','LineWidth',1.5,'FontWeight','bold');
% xline(tones(1,7),'-',{'End of', 'Sequence'},'LabelVerticalAlignment','bottom','FontWeight','bold','LineWidth',1.5)
% legend({'Attn Hi', 'Attn Lo'},'FontSize',12)
% hold off
% ax = gca;
% ax.FontSize = 14;
% ax.FontWeight="bold";
% ylabel('Amplitude (microV)','FontSize',12,'FontWeight','bold')
% xlabel('Time (ms)','FontSize',12,'FontWeight','bold')
% title(['SASA ' num2str(subject_id) ': Attending Hi & Lo (SEM across trials)'])
% % title('Worst Performer')
% xlim([-100 2000])

%exportgraphics(f1,'Cortical.png','Resolution',1000)

%% Getting ITPC 
addpath Functions/
if subject_id == "P017"
    nTrials = 477;
    num_analyze_sweeps = nTrials;
    [fftmag_lo,itpc_lo,phase_lo] = ITPC_analysis_final(num_analyze_sweeps,ch_lo,nChan,nTrials,fs,subject_id);
else
    nTrials = 480;
    num_analyze_sweeps = nTrials;
    [fftmag_lo,itpc_lo,phase_lo] = ITPC_analysis_final(num_analyze_sweeps,ch_lo,nChan,nTrials,fs,subject_id);
end
[fftmag_hi,itpc_hi,phase_hi] = ITPC_analysis_final(num_analyze_sweeps,ch_hi,nChan,nTrials,fs,subject_id);

fft_mag_difference = fftmag_hi - fftmag_lo;


avg_neural_Hi_phase_across_channels = mean(phase_hi,1);
avg_neural_Lo_phase_across_channels = mean(phase_lo,1);

% avg_neural_Hi_itpc_across_channels = mean(itpc_hi,1);
% avg_neural_Lo_itpc_across_channels = mean(itpc_lo,1);


phase_both_cond = [avg_neural_Hi_phase_across_channels(:,5) (avg_neural_Lo_phase_across_channels(:,5)+pi)];
avg_phase_both_cond = mean(phase_both_cond); 

phase_shift = avg_neural_Hi_phase_across_channels(:,5) - avg_neural_Lo_phase_across_channels(:,5);

%% This script analyzes left vs right ear for both and across conditions
function [LEar_LDist_data,LEar_LAttn_data,REar_RDist_data,REar_RAttn_data,avg_LEar_combo_data,avg_REar_combo_data,t_LEar_combo,avg_LEar_LDist_data,avg_LEar_LAttn_data,avg_REar_RDist_data,avg_REar_RAttn_data] = abr_ear_analysis_function(subject_id,blip)

% clear all
% close all

subject_id = 'P095';
blip=1;

% Loading in Ear data
if blip == 1
    LEar_combo = load(['/Users/victoriafigarola/Documents/1_CMU/Barb_Lab/Projects/SASA/Data/SASA_' num2str(subject_id) '/ABR_LeftEarCombined.mat']);
    REar_combo = load(['/Users/victoriafigarola/Documents/1_CMU/Barb_Lab/Projects/SASA/Data/SASA_' num2str(subject_id) '/ABR_RightEarCombined.mat']);

    LEar_combo_data = squeeze(LEar_combo.EEG.data(1,:,:)); %43787
    REar_combo_data = squeeze(REar_combo.EEG.data(1,:,:)); %43596
    
    % if size(LEar_combo_data,1) < 328
    %     LEar_combo_data(82:328,:) = 0;
    % end
    % 
    % if size(REar_combo_data,1) < 328
    %     REar_combo_data(82:328,:) = 0;
    % end

    if size(LEar_combo_data,2) >= 32576
        LEar_combo_data = LEar_combo_data(:,1:32576);
    end

    if size(REar_combo_data,2) >= 54720
        REar_combo_data = REar_combo_data(:,1:54720);
    end

    %%%%%%%% L vs R ear
    LEar_LDist = load(['/Users/victoriafigarola/Documents/1_CMU/Barb_Lab/Projects/SASA/Data/SASA_' num2str(subject_id) '/ABR_LeftEarLowDist.mat']);
    LEar_LAttn = load(['/Users/victoriafigarola/Documents/1_CMU/Barb_Lab/Projects/SASA/Data/SASA_' num2str(subject_id) '/ABR_LeftEarLowAttn.mat']);
    REar_RDist = load(['/Users/victoriafigarola/Documents/1_CMU/Barb_Lab/Projects/SASA/Data/SASA_' num2str(subject_id) '/ABR_RightEarHighDist.mat']);
    REar_RAttn = load(['/Users/victoriafigarola/Documents/1_CMU/Barb_Lab/Projects/SASA/Data/SASA_' num2str(subject_id) '/ABR_RightEarHighAttn.mat']);
   
    LEar_LDist_data = squeeze(LEar_LDist.EEG.data(1,:,:)); %16288
    LEar_LAttn_data = squeeze(LEar_LAttn.EEG.data(1,:,:)); %

    
    if size(LEar_LDist_data,2) >= 16288 
        LEar_LDist_data = LEar_LDist_data(:,1:16288);
    end
    if size(LEar_LAttn_data,2) >= 16288
        LEar_LAttn_data = LEar_LAttn_data(:,1:16288);
    end

    % if size(LEar_LAttn_data,1) < 328
    %     LEar_LAttn_data(82:328,:) = 0;
    % end
    % 
    % if size(LEar_LDist_data,1) < 328
    %     LEar_LDist_data(82:328,:) = 0;
    % end

    REar_RDist_data = squeeze(REar_RDist.EEG.data(1,:,:)); %27308
    REar_RAttn_data = squeeze(REar_RAttn.EEG.data(1,:,:)); %27499

    if size(REar_RDist_data,2) >= 27304
        REar_RDist_data = REar_RDist_data(:,1:27304);
    end
    if size(REar_RAttn_data,2) >= 27391
        REar_RAttn_data = REar_RAttn_data(:,1:27391);
    end

    if size(REar_RDist_data,1) < 328
        REar_RDist_data(82:328,:) = 0;
    end

    if size(REar_RAttn_data,1) < 328
        REar_RAttn_data(82:328,:) = 0;
    end


elseif blip == 2
    LEar_combo = load(['/Users/victoriafigarola/Documents/1_CMU/Barb_Lab/Projects/SASA/Data/ABR_redoing_BPF/' num2str(subject_id) '-ABR_LeftEarCombined.mat']);
    REar_combo = load(['/Users/victoriafigarola/Documents/1_CMU/Barb_Lab/Projects/SASA/Data/ABR_redoing_BPF/' num2str(subject_id) '-ABR_RightEarCombined.mat']);

    LEar_combo_data = squeeze(LEar_combo.EEG.data(1,:,:)); %43787
    REar_combo_data = squeeze(REar_combo.EEG.data(1,:,:)); %43596


    if size(LEar_combo_data,2) >= 32576
        LEar_combo_data = LEar_combo_data(:,1:32576);
    end

    if size(REar_combo_data,2) >= 54709
        REar_combo_data = REar_combo_data(:,1:54709);
    end

    %%%%%%%% L vs R ear
    LEar_LDist = load(['/Users/victoriafigarola/Documents/1_CMU/Barb_Lab/Projects/SASA/Data/ABR_redoing_BPF/' num2str(subject_id) '-ABR_LeftEarLowDist.mat']);
    LEar_LAttn = load(['/Users/victoriafigarola/Documents/1_CMU/Barb_Lab/Projects/SASA/Data/ABR_redoing_BPF/' num2str(subject_id) '-ABR_LeftEarLowAttn.mat']);
    REar_RDist = load(['/Users/victoriafigarola/Documents/1_CMU/Barb_Lab/Projects/SASA/Data/ABR_redoing_BPF/' num2str(subject_id) '-ABR_RightEarHighDist.mat']);
    REar_RAttn = load(['/Users/victoriafigarola/Documents/1_CMU/Barb_Lab/Projects/SASA/Data/ABR_redoing_BPF/' num2str(subject_id) '-ABR_RightEarHighAttn.mat']);

    LEar_LDist_data = squeeze(LEar_LDist.EEG.data(1,:,:)); %
    LEar_LAttn_data = squeeze(LEar_LAttn.EEG.data(1,:,:)); %

    if size(LEar_LDist_data,2) >= 16288
        LEar_LDist_data = LEar_LDist_data(:,1:16288);
    end
    if size(LEar_LAttn_data,2) >= 16288
        LEar_LAttn_data = LEar_LAttn_data(:,1:16288);
    end

    REar_RDist_data = squeeze(REar_RDist.EEG.data(1,:,:)); %
    REar_RAttn_data = squeeze(REar_RAttn.EEG.data(1,:,:)); %

    if size(REar_RDist_data,2) >= 27308 
        REar_RDist_data = REar_RDist_data(:,1:27308);
    end
    if size(REar_RAttn_data,2) >= 27391
        REar_RAttn_data = REar_RAttn_data(:,1:27391);
    end


end


t_LEar_combo = LEar_combo.EEG.times; %in ms
t_REar_combo = REar_combo.EEG.times; %in ms

LEar_combo_data = squeeze(LEar_combo.EEG.data(1,:,:)); %43787
REar_combo_data = squeeze(REar_combo.EEG.data(1,:,:)); %43596


%Indexing window of wave V (between 5 and 6.5ms, defined at top of script)
% t_indx = t_REar_combo(:,start_waveV:end_waveV);
% LEar_combo_data_waveV = avg_LEar_combo_data(:,start_waveV:end_waveV);
% REar_combo_data_waveV = avg_REar_combo_data(:,start_waveV:end_waveV);


% avg_LEar_combo_waveV = mean(LEar_combo_data_waveV);
% avg_REar_combo_waveV = mean(REar_combo_data_waveV);


avg_LEar_LDist_data = mean(LEar_LDist_data,2)';
avg_LEar_LAttn_data = mean(LEar_LAttn_data,2)';

avg_REar_RDist_data = mean(REar_RDist_data,2)';
avg_REar_RAttn_data = mean(REar_RAttn_data,2)';


SEM_LEar_LDist_data = std(LEar_LDist_data,[],2)./sqrt(size(LEar_LDist_data,2));
SEM_LEar_LAttn_data = std(LEar_LAttn_data,[],2)./sqrt(size(LEar_LAttn_data,2));

SEM_REar_RDist_data = std(REar_RDist_data,[],2)./sqrt(size(REar_RDist_data,2));
SEM_REar_RAttn_data = std(REar_RAttn_data,[],2)./sqrt(size(REar_RAttn_data,2));

%%

% addpath errorbar_files/
% figure()
% subplot(1,2,1)
% sgtitle(['SASA ' num2str(subject_id)])
% s_L_attn=shadedErrorBar(t_REar_combo,avg_LEar_LAttn_data,SEM_LEar_LAttn_data,'lineProps','r');
% hold on
% s_L_dist=shadedErrorBar(t_REar_combo,avg_LEar_LDist_data,SEM_LEar_LDist_data,'lineProps','b');
% s_L_dist.mainLine.LineWidth = 3;
% s_L_attn.mainLine.LineWidth = 3;
% xlim([-0.6 14])
% % ylim([-0.1 0.2])
% legend({'Attn', 'Dist'},'FontSize',14)
% ylabel('Amplitude (microV)','FontSize',14,'FontWeight','bold')
% % xlabel('Time (ms)','FontSize',14,'FontWeight','bold')
% set(gca,'FontSize',14,'FontWeight','bold');
% title('Left Ear: Attn vs. Dist')
% 
% subplot(1,2,2)
% s_R_attn=shadedErrorBar(t_REar_combo,avg_REar_RAttn_data,SEM_REar_RAttn_data,'lineProps','r');
% hold on
% s_R_dist=shadedErrorBar(t_REar_combo,avg_REar_RDist_data,SEM_REar_RDist_data,'lineProps','b');
% s_R_dist.mainLine.LineWidth = 3;
% s_R_attn.mainLine.LineWidth = 3;
% xlim([-0.6 14])
% % ylim([-0.1 0.2])
% % legend({'Attn', 'Dist'},'FontSize',14)
% ylabel('Amplitude (microV)','FontSize',14,'FontWeight','bold')
% xlabel('Time (ms)','FontSize',14,'FontWeight','bold')
% set(gca,'FontSize',14,'FontWeight','bold');
% title('Right Ear: Dist vs Attn')



function [attn_data, avg_attn_data, dist_data, avg_dist_data, t, fs] = abr_analysis_function(subject_id,blip)
% clear all

% subject_id = 'P088';
% blip = 1 ;

% Loading in Attn vs Dist data
if blip == 1
    attn = load(['/Users/victoriafigarola/Documents/1_CMU/Barb_Lab/Projects/SASA/Data/SASA_' num2str(subject_id) '/Attn_abr100.mat']);
    dist = load(['/Users/victoriafigarola/Documents/1_CMU/Barb_Lab/Projects/SASA/Data/SASA_' num2str(subject_id) '/Dist_abr100.mat']);

    attn_data = squeeze(attn.EEG.data(1,:,:)); %43787; with blip =2, 43596
    dist_data = squeeze(dist.EEG.data(1,:,:)); %43596


    if length(attn_data) >= 43628
        if subject_id == 'P003' | subject_id == 'P006' | subject_id == 'P007' | subject_id == 'P009' | subject_id == 'P012' | subject_id == 'P013' | subject_id == 'P016' | ...
                subject_id == 'P017' | subject_id == 'P019' | subject_id == 'P020'
            attn_data = attn_data(:,1:43785);

        else %subject_id == 'P024' | subject_id == 'P027' | subject_id == 'P028' | subject_id == 'P029' | subject_id == 'P030' | subject_id == 'P036' | subject_id == 'P037' | ...
            % subject_id == 'P039' | subject_id == 'P040'| subject_id == 'P044'| subject_id == 'P046'
            attn_data(:,43629:43785) = 0;
        end
    end

    if size(attn_data,1) < 328
            attn_data(82:328,:) = 0;
    end

    if size(dist_data,1) < 328
        dist_data(82:328,:) = 0;
    end


    if length(dist_data) > 43592
        dist_data = dist_data(:,1:43592);
    end


elseif blip == 2
    attn = load(['/Users/victoriafigarola/Documents/1_CMU/Barb_Lab/Projects/SASA/Data/ABR_redoing_BPF/' num2str(subject_id) '-Attn_abr30.mat']);
    dist = load(['/Users/victoriafigarola/Documents/1_CMU/Barb_Lab/Projects/SASA/Data/ABR_redoing_BPF/' num2str(subject_id) '-Dist_abr30.mat']);

    attn_data = squeeze(attn.EEG.data(1,:,:)); %43787; with blip =2, 43596
    dist_data = squeeze(dist.EEG.data(1,:,:)); %43596

    if length(attn_data) >= 43596
        attn_data = attn_data(:,1:43596);
    end
    if length(dist_data) >= 43596
        dist_data =dist_data(:,1:43596);
    end


end

t = attn.EEG.times; %in ms
fs = attn.EEG.srate; %samples/sec

avg_attn_data = mean(attn_data,2)';
% std_attn_data = std(attn_data,[],2)';
avg_dist_data = mean(dist_data,2)';
% std_dist_data = std(dist_data,[],2)';
SEM_attn = std(attn_data,[],2)./sqrt(size(attn_data,2));
SEM_dist = std(dist_data,[],2)./sqrt(size(dist_data,2));


% %%%%%%%%%%% Plotting Both Attn Hi & Attn Lo
% 
% addpath errorbar_files/
% f1 = figure;
% figure(f1)
% hold on;
% shadedErrorBar(t,avg_attn_data,SEM_attn,'lineProps','r')
% shadedErrorBar(t,avg_dist_data,SEM_dist,'lineProps','b')
% legend({'Attending', 'Ignoring Distractor'},'FontSize',14)
% set(gca,'FontSize',14,'FontWeight','bold');
% xlim([-0.5 15])
% hold off
% ylabel('Amplitude (microV)','FontSize',14,'FontWeight','bold')
% xlabel('Time (ms)','FontSize',14,'FontWeight','bold')
% title(['SASA ' num2str(subject_id) ' ABR: Attending vs Ignoring Distractor'])

% ---- Save a high-res version of the picture ----
% exportgraphics(f1,'Cortical.png','Resolution',1000)





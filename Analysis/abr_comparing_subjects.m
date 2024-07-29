% Victoria Figarola, 1/30/23
% This script collapses data from all subjects to find the average ABR
% figure, and average peak and latency for wave V for both conditions
% across all subjects. This script also collapses data from all subjects to
% plot the Right vs. Left ear across conditions. 
% Edit the following:
    %nSubj
    %subj_ID

nSubj = 44;
subj_ID = ["P003","P006", "P007","P009","P012","P013","P016",...
    "P017","P019","P020","P024","P027","P028","P029","P030",...
    "P036","P037","P039","P040","P044","P046","P048","P049",...
    "P051","P053","P054","P056","P062","P064","P065","P066",...
    "P067","P068","P076","P077","P078","P079","P083","P085",...
    "P086","P090","P092","P094","P095"];

% subj_ID_p088 = "P088";
% Getting the data of attn vs dist
for i = 1:nSubj
    [attn_data(i,:,:), avg_attn_data(i,:), dist_data(i,:,:), avg_dist_data(i,:), t, fs] = abr_analysis_function(subj_ID(i),1);
        %attn_data = 27x328x43785 (# of subj x samples x epochs)
        %avg_attn_data = average across all epochs; 27x328 (# of subj x samples)
    [LEar_LDist_data(i,:,:),LEar_LAttn_data(i,:,:),REar_RDist_data(i,:,:),REar_RAttn_data(i,:,:),avg_ABR_LEar(i,:),avg_ABR_REar(i,:),t_ABR,avg_ABR_LEar_LDist(i,:),avg_ABR_LEar_LAttn(i,:),avg_ABR_REar_RDist(i,:),avg_ABR_REar_RAttn(i,:)] = abr_ear_analysis_function(subj_ID(i),1);
    if i == 34
        [attn_data(i,:,:), avg_attn_data(i,:), dist_data(i,:,:), avg_dist_data(i,:), ~, fs] = abr_analysis_function(subj_ID(i),1);
        [LEar_LDist_data(i,:,:),LEar_LAttn_data(i,:,:),REar_RDist_data(i,:,:),REar_RAttn_data(i,:,:),avg_ABR_LEar(i,:),avg_ABR_REar(i,:),~,avg_ABR_LEar_LDist(i,:),avg_ABR_LEar_LAttn(i,:),avg_ABR_REar_RDist(i,:),avg_ABR_REar_RAttn(i,:)] = abr_ear_analysis_function(subj_ID(i),1);

    end

end 


attending_avg = mean(avg_attn_data);
distractor_avg = mean(avg_dist_data);

% attn_peaks_avg = mean(attn_pks);
% attn_peaks_latency_avg = mean(attn_loc_time);

% dist_peaks_avg = mean(dist_pks);
% dist_peaks_latency_avg = mean(dist_loc_time);

%% Average ABR Attn vs Dist -- RUN THIS SECTION
% addpath errorbar_files/
% f1 = figure;
% figure(f1)
% hold on;
% s_attn=shadedErrorBar(t,mean(avg_attn_data,1),std(avg_attn_data)/sqrt(nSubj),'lineProps','k');
% s_dist=shadedErrorBar(t,mean(avg_dist_data,1),std(avg_dist_data)/sqrt(nSubj),'lineProps','g');
% s_attn.mainLine.LineWidth = 3;
% s_dist.mainLine.LineWidth = 3;
% xline(5.92041,'--r','Avg. Attn. Latency','LineWidth',1.5,'LabelHorizontalAlignment','left', ...
%     'LabelVerticalAlignment','bottom','FontSize',12,'FontWeight','bold')
% xline(5.85938,'--b','Avg. Dist Latency','LineWidth',1.5,'LabelHorizontalAlignment','right', ...
%     'LabelVerticalAlignment','bottom','FontSize',12,'FontWeight','bold')
% yline(0.113514,'--r','Avg. Attn. Peak','LineWidth',1.5,'LabelHorizontalAlignment','left', ...
%     'LabelVerticalAlignment','bottom','FontSize',12,'FontWeight','bold')
% yline(0.111705,'--b','Avg. Dist. Peak','LineWidth',1.5,'LabelHorizontalAlignment','left', ...
%     'LabelVerticalAlignment','top','FontSize',12,'FontWeight','bold')
% xlim([-0.6 15])
% legend({'Attending', 'Ignoring Distractor'},'FontSize',14)
% set(gca,'FontSize',14,'FontWeight','bold');
% hold off
% ylabel('Amplitude (microV)','FontSize',14,'FontWeight','bold')
% xlabel('Time (ms)','FontSize',14,'FontWeight','bold')
% title(['ABR: Attending vs Ignoring Distractor (n=' num2str(nSubj) ')'])

% exportgraphics(f1,'Subcortical.png','Resolution',1000)

% figure
% plot(t,mean(avg_attn_data,1),'r');
% hold on
% plot(t,mean(avg_dist_data,1),'b')
% legend('attn','dist')

%% Plotting Right vs. Left Ears --- RUN THIS SECITON
% addpath errorbar_files/
% f2 = figure;
% figure(f2)
% hold on;
% shadedErrorBar(t_ABR,mean(avg_ABR_LEar),std(avg_ABR_LEar)/sqrt(nSubj),'lineProps','r')
% shadedErrorBar(t_ABR,mean(avg_ABR_REar),std(avg_ABR_REar)/sqrt(nSubj),'lineProps','b')
% legend({'Left Ear', 'Right Ear'},'FontSize',14)
% set(gca,'FontSize',14,'FontWeight','bold');
% hold off
% ylabel('Amplitude (microV)','FontSize',14,'FontWeight','bold')
% xlabel('Time (ms)','FontSize',14,'FontWeight','bold')
% title(['ABR: Left Ear vs Right Ear (n=' num2str(nSubj) ')'])

addpath errorbar_files/
f3 = figure;
figure(f3)
subplot(1,2,1)
% sgtitle(['ABR: Right & Left Ear with conditions (n=' num2str(nSubj) ')'])
hold on 
s_L_attn=shadedErrorBar(t_ABR,mean(avg_ABR_LEar_LAttn),std(avg_ABR_LEar_LAttn)/sqrt(nSubj),'lineProps','r');
s_L_dist = shadedErrorBar(t_ABR,mean(avg_ABR_LEar_LDist),std(avg_ABR_LEar_LDist)/sqrt(nSubj),'lineProps','b');

s_L_attn.mainLine.LineWidth = 3;
s_L_dist.mainLine.LineWidth = 3;
xlim([-0.6 15])
ylim([-0.1 0.2])
legend({'Attending','Ignoring Distractor'},'FontSize',12)
ylabel('Amplitude (microV)','FontSize',14,'FontWeight','bold')
xlabel('Time (ms)','FontSize',14,'FontWeight','bold')
set(gca,'FontSize',14,'FontWeight','bold');
title('Left Ear (Attn/Ignore Low)')

subplot(1,2,2)
% figure()
s_R_attn=shadedErrorBar(t_ABR,mean(avg_ABR_REar_RAttn),std(avg_ABR_REar_RAttn)/sqrt(nSubj),'lineProps','r');
hold on
s_R_dist=shadedErrorBar(t_ABR,mean(avg_ABR_REar_RDist),std(avg_ABR_REar_RDist)/sqrt(nSubj),'lineProps','b');
s_R_attn.mainLine.LineWidth = 3;
s_R_dist.mainLine.LineWidth = 3;
xlim([-0.5 15])
ylim([-0.1 0.2])
% legend({'Attending','Ignoring Distractor'},'FontSize',14)
ylabel('Amplitude (microV)','FontSize',14,'FontWeight','bold')
xlabel('Time (ms)','FontSize',14,'FontWeight','bold')
set(gca,'FontSize',14,'FontWeight','bold');
title('Right Ear (Attn/Ignore High)')



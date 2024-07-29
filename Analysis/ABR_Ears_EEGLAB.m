% This script calls in pre-processed ABR data in order to epoch by Left &
% Right Ears, as well as Attn vs Dist for both ears. It is saved to
% corresponding subject folder
% ------------------------------------------------
% Change lines 7 and 10 only 

subject_id = 'P095';

% [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
% EEG = pop_loadset('filename','p058_abr_fixedtrigs.set','filepath','/Users/victoriafigarola/Library/CloudStorage/Box-Box/Frog-quency_Pt2/Datasets/SASA_P058/ABR/');
% [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
% EEG = eeg_checkset( EEG );


% [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
%     % EEG = pop_loadset('filename','P071_ABR_100BPF.set','filepath','/Users/victoriafigarola/Library/CloudStorage/Box-Box/Frog-quency_Pt2/Datasets/SASA_P071/ABR/');
% EEG = pop_loadset('filename','p086_abr_fixedtrigs.set','filepath','/Users/victoriafigarola/Library/CloudStorage/Box-Box/Frog-quency_Pt2/Datasets/SASA_P086/ABR/');
%     % EEG = pop_loadset('filename','p007_abr30hz_fixedtrigs.set','filepath','/Users/victoriafigarola/Library/CloudStorage/Box-Box/Frog-quency_Pt2/ABR_redo_filtering/');
% [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
% EEG = eeg_checkset( EEG );


% Left Ear (attn + dist)
EEG = pop_epoch( EEG, {  '1'  '9'  }, [-0.005       0.015], 'newname', 'Left Ear Combined', 'epochinfo', 'yes');
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'gui','off'); 
EEG = eeg_checkset( EEG );
EEG = pop_rmbase( EEG, [],[]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 2,'gui','off'); 
EEG = eeg_checkset( EEG );
save(['/Users/victoriafigarola/Documents/1_CMU/Barb_Lab/Projects/SASA/Data/SASA_' subject_id '/ABR_LeftEarCombined.mat'],"EEG")
% save(['/Users/victoriafigarola/Documents/1_CMU/Barb_Lab/Projects/SASA/Data/ABR_redoing_BPF/' subject_id '-ABR_LeftEarCombined.mat'],"EEG")

% Right ear (attn + dist)
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 3,'retrieve',1,'study',0); 
EEG = eeg_checkset( EEG );
EEG = pop_epoch( EEG, {  '7'  '15'  }, [-0.005       0.015], 'newname', 'Right Ear Combined', 'epochinfo', 'yes');
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'gui','off'); 
EEG = eeg_checkset( EEG );
EEG = pop_rmbase( EEG, [],[]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 4,'gui','off'); 
save(['/Users/victoriafigarola/Documents/1_CMU/Barb_Lab/Projects/SASA/Data/SASA_' subject_id '/ABR_RightEarCombined.mat'],"EEG")
% save(['/Users/victoriafigarola/Documents/1_CMU/Barb_Lab/Projects/SASA/Data/ABR_redoing_BPF/' subject_id '-ABR_RightEarCombined.mat'],"EEG")

% Left ear, Low Dist
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 5,'retrieve',1,'study',0); 
EEG = eeg_checkset( EEG );
EEG = pop_epoch( EEG, {  '1'  }, [-0.005       0.015], 'newname', 'Left Ear Low Dist', 'epochinfo', 'yes');
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'gui','off'); 
EEG = eeg_checkset( EEG );
EEG = pop_rmbase( EEG, [],[]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 6,'gui','off'); 
save(['/Users/victoriafigarola/Documents/1_CMU/Barb_Lab/Projects/SASA/Data/SASA_' subject_id '/ABR_LeftEarLowDist.mat'],"EEG")
% save(['/Users/victoriafigarola/Documents/1_CMU/Barb_Lab/Projects/SASA/Data/ABR_redoing_BPF/' subject_id '-ABR_LeftEarLowDist.mat'],"EEG")


% Left Ear, Low Attn
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 7,'retrieve',1,'study',0); 
EEG = eeg_checkset( EEG );
EEG = pop_epoch( EEG, {  '9'  }, [-0.005       0.015], 'newname', 'Left ear Low attn', 'epochinfo', 'yes');
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'gui','off'); 
EEG = eeg_checkset( EEG );
EEG = pop_rmbase( EEG, [],[]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 8,'gui','off'); 
save(['/Users/victoriafigarola/Documents/1_CMU/Barb_Lab/Projects/SASA/Data/SASA_' subject_id '/ABR_LeftEarLowAttn.mat'],"EEG")
% save(['/Users/victoriafigarola/Documents/1_CMU/Barb_Lab/Projects/SASA/Data/ABR_redoing_BPF/' subject_id '-ABR_LeftEarLowAttn.mat'],"EEG")


% Right Ear, High Attn
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 9,'retrieve',1,'study',0); 
EEG = eeg_checkset( EEG );
EEG = pop_epoch( EEG, {  '15'  }, [-0.005       0.015], 'newname', 'Right Ear High Attn', 'epochinfo', 'yes');
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'gui','off'); 
EEG = eeg_checkset( EEG );
EEG = pop_rmbase( EEG, [],[]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 10,'gui','off'); 
save(['/Users/victoriafigarola/Documents/1_CMU/Barb_Lab/Projects/SASA/Data/SASA_' subject_id '/ABR_RightEarHighAttn.mat'],"EEG")
% save(['/Users/victoriafigarola/Documents/1_CMU/Barb_Lab/Projects/SASA/Data/ABR_redoing_BPF/' subject_id '-ABR_RightEarHighAttn.mat'],"EEG")


% Right Ear, Low Attn 
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 11,'retrieve',1,'study',0); 
EEG = eeg_checkset( EEG );
EEG = pop_epoch( EEG, {  '7'  }, [-0.005       0.015], 'newname', 'R Ear H Dist', 'epochinfo', 'yes');
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'gui','off'); 
EEG = eeg_checkset( EEG );
EEG = pop_rmbase( EEG, [],[]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 12,'gui','off'); 
save(['/Users/victoriafigarola/Documents/1_CMU/Barb_Lab/Projects/SASA/Data/SASA_' subject_id '/ABR_RightEarHighDist.mat'],"EEG")
% save(['/Users/victoriafigarola/Documents/1_CMU/Barb_Lab/Projects/SASA/Data/ABR_redoing_BPF/' subject_id '-ABR_RightEarHighDist.mat'],"EEG")


% Dist -- both ears
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 13,'retrieve',1,'study',0); 
EEG = eeg_checkset( EEG );
EEG = pop_epoch( EEG, {  '1'  '7'  }, [-0.005       0.015], 'newname', 'Dist', 'epochinfo', 'yes');
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'gui','off'); 
EEG = eeg_checkset( EEG );
EEG = pop_rmbase( EEG, [],[]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 14,'gui','off'); 
EEG = eeg_checkset( EEG );
save(['/Users/victoriafigarola/Documents/1_CMU/Barb_Lab/Projects/SASA/Data/SASA_' subject_id '/Dist_abr100.mat'],"EEG")
% save(['/Users/victoriafigarola/Documents/1_CMU/Barb_Lab/Projects/SASA/Data/ABR_redoing_BPF/' subject_id '-Dist_abr30.mat'],"EEG")


% Attn -- both ears
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 15,'retrieve',1,'study',0); 
EEG = pop_epoch( EEG, {  '9'  '7'  }, [-0.005       0.015], 'newname', 'Dist', 'epochinfo', 'yes');
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'gui','off'); 
EEG = eeg_checkset( EEG );
EEG = pop_rmbase( EEG, [],[]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 16,'gui','off'); 
EEG = eeg_checkset( EEG );
save(['/Users/victoriafigarola/Documents/1_CMU/Barb_Lab/Projects/SASA/Data/SASA_' subject_id '/Attn_abr100.mat'],"EEG")
% save(['/Users/victoriafigarola/Documents/1_CMU/Barb_Lab/Projects/SASA/Data/ABR_redoing_BPF/' subject_id '-Attn_abr30.mat'],"EEG")




eeglab redraw;

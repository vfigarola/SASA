% Wrapper script to make a full SASA block
% Requires sasaWithPips function

%% SETTINGS
addpath("combinator_update/combinator/")
basedir = '/Users/victoriafigarola/Documents/1_CMU/Barb_Lab/Projects/SASA/Code/StimGen/Stim_1Hz_within_band_sample/';
basedir_freq = '/Users/victoriafigarola/Documents/1_CMU/Barb_Lab/Projects/SASA/Code/StimGen/Stim_1Hz_within_band_sample/freq_order/';
basedir_trig = '/Users/victoriafigarola/Documents/1_CMU/Barb_Lab/Projects/SASA/Code/StimGen/Stim_1Hz_within_band_sample/'; 
file_prefix = 'attendLo';

% block can begin with single-band stimuli to orient listener
nr_singleBandTrials = 2;
nr_singleBandRepeats = 1;

nrBlocks = 1; % how many blocks per block?
fs = 48000; % sampling rate
samples_per_tone = 10800; %4 Hz rate
samples_after_tone = 1200; % how long should the gap between tones be (in samples)?
samples_after_sequence = 24000;  % how long should the silent gap after each sequence be (in samples)?  
    %length of 2 tones = 500 ms

%%%%%%%% Making it 1 and 2Hz rate
% nrBlocks = 1; % how many blocks per block?
% fs = 48000; % sampling rate
% samples_per_tone = ceil(0.34*fs); %1 Hz rate within band (2Hz across band)
% samples_after_tone = 0.01*fs; % how long should the gap between tones be (in samples)?
% samples_after_sequence = 0.7*fs;  % how long should the silent gap after each sequence be (in samples)?  
    %length of 2 tones = 666 ms


%%%%%%%%% What I'm actually using
lo_freq_band = [43.65, 49, 55]; %F1 G1 A1
hi_freq_band = [73.42, 82.41, 92.50]; %D2, E2, F#2

target_band = 'lo'; %'hi' or 'lo'
seq_length = 3; % number of tones per sequence

% trials_per_block = 30; % how many trials are there per block?
trials_per_block = 20; % how many trials are there per block?

nrRepTrials_min = 3; % how many trials should have pattern repeats in the target band? there'll be an equal number of trials with pattern repeats in the distractor band
nrRepTrials_max = 6; % how many trials should have pattern repeats in the target band? there'll be an equal number of trials with pattern repeats in the distractor band

interleaved = 2; % the target and distractor are interleaved and target comes first (1), or target and distractor are interleaved and target comes second (2) in sequence



%% Make a block of blocks

% We begin the block by presenting one trial with *just* the target band

[singleband, ~, block_seqs_freq, dist_block_seqs_freq, ~, ~, ~,~] = sasaWithPips_new(fs, samples_per_tone, samples_after_sequence, samples_after_tone, lo_freq_band, hi_freq_band, target_band, seq_length, nr_singleBandTrials, nr_singleBandRepeats, nr_singleBandRepeats, interleaved, 0);
% dur_singleband = length(singleband)/fs;

% And get the length of the block so far
block_length = size(singleband,2)/fs;

% Create a file where we'll save all the repeat times for the block
block_repeat_times = [basedir file_prefix '-block_repeat_times_1:24.txt'];


% Now, we can start adding each block to this
for i = 1:nrBlocks

    % Generate the block
    [sasa_block_waveform, dist_sasa_block_waveform, block_seqs_freq, dist_block_seqs_freq, rep_pos_times, pip_times_attn, pip_times_dist,impulse_train] = sasaWithPips_new(fs, samples_per_tone, samples_after_sequence, samples_after_tone, lo_freq_band, hi_freq_band, target_band, seq_length, trials_per_block, nrRepTrials_min, nrRepTrials_max, interleaved, 1);

    sasa_block_waveform = [singleband sasa_block_waveform];

    dist_sasa_block_waveform = [zeros(1, size(singleband,2)) dist_sasa_block_waveform];


    if strcmp(target_band,'hi')
        wavename = [basedir file_prefix, '-block',num2str(i),'-right.wav'];
        sasa = [dist_sasa_block_waveform(:),sasa_block_waveform(:)];
        % audiowrite(wavename,sasa,fs);
        
        blockinfo_attn = [basedir_freq file_prefix, '-block',num2str(i),'-attend-block-freqs.xlsx'];
        % writematrix(round(block_seqs_freq), blockinfo_attn);

        blockinfo_dist = [basedir_freq file_prefix, '-block',num2str(i),'-dist-block-freqs.xlsx'];
        % writematrix(round(dist_block_seqs_freq), blockinfo_dist);

    elseif ~strcmp(target_band,'hi')
        wavename = [basedir file_prefix, '-block',num2str(i),'-left.wav'];
        sasa = [sasa_block_waveform(:),dist_sasa_block_waveform(:)];
        % audiowrite(wavename,sasa,fs)

        blockinfo_attn = [basedir_freq file_prefix, '-block',num2str(i),'-attend-block-freqs.xlsx'];
        % writematrix(round(block_seqs_freq), blockinfo_attn);
        blockinfo_dist = [basedir_freq file_prefix, '-block',num2str(i),'-dist-block-freqs.xlsx'];
        % writematrix(round(dist_block_seqs_freq), blockinfo_dist);
    end

    %Setting up triggers
    triggers_attn = [zeros(1, size(singleband,2)) pip_times_attn];
    triggers_attn = triggers_attn * (1 + 2*(strcmp(target_band,'hi')) + 8*1); %1/2/8 means the bit, 8*1 (1= attended)
    wavname = [basedir_trig file_prefix, '-block', num2str(i), '-attn_trig.wav'];
    % audiowrite(wavname, triggers_attn/100, fs); %will need to multiply by 100 in experiment; dividing to avoid clipping

    triggers_dist = [zeros(1, size(singleband,2)) pip_times_dist];
    triggers_dist = triggers_dist * (1 + 2*(~strcmp(target_band,'hi')) + 8*0);
    wavname = [basedir_trig file_prefix, '-block', num2str(i), '-dist_trig.wav'];
    % audiowrite(wavname, triggers_dist/100, fs);

    %Getting target windows
        rep_pos_times = block_length*1000 + rep_pos_times; %in ms
        % Create a string with the repetition times for this block
        block_rep_times = mat2str(rep_pos_times);
        block_rep_times = strrep(block_rep_times, ' ', ',');
        block_rep_times = strrep(block_rep_times, ';', '],[');
        block_rep_times = ['[' block_rep_times ']' newline];
    
        % And then save the repetition times to a master text file for the block
        % (But if we're on the first block, create a new file)
        if i == 1
            fid = fopen(block_repeat_times, 'w');
            fprintf(fid, block_rep_times);
            fclose(fid);
        else
            fid = fopen(block_repeat_times, 'a+');
            fprintf(fid, block_rep_times);
            fclose(fid);
        end


end


% figure; plot(sasa(:,1)); hold on; plot(triggers_attn); plot(triggers_dist); legend('sasa','attn','dist')


%% GET FILE OF EXACT FREQ IN BAND for both attn and dist

% sound(sasa,fs)
% 
% t = 0:1/fs:(length(lo_tones)-1)/fs;
% figure
% plot(t,lo_tones(2,:),'b','LineWidth',2.5);
% xlabel('Time (sec)')
% ylabel('Amplitude (a.u.)')
% set(gca,'FontSize',12,'FontWeight','bold')


% loenv_sasa = envelope(sasa(:,1));
% hienv_sasa = envelope(sasa(:,2));


% basedir_training = '/Users/victoriafigarola/Documents/1_CMU/Barb_Lab/Code/Fred/Frog-quency/Stim-11:23/Dichotic/';
%     diotic = zeros(1,length(sasa_block_waveform));
%     training_sasa = [diotic(:),sasa_block_waveform(:)];
%     wavename = [basedir_training file_prefix, '-block', num2str(i),'.wav'];
%     audiowrite(wavename,training_sasa,fs)





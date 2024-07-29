function [sasa_block_waveform, dist_sasa_block_waveform, block_seqs_freq, dist_block_seqs_freq, rep_pos_times, pip_times_attn, pip_times_dist,impulse_train] = sasaWithPips_new(fs, samples_per_tone, samples_after_sequence, samples_after_tone, lo_freq_band, hi_freq_band, target_band, seq_length, trials_per_block, nrRepTrials_min, nrRepTrials_max, interleaved, distractorRepeats)

% SETUP -----
% Read in gap parameters
sil_gap = zeros(1,samples_after_sequence); %gap between each trial
tone_gap = zeros(1, samples_after_tone); %gap between each tone


% CREATE ALL TONES -----
% Parameters for hi/lo pips
mpip_lo = maddox_pip(3500,fs, 1)'; %calling in maddox pip script
mpip_hi = maddox_pip(4500,fs, 1)';


end_buffer = 0.005;
end_buffer_samples = end_buffer * fs;


% Create lo tones
for i = 1:size(lo_freq_band,2)
    Hz = lo_freq_band(i);
    samples_per_cycle = round(fs / Hz,0);
%     cycles_per_tone = round(samples_per_tone/samples_per_cycle,0);

    %   Create an impulse train
    impulse_train = zeros(1, samples_per_tone);
    impulse_train(1:samples_per_cycle:end-end_buffer_samples) = 1;

    %     Find all indexes of the impulse train and set every other index to -1
    idx = find(impulse_train);
    impulse_train(idx(1:2:end)) = -1;

    % Figure out if we'll need to add/subtract extra samples

    % Create tone
    tones = conv(impulse_train, mpip_lo); %convolve impulse train with tone pip

    % Trim
    lo_tones(i, :) = tones(1:samples_per_tone);

end

% t_train = 0:1/fs:(length(impulse_train)-1)/fs;
% figure()
% plot(t_train,impulse_train,'b','LineWidth',3);
% set(gca,'FontSize',12,'FontWeight','bold'); xlabel('Time (sec)'); ylabel('Amplitude')
% 
% t_pip = 0:1/fs:(length(lo_tones(1,:))-1)/fs;
% figure()
% plot(t_pip,lo_tones(1,:),'b','LineWidth',3);
% set(gca,'FontSize',12,'FontWeight','bold'); xlabel('Time (sec)'); ylabel('Amplitude')



% Create hi tones
for i = 1:size(hi_freq_band,2)
    Hz = hi_freq_band(i);
    samples_per_cycle = round(fs / Hz,0);
%     cycles_per_tone = round(samples_per_tone/samples_per_cycle,0);


     % Create an impulse train
    impulse_train = zeros(1, samples_per_tone); 
    impulse_train(1:samples_per_cycle:end-end_buffer_samples) = 1;

    % Find all indexes of the impulse train and set every other index to -1
    idx = find(impulse_train);
    impulse_train(idx(1:2:end)) = -1;

%     tones = [1, zeros(1,samples_per_cycle-size(mpip_lo,2))];
    tones = conv(impulse_train, mpip_hi);

    % Trim
    hi_tones(i, :) = tones(1:samples_per_tone);


end

% t_pip = 0:1/fs:(length(hi_tones(1,:))-1)/fs;
% figure()
% plot(t_pip,hi_tones(1,:),'r','LineWidth',3);
% set(gca,'FontSize',12,'FontWeight','bold'); xlabel('Time (sec)'); ylabel('Amplitude')


% CREATE MINISEQUENCES -----

% What are all the possible minisequences we could present?
% Suppose we have a 3-tone sequence, and there are four possible tones in
% the frequency band of interest. There are thus 4 possible tones for
% each position. First, we make a matrix of all the possible sequences.
% Each row is a different possible 3-tone sequence, and the numbers
% (1,2,...,4) indicate which particular tone from the band we're using.
% This uses the combinator script, which must be on the path
all_lo_seqs = combinator(size(lo_freq_band,2),seq_length,'p','r'); %set of possible digits x length (6 digits/seq) 
all_hi_seqs = combinator(size(hi_freq_band,2),seq_length,'p','r');

% all_lo_seqs = combinator(4,seq_length,'p','r');
% all_hi_seqs = combinator(4,seq_length,'p','r');


% Weed out minisequences where all notes are repeated (so that we don't
% present, e.g., 1-1-1 as a minisequence.)
% To do so, we use a "while" loop, where we test if all the elements in
% a row are the same, and if they are, we get rid of the row.

% We do this for the low-pitch sequences
i=1;
while i <= size(all_lo_seqs,1)
    if all_lo_seqs(i,1) == all_lo_seqs(i,:)
        all_lo_seqs(i,:) = [];
    end
    i = i+1;
end

% And then for the high-pitch sequences
i=1;
while i <= size(all_hi_seqs,1)
    if all_hi_seqs(i,1) == all_hi_seqs(i,:)
        all_hi_seqs(i,:) = [];
    end
    i = i+1;
end


% So let's say we are putting 30 sequences in a block.
% We pull 30 different row numbers from each matrix we made above.
% Each one corresponds to a different possible minisequence,
% so we don't currently have any repeats.

if target_band == 'hi'
    block_seqs = randsample(size(all_hi_seqs,1),trials_per_block,true);
    dist_block_seqs = randsample(size(all_lo_seqs,1),trials_per_block,true);
elseif target_band == 'lo'
    block_seqs = randsample(size(all_lo_seqs,1),trials_per_block,true);
    dist_block_seqs = randsample(size(all_hi_seqs,1),trials_per_block,true);
end

% while find(diff(block_seqs) == 0 | diff(dist_block_seqs) == 0)
while or(ismember(0,diff(block_seqs)),ismember(0,diff(dist_block_seqs)))
    if target_band == 'hi'
        block_seqs = randsample(size(all_hi_seqs,1),trials_per_block,true);
        dist_block_seqs = randsample(size(all_lo_seqs,1),trials_per_block,true);
    elseif target_band == 'lo'
        block_seqs = randsample(size(all_lo_seqs,1),trials_per_block,true);
        dist_block_seqs = randsample(size(all_hi_seqs,1),trials_per_block,true);
    end
end

block_seqs = block_seqs';
dist_block_seqs = dist_block_seqs';

% Now, depending on how many repeats we've specified (e.g., 4 minisequence
% repeats per block), we can choose which minisequences will get repeated.
% We'll space them out roughly evenly (e.g., first repeat in the first fourth
% of the block, second  in the second fourth, ..., and fourth in the
% final fourth).

% First, figure out how many repeats we need for this block

rep_range = nrRepTrials_max - nrRepTrials_min;

if rep_range == 0
    nrRepTrials = nrRepTrials_min;
else
    nrRepTrials = randperm(rep_range, 1) + nrRepTrials_min - 1;
end

% And then figure out which ones will be repeats
switch nrRepTrials

    % If we don't have any repeats, we can just go with what we have
    case 0
        disp('No repeats specified, moving on');

    otherwise

        % For the single-band stimuli at the start of a run, we may not
        % want to have any distractor stimuli
        if distractorRepeats == 0
            repeat_cnt = ones(1,nrRepTrials);
            reps = [repeat_cnt];
            distrep_pos = [];


            [r,targ_count] = size(reps);
            idx = randperm(targ_count);
            reps_rand = reps;
            reps_rand(1,idx) = reps(1,:);

            j = 1; % for rep_pos counter

            for i = 1:(targ_count)

                target_pos(i) = randperm(floor(trials_per_block/targ_count), 1)  + floor(trials_per_block/targ_count)*(i-1);

                % but we can't have a target in the first position, since we'll
                % be setting this element to be the same as the one before it
                while target_pos(i) == 1
                    target_pos(i) = randperm(floor(trials_per_block/targ_count), 1)  + floor(trials_per_block/targ_count)*(i-1);
                end

                if reps_rand(1,i) == 1
                    block_seqs(1, target_pos(i)) = block_seqs(1, target_pos(i)-1);

                    rep_pos(j) = target_pos(i);
                    j = j+1;

                end
            end

        else

            % Otherwise, we figure out where the repeats will be, doing so
            % separately for the target band and distractor band

            repeat_cnt = ones(1,nrRepTrials);
            distband_repeat_cnt = (ones(1, nrRepTrials))*2;

            reps = [repeat_cnt, distband_repeat_cnt];

            [r,targ_count] = size(reps);
            idx = randperm(targ_count);
            reps_rand = reps;
            reps_rand(1,idx) = reps(1,:);

            j = 1; % for rep_pos counter
            k = 1; % for distrep_pos  counter


            for i = 1:(targ_count)

                target_pos(i) = randperm(floor(trials_per_block/targ_count), 1)  + floor(trials_per_block/targ_count)*(i-1);

                % but we can't have a target in the first position, since we'll
                % be setting this element to be the same as the one before it
                while target_pos(i) == 1
                    target_pos(i) = randperm(floor(trials_per_block/targ_count), 1)  + floor(trials_per_block/targ_count)*(i-1);
                end

                if reps_rand(1,i) == 1
                    block_seqs(1,target_pos(i)) = block_seqs(1,target_pos(i)-1);

                    rep_pos(j) = target_pos(i);
                    j = j+1;

                else
                    if reps_rand(1,i) == 2
                        dist_block_seqs(1,target_pos(i)) = dist_block_seqs(1,target_pos(i)-1);
                        distrep_pos(k) = target_pos(i);
                        k = k+1;
                    end

                end
            end
        end

        % But if there are repeats in consecutive trials for either
        % frequency band, we have to redo this.

        while or(ismember(1,diff(rep_pos)),ismember(1,diff(distrep_pos)))

            % Starting with pulling new sequences

            if target_band == 'hi'
                block_seqs = randsample(size(all_hi_seqs,1),trials_per_block,true);
                dist_block_seqs = randsample(size(all_lo_seqs,1),trials_per_block,true);
            elseif target_band == 'lo'
                block_seqs = randsample(size(all_lo_seqs,1),trials_per_block,true);
                dist_block_seqs = randsample(size(all_hi_seqs,1),trials_per_block,true);
            end

            while or(ismember(0,diff(block_seqs)),ismember(0,diff(dist_block_seqs)))
                if target_band == 'hi'
                    block_seqs = randsample(size(all_hi_seqs,1),trials_per_block,true);
                    dist_block_seqs = randsample(size(all_lo_seqs,1),trials_per_block,true);
                elseif target_band == 'lo'
                    block_seqs = randsample(size(all_lo_seqs,1),trials_per_block,true);
                    dist_block_seqs = randsample(size(all_hi_seqs,1),trials_per_block,true);
                end
            end
            block_seqs = block_seqs';
            dist_block_seqs = dist_block_seqs';
            
            if distractorRepeats == 0
                repeat_cnt = ones(1,nrRepTrials);
                reps = [repeat_cnt];


                [r,targ_count] = size(reps);
                idx = randperm(targ_count);
                reps_rand = reps;
                reps_rand(1,idx) = reps(1,:);

                j = 1; % for rep_pos counter

                for i = 1:(targ_count)

                    target_pos(i) = randperm(floor(trials_per_block/targ_count), 1)  + floor(trials_per_block/targ_count)*(i-1);

                    % but we can't have a target in the first position, since we'll
                    % be setting this element to be the same as the one before it
                    while target_pos(i) == 1
                        target_pos(i) = randperm(floor(trials_per_block/targ_count), 1)  + floor(trials_per_block/targ_count)*(i-1);
                    end

                    if reps_rand(1,i) == 1
                        block_seqs(1, target_pos(i)) = block_seqs(1, target_pos(i)-1);

                        rep_pos(j) = target_pos(i);
                        j = j+1;

                    end
                end

            else

                repeat_cnt = ones(1,nrRepTrials);
                distband_repeat_cnt = (ones(1, nrRepTrials))*2;

                reps = [repeat_cnt, distband_repeat_cnt];

                [r,targ_count] = size(reps);
                idx = randperm(targ_count);
                reps_rand = reps;
                reps_rand(1,idx) = reps(1,:);

                j = 1; % for rep_pos counter
                k = 1; % for distrep_pos  counter


                for i = 1:(targ_count)

                    target_pos(i) = randperm(floor(trials_per_block/targ_count), 1)  + floor(trials_per_block/targ_count)*(i-1);

                    while target_pos(i) == 1
                        target_pos(i) = randperm(floor(trials_per_block/targ_count), 1)  + floor(trials_per_block/targ_count)*(i-1);
                    end

                    if reps_rand(1,i) == 1
                        block_seqs(1, target_pos(i)) = block_seqs(1, target_pos(i)-1);

                        rep_pos(j) = target_pos(i);
                        j = j+1;

                    else
                        if reps_rand(1,i) == 2
                            dist_block_seqs(1, target_pos(i)) = dist_block_seqs(1, target_pos(i)-1);
                            distrep_pos(k) = target_pos(i);
                            k = k+1;
                        end

                    end
                end
            end
        end

        % We now have our repeat positions, and we won't ever have two
        % consecutive trials where there are repeats in the same frequency
        % band.

        disp("Attended band repeats occur in trials:     " + int2str(rep_pos));
        disp("Distractor band repeats occur in trials:     " + int2str(distrep_pos));

end

% So now the block_seqs and dist_block_seqs matrices contain the repeats at
% the correct spots. Each cell has an index, which corresponds to a unique
% minisequence from the all_hi_seqs/all_lo_seqs matrices.

% Let's now create two matrices (block_seqs_freq and dist_block_seqs_freq)
% that specify the actual minisequences, with the actual frequency values.

% ATTENDED BAND
i=1;
if target_band == 'hi'
    block_seqs_freq = zeros(size(block_seqs,2),size(all_hi_seqs,2));

    while i <= size(block_seqs,2)
        block_seqs_freq(i,:) = hi_freq_band(all_hi_seqs(block_seqs(1,i),:));
        i = i+1;
    end

elseif target_band == 'lo'
    block_seqs_freq = zeros(size(block_seqs,2),size(all_lo_seqs,2));
    while i <= size(block_seqs,2)
        block_seqs_freq(i,:) = lo_freq_band(all_lo_seqs(block_seqs(1,i),:));
        i = i+1;
    end
end

% DISTRACTOR BAND
i=1;
if target_band == 'hi'
    dist_block_seqs_freq = zeros(size(dist_block_seqs,2),size(all_lo_seqs,2));
    while i <= size(dist_block_seqs,2)
        dist_block_seqs_freq(i,:) = lo_freq_band(all_lo_seqs(dist_block_seqs(1,i),:));
        i = i+1;
    end

elseif target_band == 'lo'
    dist_block_seqs_freq = zeros(size(dist_block_seqs,2),size(all_hi_seqs,2));
    while i <= size(dist_block_seqs,2)
        dist_block_seqs_freq(i,:) = hi_freq_band(all_hi_seqs(dist_block_seqs(1,i),:));
        i = i+1;
    end
end


% CREATE AUDIO FILES -----

% First, we'll create a silent gap that is as long as a tone. We'll use
% that to ensure that the target/distractor tones are interleaved.
silent_tone = zeros(1,samples_per_tone);

% Now, we go through and construct the audios for each trial, adding the
% constructed audio into a waveform for the block

% ATTEND BAND
% initialize some variables
miniseq = [];
sasa_block_waveform = [];

i=1; % sequence number within block
while i <= size(block_seqs_freq,1)
    j = 1; % pos within sequence

    while j <= size(block_seqs_freq,2)

        % Choose the right tone, depending on which band is the target

        if target_band == 'hi'
            tone = hi_tones(all_hi_seqs(block_seqs(1,i),j),:);
        elseif target_band == 'lo'
            tone = lo_tones(all_lo_seqs(block_seqs(1,i),j),:);
        end

        switch interleaved
            case 1
                miniseq = [miniseq, tone_gap, tone, tone_gap, silent_tone];
            case 2
                miniseq = [miniseq, tone_gap, silent_tone, tone_gap, tone];
        end

        j = j+1;
    end

    % Put in a silent gap after each sequence
    miniseq = [miniseq, sil_gap];
    % Note how long the trial was
    trial_length = size(miniseq,2)/fs;

    % Save this trial into the target-band waveform for the block
    sasa_block_waveform = [sasa_block_waveform, miniseq];

    % Clear variables
    miniseq = [];

    % And go on to the next trial
    i = i+1;

end

% DISTRACT BAND
% initialize some variables
miniseq = [];
dist_sasa_block_waveform = [];

i=1; % sequence number within block
while i <= size(dist_block_seqs_freq,1)
    j = 1; % pos within sequence

    while j <= size(dist_block_seqs_freq,2)

        % Choose the right tone, depending on which band is the target

        if target_band == 'hi'
            tone = lo_tones(all_lo_seqs(dist_block_seqs(1,i),j),:);
        elseif target_band == 'lo'
            tone = hi_tones(all_hi_seqs(dist_block_seqs(1,i),j),:);
        end

        switch interleaved
            case 1
                miniseq = [miniseq, tone_gap, silent_tone, tone_gap, tone];
            case 2
                miniseq = [miniseq, tone_gap, tone, tone_gap, silent_tone];
        end

        j = j+1;
    end

    % Put in a silent gap
    miniseq = [miniseq, sil_gap];

    % Put this into the distractor-band waveform for this block

    dist_sasa_block_waveform = [dist_sasa_block_waveform, miniseq];

    % Clear variables
    miniseq = [];

    i = i+1;

end

% Now combine the target and distractor-band waveforms
% sasa_attn_dist_block_waveform = sasa_block_waveform + dist_sasa_block_waveform;
% sasa_waveform_length = length(sasa_attn_dist_block_waveform)/fs;

%Creating attn and dist pip locations which will be used for triggers in
%other script
pip_attn_times_diff= diff(sasa_block_waveform~=0);
pip_attn_times_loc = find(pip_attn_times_diff==1); %finding the location of tone pips
% pip_attn_times_loc_everyother = pip_attn_times_loc(:,1:2:end); %indexing every other tone pip 
pip_times_attn = zeros(1,length(pip_attn_times_loc));
pip_times_attn(pip_attn_times_loc) = 1; %setting location of tone pip = 1 

pip_dist_times_diff= diff(dist_sasa_block_waveform~=0);
pip_dist_times_loc = find(pip_dist_times_diff==1);
% pip_dist_times_loc_everyother = pip_dist_times_loc(:,1:2:end);
pip_times_dist = zeros(1,length(pip_dist_times_loc));
pip_times_dist(pip_dist_times_loc) = 1;


% CREATE TIMING FILES -----

switch nrRepTrials
    % If there aren't any repeats, we don't need timing files
    case 0
        %         disp('No repeats, all done');
        %         pip_times = find(diff(tone~=0)==1);

        % But if there are repeats, the earliest a participant can respond is
        % when the final tone of the repeated sequences starts to play
    otherwise

        % Which will depend on whether the target tone is presented before
        % or after the distractor tone
        switch interleaved

            case 1

                for i = 1:length(rep_pos)

                    % start window is from onset of final target tone in
                    % trial
                    rep_pos_times(i,1) = trial_length*(rep_pos(i)-1) + ((seq_length-1)*2)*(size(tone_gap,2) + size(tone,2))/fs;

                    % participants have the length of a trial to respond
                    rep_pos_times(i,2) = rep_pos_times(i,1) + 1.25; %4/22 modified with andrew & chris
                    %goal: get window to be 1.25 sec

                end

            case 2

                for i = 1:length(rep_pos)

                    % start window is from onset of final target tone in
                    % trial
                    rep_pos_times(i,1) = trial_length*(rep_pos(i)-1) + ((seq_length-1)*2+1)*(size(tone_gap,2) + size(tone,2))/fs;

                    % participants have the length of a trial to respond
                    %                     rep_pos_times(i,2) = trial_length*rep_pos(i);
                    rep_pos_times(i,2) = rep_pos_times(i,1) + 1.25;

                end
        end
end

rep_pos_times = rep_pos_times * 1000; % convert timing to ms
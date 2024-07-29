%% Adam Tierney & Victoria Figarola

%Outputs:
%itpc = ITPC value (1x1000; 2Hz is at idx=5)
%phase (1x1000)
function [fftmag,itpc,phase] = ITPC_analysis_final(num_analyze_sweeps,data,nChan,nTrials,fs,subject_id)


% data = ch_hi;
for p = 1:nChan
    index = 1;
    for n = 1:nTrials

        Epoch = data(p,51:end,n);
        [~,b,~] = size(Epoch);
        fftsignalsize = b;

        ramptime = fftsignalsize.*(1/fs);
        ramp = hann(floor(fs*0.04));

        if rem(ramptime, 2)==0   % if even
            % find the mid-point of the ramp.
            z = length(ramp)/2;
            ramp = [ramp(1:z);ones(fftsignalsize-size(ramp,1),1); ramp(z+1:length(ramp))];
        else
            % find the mid-point of the ramp.
            z = floor(length(ramp)/2);
            % ramps up for 7.5 ms  then flat until last 7.5 ms
            ramp = [ramp(1:z);ones(fftsignalsize-size(ramp,1),1);ramp(z+1:length(ramp))];
        end
        clear z
        signal_detr= detrend(Epoch, 'constant');
        signal_ramp =  signal_detr'.*ramp(:,1);
        signal = signal_ramp;
        nowfft = fft(signal,fs*2);
        fftmag = abs(nowfft);

        if index == 1
            all_vectors = zeros(num_analyze_sweeps,length(nowfft));
        end

        all_vectors(index,:) = nowfft./fftmag;
        index = index + 1;
    end

    [one,two] = size(all_vectors);
    final_vector = zeros(1,two);
    for n = 1:one
        final_vector = final_vector + all_vectors(n,:); %summing all trials
    end
    new_vector = final_vector/num_analyze_sweeps;
    % new_vector = new_vector(1:round(length(new_vector)/2)+1);

    itpc(p,:) = abs(new_vector);
    phase(p,:) = angle(new_vector);


end




% figure
% polarplot(test_phase2(5),test_itpc2(5),'r*');
% hold on;
% polarplot(test_phase(5),test_itpc(5),'b*')
%
% phase_both_cond = [test_phase2(5) (test_phase(5)+pi)];
% avg_phase_both_cond = mean(phase_both_cond);
%
% scatter(avg_phase_both_cond,4.78785)
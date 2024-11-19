function [t0,rt] = PresentStimv2(STIM,trig,deviceIndex)

stimchanlist = [1,2,14];
timesPressed = [];
fs = 48000;
StimDur = length(STIM)/fs;

while KbCheck; end % Wait until all keys are released.

% If the user has pressed a key, then display its code number and name.
KbQueueStart(deviceIndex); % starts delivering keypresses to the queue
t0 = GetSecs;
TimeElapsed = 0;
pg = playrec('play',[STIM,trig],stimchanlist); % start auditory playback

while TimeElapsed < StimDur
    %     timePressed = KbQueueWait(deviceIndex,0);
    [pressed,timePressed] = KbQueueCheck(deviceIndex);
    if pressed
        timesPressed= [timesPressed,timePressed];
    end
    %     playrec('block',pg);
    
    % obtains keypress data since queuestart or queuewait
    %     checkBit = 0;
    TimeElapsed = GetSecs-t0;
end

playrec('block',pg);
rt = timesPressed(timesPressed~=0)-t0;
% rt = timesPressed-t0;

end
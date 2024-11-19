                                                                                                                                                                     clear all; close all; clc;
commandwindow;
tic;

%% Initialize sound card
% Step One: Connect to and properly initialize sound card
% pr = genpath('C:\Program Files\PlayrecForMatlab-playrec-7c15bf6');
% addpath(pr);
addpath(genpath('..\Frog-quency'));
fprintf('Initializing connection to sound card...\n')

%% Initialize sound card (RME)
% Step One: Connect to and properly initialize RME sound card
fprintf('Initializing connection to sound card...\n')
Devices=playrec('getDevices');
if isempty(Devices)
    error(sprintf('There are no devices available using the selected host APIs.\nPlease make sure the RME is powered on!')); %#ok<SPERR>
else
    i=1;
    while ~strcmp(Devices(i).name,'ASIO MADIface USB') && i <= length(Devices)
        i=i+1;
    end
end
stimchanList=[1,2];
% fs = Devices(i).defaultSampleRate;
fs = 48000;
playDev = Devices(i).deviceID;
playrec('init',fs,playDev,-1,14,-1);
fprintf('Success! Connected to %s.\n', Devices(i).name);

todayStr = datestr(now,'yyyymmdd');
inputInfo=inputdlg({'subject ID: '});
sID = inputInfo{1};

streamflag = 1; %Define whether the deviant tones will be time-displaced or have different frequencies
while streamflag == 1
    stream_type = input('Please enter which stream the subject will be attending to first (b for bull or f for frog):', 's');
    switch stream_type
        case {'b','B','bull','Bull','BULL'}
            stream_type = 1;
            streamflag = 0;
            conds = ["Lo","Hi";"left","right"];

%             HighorLow = ['Hi','Hi','Lo','Lo'];
%             Side = ['right','left','right','left'];
        case {'f','F','frog','Frog','FROG'}
            stream_type = 2;
            streamflag = 0;
            conds = ["Hi","Lo";"right","left"];

%             HighorLow = ['Hi','Lo','Hi','Lo'];
%             Side = ['right','right','left','left'];
        otherwise
            fprintf(2, 'Unrecognized answer! Try again!');
    end
end

AVsave = 1;

while AVsave == 1
    prompt = {'IS ACTIVIEW SAVING?!'};
    dlgtitle = 'Input';
    dims = [1 50];
    definput = {'y/n'};
    answer = inputdlg(prompt,dlgtitle,dims,definput);
    Actiview = answer{1};
    
    switch Actiview
        case {'Y', 'y', 'yes', 'Yes', 'YES'}
            AVsave = 0;
        otherwise
            fprintf(2, 'Unrecognized answer! Try again!');
    end
end


%% Define experimental parameters
blocknum = 16; % Number of blocks
num_cond = 2;
% blocknum = 1; % Number of blocks
% num_cond = 1;
Subj_resp = cell(num_cond,blocknum); % Initialize the subject response array
Responses = cell(num_cond,blocknum);
Hit_Boxes = cell(num_cond,blocknum);
t0 = zeros(num_cond,blocknum); % Initialize the start times
elapsedtime = [];

audiofiles = struct2cell(dir('WavFiles\'))';
% audiofiles = struct2cell(dir('WavFiles_testing\'))';
audionames = audiofiles(3:end,1);

for condition_counter = 1:num_cond
    audio_freq = audionames(contains(audionames,conds(1,condition_counter)));
    audio_side = audionames(contains(audionames,conds(2,condition_counter)));
    desired_audio(:,condition_counter) = intersect(audio_side,audio_freq);
end

%% Initialize Psychtoolbox
PsychDefaultSetup(2);
Screen('Preference','DefaultFontSize',50);
Screen('Preference','SyncTestSettings', 0.002); % Rather than skipping the sync test, we increase the tolerance to ±2 ms.
Screen('Preference','VisualDebugLevel',1);
Screen('Preference', 'TextAntiAliasing', 2);
Screen('Preference', 'TextRenderer', 1);
sn = max(Screen('Screens'));

% % Define black, white and grey
black = BlackIndex(sn);
white = WhiteIndex(sn);
grey = [0.6 0.6 0.6];
HideCursor;

[window, winRect] = PsychImaging('OpenWindow', sn, black);
Screen('BlendFunction', window, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
[width,height] = Screen('WindowSize',window);
[scrX,scrY] = RectCenter(winRect);
Screen('TextFont', window, 'Helvetica');

% Set up keyboard
deviceIndex = GetKeyboardIndices(1);
deviceIndex = [];
KbName('UnifyKeyNames');

keysOfInterest=zeros(1,256);
keysOfInterest(KbName('space'))=1;
keysOfInterest(KbName('ESCAPE'))=1;

KbQueueCreate(deviceIndex,keysOfInterest);

curText = ['<color=ffffff>In this task, you will be listening to streams of'...
    ' six-tone patterns. You will either be attending to a high frequency band '...
    ' or a low frequency band, and will have to determine whether they contain '...
    ' the same pair of 3-tones than the previous pattern. The band you will be'...
    ' attending to will occur in isolation at the beginning of the block, and will'...
    ' determine which ear you will attending that band in. There are 16 blocks'...
    ' for attending to the low band, and 16 blocks for attending to the high band.'...
    ' Press [SPACE] when you hear a repeated pattern.'...
    ' \n\nWave your hand if you need the instructions clarified.'...
    ' Press <b>[Space]<b> to continue.'];

DrawFormattedText2(curText,'win',window,'sx',140,'sy','center','xalign','left','yalign','center','wrapat',75);
Screen('Flip',window);
KbWait([],2); %Wait for a keystroke

triggers_attn = cell(num_cond,blocknum);
triggers_dist = cell(num_cond,blocknum);
alltriggeraudio = cell(num_cond,blocknum);

%% Change screens and present stimuli
for ncon = 1:num_cond
    for block = 1: blocknum
        % Open an on screen window and color it grey
        % Draw text in the upper portion of the screen with the default font in red
        % On-screen instructions
        Screen('TextSize',window, 60);
        curText = strcat(sprintf('Block # %d ',block));
        DrawFormattedText(window,curText,'center',scrY*0.6,0.6*ones(1,3));
        Screen('TextSize',window,45);
        DrawFormattedText(window,'Press [SPACE] when you hear a repeated pattern','center',scrY*0.76,0.6*ones(1,3));
        Screen('Flip',window);
        
        %% Generate stimulus audio
        fprintf('Generating stimulus audio...');
        
         audiofilename = char(desired_audio(block,ncon));
        [stimaudio,fs] = audioread(['WavFiles/' audiofilename]);
%         [stimaudio,fs] = audioread(['WavFiles_testing/' audiofilename]);
        trigfilename_attn = strcat(audiofilename(1:14),num2str(block),'-attn_trig.wav');
        trigfilename_dist = strcat(audiofilename(1:14),num2str(block),'-dist_trig.wav');
%         [trigaudio_attn,fs] = audioread(['WavFiles_testing/' trigfilename_attn]);
%         [trigaudio_dist,fs] = audioread(['WavFiles_testing/' trigfilename_dist]);
        [trigaudio_attn,fs] = audioread(['WavFiles/' trigfilename_attn]);
        [trigaudio_dist,fs] = audioread(['WavFiles/' trigfilename_dist]);

        trigaudio_dist = trigaudio_dist*100;
        trigaudio_attn = trigaudio_attn*100;

        % make the first trigger of the block show the block number and all the
        % others after that just a 1 (for each trial)
        % filter with ones to give the triggers sufficient duration
        side_attn = isempty(strfind(audiofilename,'left'));
        side_dist = ~side_attn;
        
        %Adds in the trigger values for 5,7,13,15
        trigaudio_attn(trigaudio_attn~=0) = trigaudio_attn(trigaudio_attn~=0)+4*side_attn;
        trigaudio_dist(trigaudio_dist~=0) = trigaudio_dist(trigaudio_dist~=0)+4*side_dist;

        trigaudio_attn = round(trigaudio_attn);
        trigaudio_dist = round(trigaudio_dist);

        triggers_attn{ncon,block} = trigaudio_attn;
        triggers_dist{ncon,block} = trigaudio_dist;
               
        trigaudio_attn(trigaudio_attn~=0) = trignum2scalar(trigaudio_attn(find(trigaudio_attn ~= 0,1))); %
        trigaudio_dist(trigaudio_dist~=0) = trignum2scalar(trigaudio_dist(find(trigaudio_dist ~= 0,1))); %
        trigaudio_attn = [trigaudio_attn;zeros(length(stimaudio)-length(trigaudio_attn),1)];
        trigaudio_dist = [trigaudio_dist;zeros(length(stimaudio)-length(trigaudio_dist),1)];

        trigaudio = trigaudio_attn + trigaudio_dist;

        trigaudio = filter(ones(480, 1), 1, trigaudio); % Recording at 10kHz, as long as were recording >2kHz, we're good
        
%         repeattimes = table2cell(readtable(strcat('attend',conds(1,ncon),'-block_repeat_times.txt')));
%         repeattimes = strrep(repeattimes,char(91),'');
        
        fprintf('DONE\n');
        
        %% Experimental trial mock up
        
        startTime = GetSecs;
        
%         [T0,PressT] = PresentStim(stimaudio,trigaudio,deviceIndex);
        [T0,PressT] = PresentStimv2(stimaudio,trigaudio,deviceIndex);
        elapsedtime(block) = GetSecs; % Duration of each trial
        t0(ncon,block) = T0; % Time of the beginning of each trial
        Subj_resp{ncon,block} = PressT; % Time of the click
        
        alltriggeraudio{ncon,block} = trigaudio;
        
        endTime = GetSecs;
        KbQueueStart(deviceIndex);
        runTime = endTime-startTime;
        
        %     BaseTime = t0+t(Dev_locs{diff,block}(end),dev_pos); % Adds the time for the deviant to be played
        
        %     Hit_Boxes{block}(:,1) = t(Dev_locs{block}); % Defines the beginning of the hitbox
        %     Hit_Boxes{block}(:,2) = t(Dev_locs{block})+response_window; % Defines the end of the hitbox
        
        Responses{ncon,block} = Responses{ncon,block}-t0(ncon,block);

        if ncon == num_cond && block == blocknum
            curText = ['<color=ffffff>Experiment complete!'];
            DrawFormattedText2(curText,'win',window,'sx',140,'sy','center','xalign','left','yalign','center','wrapat',75);
            Screen('Flip',window);
            KbWait([],2); %Wait for a keystroke
            KbQueueWait(deviceIndex,1); %Returns once no keys are down
        elseif block == blocknum
            curText = ['<color=ffffff>Condition complete. Feel free to take a break. \nPress <b>[SPACE]<b> to continue'];
            DrawFormattedText2(curText,'win',window,'sx',140,'sy','center','xalign','left','yalign','center','wrapat',75);
            Screen('Flip',window);
            KbWait([],2); %Wait for a keystroke
            KbQueueWait(deviceIndex,1); %Returns once no keys are down
        else
            curText = ['<color=ffffff>Block complete. Feel free to take a break. \nPress <b>[SPACE]<b> to continue'];
            DrawFormattedText2(curText,'win',window,'sx',140,'sy','center','xalign','left','yalign','center','wrapat',75);
            Screen('Flip',window);
            KbWait([],2); %Wait for a keystroke
            KbQueueWait(deviceIndex,1); %Returns once no keys are down
        end

        save(strcat(sID,'_EEG_behav'))
% toc;
%     toc;
    end
end

%%

Screen('FillRect',window,black);
Screen('Flip',window);

KbQueueRelease(deviceIndex);
WaitSecs(2);
ShowCursor;
ListenChar(0);
sca;

%% Save

runtime=toc;
save(strcat(sID,'_EEG_behav'))
%         save(strcat(sID,'_condition_and_block_order_',todayStr))

fprintf('Total time elapsed %d min %d sec.\n',floor(runtime/60),floor(mod(runtime,60)));

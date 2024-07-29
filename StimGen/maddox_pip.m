function mpip = maddox_pip(f_pip,fs, amp)
% f_pip = 3500;  % tone frequency
% fs=48000;
% amp = 1;
n_cycles = 5; % this determines length of pip along w freq, obviously
sign_flip = false; % no need to invert here

len_pip = round(n_cycles / f_pip * fs);
t = (0:len_pip - 1).' / fs;

mpip = amp* blackman(len_pip) .* cos(2 * pi * f_pip * t);
mpip = amp .* cos(2 * pi * f_pip * t);

if sign_flip
    mpip = mpip * -1;
end
% 
% figure()
% plot(t,mpip,'b','LineWidth',3);
% set(gca,'FontSize',12,'FontWeight','bold')
% axis off
% 


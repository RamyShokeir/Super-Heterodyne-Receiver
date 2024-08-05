clear all;
% Load audio signals
[Signal1, Fs1] = audioread('Short_BBCArabic2.wav');
[Signal2, Fs2] = audioread('Short_FM9090.wav');

% Combine stereo channels
Signal1 = sum(Signal1, 2);
Signal2 = sum(Signal2, 2);

% Upsample signals if necessary
maxSamples = 1e5;
if length(Signal1) > maxSamples
    targetFs = 10 * Fs1;  % Set your desired new sampling frequency
    Signal1 = interp(Signal1, 10);
    Fs1 = targetFs;
end

if length(Signal2) > maxSamples
    targetFs = 10 * Fs2;  % Set your desired new sampling frequency
    Signal2 = interp(Signal2, 10);
    Fs2 = targetFs;
end

% Interpolate signals
%interpFactor = 10;
%Signal1 = interp(Signal1, interpFactor);
%Signal2 = interp(Signal2, interpFactor);

% Equalize lengths of both signals
maxLength = max(length(Signal1), length(Signal2));
Signal1 = padarray(Signal1, maxLength - length(Signal1), 0, 'post');
Signal2 = padarray(Signal2, maxLength - length(Signal2), 0, 'post');

L1=length(Signal1);
L2=length(Signal2);
Signal1_truncated = Signal1(1:L1);
Signal2_truncated = Signal2(1:L2);
%Equalize legnths of both signals
maxLength = max(length(Signal1), length(Signal2));
Signal1 = padarray(Signal1, maxLength - length(Signal1) ,0,'post');
Signal2 = padarray(Signal2, maxLength - length(Signal2) ,0,'post');

carrier_frequency1 = 100e3; % Carrier frequency (100 kHz)
carrier_frequency2 = 155e3;
% Modulate the signal onto the carrier
Fs_c = 800e3;
t = (0:1/Fs_c:(length(Signal1)-1)/Fs_c).';
%modulated_signal1 = Signal1 .* Carrier_signal;
carrierSignal1 = cos(2 * pi * carrier_frequency1 * t);
carrierSignal2 = cos(2 * pi * carrier_frequency2 * t);
Nc1 = length(carrierSignal1); % Length of the signal
frequencies_c1 = Fs_c*(-Nc1/2:Nc1/2-1)/Nc1; % Frequency axis
fft_result_c1 = fftshift(fft(carrierSignal1));
magnitude_spectrum_c1 = abs(fft_result_c1);

Nc2 = length(carrierSignal2); % Length of the signal
frequencies_c2 = Fs_c*(-Nc2/2:Nc2/2-1)/Nc2; % Frequency axis
fft_result_c2 = fftshift(fft(carrierSignal2));
magnitude_spectrum_c2 = abs(fft_result_c2);

modulatedSignal1 = Signal1 .* carrierSignal1;
modulatedSignal2 = Signal2 .* carrierSignal2;

N_mod1 = length(modulatedSignal1); % Length of the signal
frequencies_mod1 = Fs_c*(-N_mod1/2:N_mod1/2-1)/N_mod1; % Frequency axis
fft_result_mod1 = fftshift(fft(modulatedSignal1));
magnitude_spectrum_mod1 = abs(fft_result_mod1);

N_mod2 = length(modulatedSignal2); % Length of the signal
frequencies_mod2 = Fs_c*(-N_mod2/2:N_mod2/2-1)/N_mod2; % Frequency axis
fft_result_mod2 = fftshift(fft(modulatedSignal2));
magnitude_spectrum_mod2 = abs(fft_result_mod2);

%FDM Signal Construction
fdmSignal = modulatedSignal1 + modulatedSignal2;
N_fdm = length(fdmSignal); % Length of the signal
frequencies_fdm = Fs_c*(-N_fdm/2:N_fdm/2-1)/N_fdm; % Frequency axis
fft_result_fdm = fftshift(fft(fdmSignal));
magnitude_spectrum_fdm = abs(fft_result_fdm);

subplot(7,1,3);
plot(frequencies_c1,magnitude_spectrum_c1);
title('Carrier Signal1');
xlabel('Time (seconds)');
ylabel('Amplitude');

subplot(7,1,4);
plot(frequencies_c2,magnitude_spectrum_c2);
title('Carrier Signal2');
xlabel('Time (seconds)');
ylabel('Amplitude');

subplot(7,1,5);

plot(frequencies_mod1, magnitude_spectrum_mod1);
title('Magnitude Spectrum Modulated1');
xlabel('Frequency (Hz)');
ylabel('Magnitude');

subplot(7,1,6);

plot(frequencies_mod2, magnitude_spectrum_mod2);
title('Magnitude Spectrum Modulated2');
xlabel('Frequency (Hz)');
ylabel('Magnitude');

subplot(7,1,7);
figure;
plot(frequencies_fdm, magnitude_spectrum_fdm);
title('Magnitude Spectrum FDM Signal');
xlabel('Frequency (Hz)');
ylabel('Magnitude');

% Step 7: RF Stage (BPF) // Signal1
centerFreqRF1 = carrier_frequency1;
bandwidthRF1 = 20e3; % Adjust as needed
bpfRF1 = design(fdesign.bandpass('N,Fc1,Fc2', 100, centerFreqRF1 - bandwidthRF1/2, centerFreqRF1 + bandwidthRF1/2, Fs_c));

rfOutput1 = filter(bpfRF1, fdmSignal);

%Signal2
%centerFreqRF2 = carrier_frequency2;
%bandwidthRF2 = 18e3; % Adjust as needed
%bpfRF2 = design(fdesign.bandpass('N,Fc1,Fc2', 100, centerFreqRF2 - bandwidthRF2/2, centerFreqRF2 + bandwidthRF2/2, Fs_c));

%rfOutput2 = filter(bpfRF2, fdmSignal);

Nrf2 = length(rfOutput1); % Length of the signal
frequencies_rf2 = Fs_c*(-Nrf2/2:Nrf2/2-1)/Nrf2; % Frequency axis
fft_result_rf2 = fftshift(fft(rfOutput1));
magnitude_spectrum_rf2 = abs(fft_result_rf2);

figure;
plot(frequencies_rf2, magnitude_spectrum_rf2);
title('Magnitude Spectrum RF2');
xlabel('Frequency (Hz)');
ylabel('Magnitude');


% Step 8: Oscillator and Mixer // Signal1
IFFreq1 = 27.5e3;
tIFF1= ((1:length(fdmSignal)) / Fs_c).';
oscillatorSignal1 = cos(2 * pi * (carrier_frequency1 + IFFreq1) *tIFF1);
mixedSignal1 = fdmSignal .* oscillatorSignal1;

%Signal2
%IFFreq2 = 27.5e3;
%tIFF2= ((1:length(rfOutput2)) /Fs_c ).';
%oscillatorSignal2 = cos(2 * pi * (carrier_frequency2 + IFFreq2) *tIFF2);
%mixedSignal2 = rfOutput2 .* oscillatorSignal2;

Nmix2 = length(mixedSignal1); % Length of the signal
frequencies_mix2 = Fs_c*(-Nmix2/2:Nmix2/2-1)/Nmix2; % Frequency axis
fft_result_mix2 = fftshift(fft(mixedSignal1));
magnitude_spectrum_mix2 = abs(fft_result_mix2);
figure;
plot(frequencies_mix2, magnitude_spectrum_mix2);
title('Magnitude Spectrum MIX2');
xlabel('Frequency (Hz)');
ylabel('Magnitude');

%IF Stage BPF
centerFreqIF = IFFreq1;
bandwidthRF2 = 18e3; % Adjust as needed
bpfIF = design(fdesign.bandpass('N,Fc1,Fc2', 100, centerFreqIF - bandwidthRF2/2, centerFreqIF + bandwidthRF2/2, Fs_c));
IfOutput = filter(bpfIF, mixedSignal1);
NIf = length(IfOutput); % Length of the signal
frequencies_If = Fs_c*(-NIf/2:NIf/2-1)/NIf; % Frequency axis
fft_result_If = fftshift(fft(IfOutput));
magnitude_spectrum_If = abs(fft_result_If);

figure;
plot(frequencies_If, magnitude_spectrum_If);
title('Magnitude Spectrum IF');
xlabel('Frequency (Hz)');
ylabel('Magnitude');

%Final Stage Signal to BaseBand
% Step 10: Baseband Detection
tBB2=(0:1/Fs_c:(length(IfOutput)-1)/Fs_c).';
basebandCarrier2 = cos(2 * pi * IFFreq1 * tBB2);
basebandOutput2 = IfOutput .* basebandCarrier2;

lpf = design(fdesign.lowpass('N,Fc', 100, 22.5e3, Fs_c));
finalOutput = filter(lpf, basebandOutput2);
finalOutput = finalOutput / max(abs(finalOutput));
%finalOutput=resample(finalOutput,44100,441000);
Nfinal = length(finalOutput); % Length of the signal
frequencies_finalOutput = Fs_c*(-Nfinal/2:Nfinal/2-1)/Nfinal; % Frequency axis
fft_result_final = fftshift(fft(finalOutput));
magnitude_spectrum_final = abs(fft_result_final);
plot(frequencies_finalOutput, magnitude_spectrum_final);
title('Magnitude Spectrum FinalOutput');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
sound(finalOutput,44100);
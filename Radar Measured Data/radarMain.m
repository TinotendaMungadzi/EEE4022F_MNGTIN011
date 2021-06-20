
%% Processing recordings

clear all;
close all;
                

wavFile_CW_All = {'Audi_A1_Driving_Away_30KPH.wav'; 
                  'Audi_A1_Driving_Away_45KPH.wav';
                  'Audi_A1_Driving_Towards_15KPH_No_Slowing.wav'};
              
RecordingNo2Process = 1;              

wavFile = wavFile_CW_All{RecordingNo2Process};
% cantenna_dop_v3_yunus(wavFile_CW);
 
[rows columns] = size(wavFile)
% ---- setup constants and parameters ----
c = 299e6; % (m/s) speed of light
cpi = 0.067; % (s) coherent processing interval 
fc = 2590e6; % (Hz) Center frequency (connect VCO Vtune to +5)
maxSpeed = 50; % (m/s) maximum speed to display
overlapFactor = 3; % (unitless) number of overlapped pulse windows (1 for no overlap)
ovsDop = 5; % (unitless) oversample factor for Doppler axis
% ----- end constants and parameters -----

% use a default filename if none is given
if ~exist('wavFile','var')
    wavFile = 'radar_test2.wav';
end

% read the raw wave data
fprintf('Loading WAV file...\n');
[Y,Fs] = audioread(wavFile,'native');

% derived parameters
N = round(cpi*Fs) % # of samples per pulse

% the input appears to be inverted
x = -Y(:,2); % Received signal at baseband
clear Y;

% grab an integer number of overlapped frames
M = floor(numel(x) / N * overlapFactor) - (overlapFactor) + 1;

% compute axes parameters for the plot
% Note: the Doppler data is oversampled by ovsDop
delta_f = (0:ovsDop*N/2-1).' / (ovsDop*N) * Fs; % Doppler freq. (Hz)
lambda = c / fc; % wavelength (m)
speed = delta_f * lambda / 2; % Doppler -> speed
time = (1:M) * cpi / overlapFactor; % collection time (sec)

% limit the speed axis to a reasonable range
speed = speed(speed <= maxSpeed);
nSpeed = numel(speed);

% compute a Doppler window
dopWin = hann_window(N);

% compute the Doppler vs. time plot
fprintf('Processing...\n');
dti = zeros(nSpeed, M);
DataAfterPowerLawDetector = zeros(nSpeed, M); 

for mIdx = 1:M
    xInds = (1:N).' + (mIdx-1)*floor(N/overlapFactor); % compute indices
    tmp = double(x(xInds)) .* dopWin; % apply Doppler window
    subplot(2,1,1)
    plot(double(x(xInds)))
    subplot(2,1,2)
    plot(tmp)
    tmp = tmp - mean(tmp); % remove DC component if it exists
    tmp = fft(tmp, ovsDop*N); % compute oversampled Doppler spectrum
    dti(:,mIdx) = 20*log10( abs(tmp(1:nSpeed)) ); % grab result in dB
    DataAfterPowerLawDetector(:,mIdx) =  abs(tmp(1:nSpeed)).^2 ; % result in linear for processing
end
clear x;
dti = dti.';

% make Doppler vs. time plot
figure;
imagesc(time, speed, dti');
%imagesc(speed,time,dti);
colormap(jet(256));
caxis(max(dti(:)) + [-60 0]); % show 60 dB dynamic range
colorbar;
ylabel('Speed (m/sec)');
xlabel('Time (sec)');
axis xy;

%% Detection
% CA-CFAR DETECTION CODE

% Parameters
PFA = 10^-3; %Prob of false alarm
% RefWindow = 12;                                       % Dr Abdul Gaffar
% N = RefWindow;                                        % Dr Abdul Gaffar 
guardLength = 2;
SignalLength = length(DataAfterPowerLawDetector);
Window_Size = 20;                                       
N = Window_Size*2;                                      % Dr Abdul Gaffar
%N = 2*(Window_Size/2+guardLength/2)+1;
                  

[rows columns] = size(DataAfterPowerLawDetector)

SFAlpha = N*(PFA^(-1/N)-1);

% z
%DataAfterPowerLawDetector = abs(y_complex).^2 % z

gca = [];
detectionPosX = [];
detectionPosY = [];
NumberOfDetections = 0;

for a = 1:columns
    colData = DataAfterPowerLawDetector(1:rows,a);
    G = [];
    
    for i = Window_Size+(guardLength/2)+1:(rows)-(Window_Size+(guardLength/2))  
       
        CUT_Power = colData(i);
        FLag = colData(i-Window_Size-(guardLength/2):i-(guardLength/2)-1);     % Lagging Window
        FLead = colData(i+1+(guardLength/2):i+(guardLength/2)+Window_Size);    % Leading Window
        
        AvgRefCells = (mean(FLag) + mean(FLead))/2;
        G = [G; AvgRefCells];

   %     T = SFAlpha.*gca;                          % Dr Abdul Gaffar
         T = SFAlpha.*AvgRefCells;                  % Dr Abdul Gaffar 
                
        if (T<CUT_Power)
            detectionPosX = [detectionPosX; i];              % Add row position of detection 
            detectionPosY = [detectionPosY; a];              % Add column position of detection 
            NumberOfDetections = NumberOfDetections + 1;     % Count number of targets 
       end

                
    end
    gca = [gca, G];
end



% make Doppler vs. time plot

detectionToPlot_PosY = detectionPosX;              % Dr Abdul Gaffar
detectionToPlot_PosX = detectionPosY;              % Dr Abdul Gaffar 

detection_speed =  speed(detectionToPlot_PosY);     % Dr Abdul Gaffar 
detection_time  =  time(detectionToPlot_PosX);      % Dr Abdul Gaffar

% Plot detections over the data
figure;
% imagesc(dti.');
imagesc(time, speed,dti.');                            % Dr Abdul Gaffar: plot with speed and time as axes
colormap(jet(256));
caxis(max(dti(:)) + [-60 0]); % show 60 dB dynamic range
colorbar;
ylabel('Speed (m/sec)');
xlabel('Time (sec)');
axis xy;
hold on;
plot(detection_time,detection_speed,'kx', 'MarkerSize',8, 'LineWidth',2);
grid on;

% plot speed estimates
figure;
plot(detection_time, detection_speed, 'x'); 
xlabel('Time (sec)');
ylabel('Speed (m/sec)');
grid on;

% PFA error
NumberOfDetections
SignalLength
PFA_expected = PFA;
PFA_Simulation = NumberOfDetections/(rows*columns)
PFA_new = (1+SFAlpha/N)^(-N)
PFA_error_percentage = (PFA_expected - PFA_Simulation)/PFA_expected * 100

%plot(hann_window(N))
% 
% T = SFAlpha.*gca
% Threshold = pow2db(T(N/2,:))
% SignalPowerdB = pow2db(DataAfterPowerLawDetector(N/2,:));
% SNR = snr(DataAfterPowerLawDetector,Fs)
%   
% size = length(Threshold);
%     
% detection = [];

% plot
% figure
% xaxis = 0:size-1;
% SignalPowerdB(size+1:SignalLength) = 0;
% plot(xaxis, SignalPowerdB(1:size), xaxis, Threshold);
% % ylim([-40 40]);
% xlabel('Samples');
% ylabel('Power (dB)');
% title('CA-CFAR');
% legend('Signal','Threshold')

%% ---- standard DSP helper functions below ----


function [w] = hann_window(N)
% create a hann (cosine squared) window
w = .5 + .5*cos(2*pi*((0:N-1).'/(N-1) - .5));

end


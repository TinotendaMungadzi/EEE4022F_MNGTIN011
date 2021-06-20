% CA-CFAR DETECTION CODE
close all;
clear;

% simulated data
x = rand(1,1000);

% Parameters
PFA = 10^-3; 
RefWindow = 12;
N = 2*RefWindow; 
guardLength = 2;
SignalLength = length(x);
Window_Size = RefWindow;

%SNR_dB = 20;
%SNR_Linear = 10^(SNR_dB/10);

SFAlpha = N*(PFA^(-1/N)-1);

% Noise
yReal =  normrnd(1,10, [1,1000]);
yImag = 1i*normrnd(1,10, [1,1000]);
y_complex = yReal + yImag %y

% z
DataAfterPowerLawDetector = abs(y_complex).^2 % z

gca = [];
test = Window_Size+(guardLength/2)+1
test2 = (SignalLength)-(Window_Size+(guardLength/2))

for i = Window_Size+(guardLength/2)+1:(SignalLength)-(Window_Size+(guardLength/2))  
    tempG = [];
    B = DataAfterPowerLawDetector;
    CUT_Power = B(i);
    
    % f = i-Window_Size-(guardLength/2)
    % l = i-(guardLength/2)-1
    
    FLag = B(i-Window_Size-(guardLength/2):i-(guardLength/2)-1);     % Lagging Window
    FLead = B(i+1+(guardLength/2):i+(guardLength/2)+Window_Size);    % Leading Window
    
    AvgRefCells = (mean(FLag) + mean(FLead))/2;
    tempG = [tempG; AvgRefCells];
    
    gca = [gca, tempG];
end

T = SFAlpha.*gca
Threshold = 10*log10(T)
SignalPowerdB = 10*log10(x + DataAfterPowerLawDetector)

size = length(Threshold)

detection = [];
NumberOfDetections = 0;
for j = 1:size
    if T(j)<SignalPowerdB(j)
        detection(j)=1;
        NumberOfDetections = NumberOfDetections+1;
    else
        detection(j)=0;
    end
end

NumberOfDetections

% PFA error
PFA_expected = PFA;
PFA_Simulation = (1+SFAlpha);
PFA_error = (PFA_expected - PFA_Simulation)/PFA_expected

figure(1)
n = 0:size-1;
SignalPowerdB(size+1:SignalLength) = 0;
plot(n, SignalPowerdB(1:size), n, Threshold);
ylim([-40 40]);
xlabel('Resolution Bin Index');
ylabel('Power (dB)');
title('CA-CFAR');
legend('Signal','Threshold')

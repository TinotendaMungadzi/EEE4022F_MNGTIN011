% CA-CFAR DETECTION CODE
close all;
clear;

% simulated data
columns = 794;
rows = 208;

x = zeros(rows,columns);
x([40, 50, 70, 200, 500, 900])=[3000 4000 5000 6000 7000 80000]
figure,plot(x);

% Parameters
PFA = 10^-3; %Prob of false alarm
RefWindow = 32;
N = 2*RefWindow; 
guardLength = 4;
SignalLength = length(x);
Window_Size = 32;

SFAlpha = N*(PFA^(-1/N)-1);

yReal =  normrnd(0,10, [rows,columns]);
yImag = 1i*normrnd(0,10, [rows,columns]);
y_complex = yReal + yImag %y

% z
DataAfterPowerLawDetector = abs(y_complex).^2 % z

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
        
        T = SFAlpha.*gca;
        
        if (T<CUT_Power)
            detectionPosX = [detectionPosX; i];      % Add row position of detection 
            detectionPosY = [detectionPosY; a];      % Add column position of detection 
            NumberOfDetections = NumberOfDetections + 1;                   % Count number of targets 
        end
                
    end
    gca = [gca, G];
end


% PFA error

PFA_expected = PFA;
PFA_Simulation = NumberOfDetections/SignalLength;
PFA_error = ((PFA_expected - PFA_Simulation)/PFA_expected)*100

Threshold = pow2db(T(N/2,:));
SignalPowerdB = pow2db(DataAfterPowerLawDetector(N/2,:));
  
size = length(Threshold);
    
% plot
figure
xaxis = 0:size-1;
SignalPowerdB(size+1:SignalLength) = 0;
plot(xaxis, SignalPowerdB(1:size), xaxis, Threshold);
ylim([-40 40]);
xlabel('Samples');
ylabel('Power (dB)');
title('CA-CFAR');
legend('Signal','Threshold')

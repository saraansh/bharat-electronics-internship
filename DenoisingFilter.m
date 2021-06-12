%% Setting up the variables before jumping into processing.
WinLenSec = 0.0050; % Window length in seconds.
HopPercent = 1; % percentage of hopping.
AROrder = 20; % Auto regressive filter order.
NumIter = 7;

%% Reading Input signal and creating noisy input as well.
[Input, Fs] = audioread('input1.wav');
Input = Input(:,1);
Noise = normrnd(0,sqrt(0.01),size(Input));
audiowrite('noise.wav', Normalize(Noise), Fs); 
NoisyInput = Input + Noise;
audiowrite('noisy_input.wav', Normalize(NoisyInput), Fs);
Time = (0:1/Fs:(length(Input)-1)/Fs)';

% Chopping session.
WinLenSamples = fix(WinLenSec * Fs);
Window = ones(WinLenSamples,1);
[ChoppedSignal, NumSegments] = Chopper(NoisyInput, WinLenSamples, Window, HopPercent);

% Matrix initializations.
H = [zeros(1,AROrder-1),1];   % Measurement matrix.
R = var(Noise);     % Variance of noise.

[FiltCoeff, Q] = lpc(ChoppedSignal, AROrder);   % Finding filter coefficients.

P = R * eye(AROrder,AROrder);   % Error covariance matrix.
Output = zeros(1,size(NoisyInput,1));   % Allocating memory for output signal.
Output(1:AROrder) = NoisyInput(1:AROrder,1)';   % Initializing output signal according to equation (13)
OutputP = NoisyInput(1:AROrder,1);

% Iterators.
i = AROrder+1;
j = AROrder+1;

% Processing.
for k = 1:NumSegments   % For every segment of chopped signal...
    jStart = j;     % Keeping track of AROrder+1 value for every iteration.
    OutputOld = OutputP;    % Keeping the first AROrder amount of samples for every iteration.
    
    for l = 1:NumIter
        A = [zeros(AROrder-1,1) eye(AROrder-1); fliplr(-FiltCoeff(k,2:end))];
        
        for ii = i:WinLenSamples
            OutputC = A * OutputP;
            Pc = (A * P * A') + (H' * Q(k) * H);
            K = (Pc * H')/((H * Pc * H') + R);
            OutputP = OutputC + (K * (ChoppedSignal(ii,k) - (H*OutputC)));
            Output(j-AROrder+1:j) = OutputP';
            P = (eye(AROrder) - K * H) * Pc;
            j = j+1;
        end
        
        i = 1;
        if l < NumIter
            j = jStart;
            OutputP = OutputOld;
        end
        
        % update lpc on filtered signal
        [FiltCoeff(k,:), Q(k)] = lpc(Output((k-1)*WinLenSamples+1:k*WinLenSamples),AROrder);
    end
end
Output = Output';
audiowrite('output.wav', Normalize(Output), Fs);

% Sound-to-Noise-Ratio
snr(NoisyInput, NoisyInput - Input)
snr(Output, Output - Input)

%% Plotting the results
figure
ax(1) = subplot(2,2,1)
plot(Time, Input)
xlabel('Time in seconds')
ylabel('Amplitude')
title('Clean speech signal')
ax(2) = subplot(2,2,2)
plot(Time, Noise)
xlabel('Time in seconds')
ylabel('Amplitude')
title('Generated noise')
ax(3) = subplot(2,2,3)
plot(Time, NoisyInput)
xlabel('Time in seconds')
ylabel('Amplitude')
title('Noisy input signal')
ax(4) = subplot(2,2,4)
plot(Time, Output)
xlabel('Time in seconds')
ylabel('Amplitude')
title('Estimated clean output signal')
axis(ax(1:4), [0 2 -1.25 1.25]);

%% Correlation Plots
[C1,lag1] = xcorr(Input, Output);
[C2,lag2] = xcorr(Noise, Output);
figure
ax(1) = subplot(2, 1, 1);
plot(lag1/Fs, C1, 'r');
ylabel('Amplitude')
title('Cross-correlation b/w input and output signal')
ax(2) = subplot(2, 1, 2);
plot(lag2/Fs, C2, 'r');
ylabel('Amplitude')
title('Cross-correlation b/w noise and output signal')
xlabel('Time(secs)');
axis(ax(1:2), [-1.5 1.5 -700 700]);
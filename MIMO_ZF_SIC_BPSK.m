% Make sure to define this function or just download the same for the code
% to run.
% Get the channel co-efficient rayleigh distributed
% E[|h|^2] = 1 ==> complex gaussian circular symmetric with variance 1
%function [H] = get_H(Nr,Nt)
%    H = sqrt(0.5)*(randn(1, Nr*Nt) + j*randn(1, Nr*Nt));
%    H = reshape(H, Nr, Nt);
%end

% For any queries or doubts or any Matlab code requirements
% please write to "sweetdahaka144@outlook.com"

% Script for calculating the BER for BPSK modulation in a
% Rayleigh fading channel with 2 Tx, 2Rx MIMO channel 
% 1)Zero forcing receiver (ZF)
% 2)Zero forcing with Successive Interference (No specific ordering in cancellation)(ZF SIC) 
% 3)Zero forcing with Successive Interference Cancellation with optimal
%   ordering (ZF SIC Optimal)
% You can change Nr, Nt, dbMax, dbMin, dbStep as you wish and the code
% still works

%Optimal Ordering means in successive interference cancellation, the symbol
%with high SNR is detected first and its effect is removed and the
%procedure continues until we are left with only one symbol. The advantage
%of this method is that the symbols which are decoded in later stages have
%an increased diversity order. 

%Variables explained 
% Nr ---> Rx Antennas
% Nt ---> Tx Antennas
% dbMax, dbMin ---> Max and Min dB of Eb/No
% dbStep ---> increment in dB in each loop
% N ---> Number of trials / Number of times experiment is performed
% data ---> Data vector 
% x ---> Symbols vector (Modulated Symbols in this case BPSK)(size-(Ntx1))
% noise ---> AWGN noise at the receiver with variance sigma^2 (size-(Nrx1))
% sigma ---> Standard deviation of noise calculated for a given Eb/No
%            considering Eb = 1 from Eb/No in dB
% H ---> Channel matrix consisting of independent rayleigh distributed co-efficients(size-Nr x Nt)
%        get_H(Nr,Nt) gives H matrix with independent rayleigh distributed
%        (Nr x Nt)matrix
% x_cap ---> Symbols detected by using ZF,ZF SIC, ZF SIC Optimal schemes
%




%clc;
tic;
Nr = 2; % Number of receive antennas
Nt = 2; % Number of transmit antennas
dbMax = 10;
dbMin = 0;
dbStep = 1;
dbLoopLength = ((dbMax-dbMin)/dbStep)+1;
N = 10000; % Number of trials

%initialisations for matlab perfomance
x = zeros(Nt, N);
data = x;
y = zeros(Nr, N);
noise = y;
x_cap_sic = zeros(1,Nt);


%Channel equation for different Eb/No
bitError_sim_zf = zeros(1, dbLoopLength);
bitError_sim_zf_sic = bitError_sim_zf;
bitError_sim_zf_sic_sort = bitError_sim_zf;
for i = 1:dbLoopLength
    EbNodB = dbMin + (i-1)*dbStep;
    sigma = 1/sqrt((10^(EbNodB/10))*2);
    Nerrs_zf = 0;
    Nerrs_zf_sic = 0;
    Nerrs_zf_sic_sort = 0;
    for trialIndex = 1:N
        
        % Information generation and modulation
        data(:,trialIndex) = randi([0 1], Nt, 1);
        x(:,trialIndex) = 1-2*data(:,trialIndex);
        
        %Channel equation
        H = get_H(Nr, Nt);
        noise(:, trialIndex) = sigma*(randn(Nr,1) + j*randn(Nr,1));
        y(:, trialIndex) = H*x(:, trialIndex) + noise(:, trialIndex);
        
        % Zero forcing receiver
        x_cap = ((H'*H)\H')*y(:, trialIndex);
        x_cap = x_cap < 0;
        Nerrs_zf = Nerrs_zf + sum(x_cap ~= data(:, trialIndex));
        
        % Zero forcing SIC receiver
        H_eff = H;
        y_eff = y(:, trialIndex);
        x_cap_sic = zeros(Nt,1);
        for zfSicIndex = 1:Nt
            x_eff = ((H_eff'*H_eff)\H_eff')*y_eff;
            
            % Decoding the first symbol
            if(x_eff(1) < 0)
                temp = -1;
            else
                temp = 1;
            end
            
            x_cap_sic(zfSicIndex) = temp;
            y_eff = y_eff - (temp*H_eff(:,1));
            H_eff(:,1) = [];
        end
        x_cap = x_cap_sic < 0;
        Nerrs_zf_sic = Nerrs_zf_sic + sum(x_cap ~= data(:, trialIndex));
        
        
        % Zero forcing SIC Sort receiver
        H_eff = H;
        y_eff = y(:, trialIndex);
        x_cap_sic = zeros(Nt,1);
        [~,sicOrder] = sort(sum(abs(H_eff).^2), 'descend');
        for zfSicIndex = 1:Nt
            garbage = sum(abs(H_eff).^2);
            [~,pos] = max(garbage(:));
            x_eff = ((H_eff'*H_eff)\H_eff')*y_eff;
            
            % Decoding the symbol with highest SNR
            if(x_eff(pos) < 0)
                temp = -1;
            else
                temp = 1;
            end
            x_cap_sic(sicOrder(zfSicIndex)) = temp;
            y_eff = y_eff - (temp*H_eff(:,pos)); % Cancelling the effect of the decoded symbol
            H_eff(:,pos) = [];
        end
        x_cap = x_cap_sic < 0;
        Nerrs_zf_sic_sort = Nerrs_zf_sic_sort + sum(x_cap ~= data(:, trialIndex));
        
    end
    
    bitError_sim_zf(i) = Nerrs_zf/(N*Nt);
    bitError_sim_zf_sic(i) = Nerrs_zf_sic/(N*Nt);
    bitError_sim_zf_sic_sort(i) = Nerrs_zf_sic_sort/(N*Nt);
end
figure
xaxis = dbMin:dbStep:dbMax;
semilogy(xaxis, bitError_sim_zf, 'bp-', xaxis, bitError_sim_zf_sic, 'kd-', xaxis, bitError_sim_zf_sic_sort, 'mo-');
legend('ZF', 'ZF SIC', 'ZF SIC Optimal');
xlabel('E_b/N_0 [dB]');
ylabel('BER');
toc;
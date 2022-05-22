% 
% Wireless Communication Systems HW4, 通訊所一年級 110064533 陳劭珩
% 
% 2.(Repeat the problem for uncorrelated Ricean fading with K = 1.) Assume 
% that there are L(L=1,2,3,4) diversity branches of uncorrelated Ricean 
% fading signals with K=1. Each branch has the same average symbol energy-
% to–noise power ratio E_s/N_0, for E_s/N_0 = 1,3,5,7,9 dB.
%   
% Simulate the QPSK bit error probability (at least to P_b=10^(-4)) for 
% the following 4 Combining techniques:
%   (a) Selective Combining
%   (b) Maximal Ratio Combining
%   (c) Equal Gain Combining
%   (d) Direct Combining (which combines all paths directly and then 
%       compensates the overall phase shift before demodulation.
% 
% (1.You may generate the fading gains via combining a Rayleigh random 
% number and a uniform random phase, or via combining two Gaussian random 
% variables (Complex Gaussian.)
% 2.For coherence detection, must equalize the phase before demodulation.)
%
% clear;
% clc;
% 
N = 1e5;               % number of simulations, say 100000 symbols
L = 4;                 % diversity branch L = 1~4
SNR = [1 3 5 7 9];     % SNR = Es/N0 = 1,3,5,7,9 dB, in our case
Es = 1;                % set symbol energy = 1 Joule 
N0 = Es*10.^(-SNR/10); % compute N0 from given SNR, Es. derived as below
% 10log(Es/N0) = SNR(dB) --> log(Es/N0) = SNR/10 --> 10^(SNR/10) = Es/N0
% --> N0 = Es*10^-(SNR/10)
%
K = 1;                 % Ricean factor K = 1
Ricean_mean = 1/2;     % Ricean mean m = 0.5
Ricean_sigma = 1/2;    % Ricean variance σ^2 = 0.25
%
%
SelectiveComb_data = zeros(1, N);
SelectiveComb_error = zeros(L, length(SNR));
SelectiveComb_BER = zeros(L, length(SNR));
%
MaxRatioComb_data = zeros(L, N);
MaxRatioComb_error = zeros(L, length(SNR));
MaxRatioComb_BER = zeros(L, length(SNR));
%
EqualGainComb_data = zeros(L, N);
EqualGainComb_error = zeros(L, length(SNR));
EqualGainComb_BER = zeros(L, length(SNR));
%
DirectComb_error = zeros(L, length(SNR));
DirectComb_BER = zeros(L, length(SNR));
%
% Ricean fading channel case
%
for N0_idx = 1 : length(N0)
    %
    % QPSK modulation
    %
    x = 2*round(rand(1, N)) - 1;
    y = 2*round(rand(1, N)) - 1;
    modulated_data = sqrt(1/2) * (x + 1i*y);
    %
    % Channel
    %
    channel_gain = Ricean_sigma .* (randn(L,N) + randn(L,N)*1i) ...
                 + Ricean_mean * (1 + 1i);
    awgn_noise = sqrt(N0(N0_idx)/2) * randn(L,N) ...
               + sqrt(N0(N0_idx)/2) * 1i*randn(L,N);
    branch_signal = channel_gain .* modulated_data;
    rx_data = branch_signal + awgn_noise;
    %
    % L = 1~4 diversity branches
    %
    for L_idx = 1 : L
        %
        % Selctive Combining
        %
        if L_idx == 1
            % equalize the phase
            SelectiveComb_data = rx_data(L_idx,:) ...
                               .* exp(-1i*angle(channel_gain(L_idx,:)));
        else
            [~, max_idx] = max(abs(channel_gain(1:L_idx, :)));
            for N_idx = 1 : N
                SelectiveComb_data(N_idx) = rx_data(max_idx(N_idx), N_idx) ...
                    .* exp(-1i*angle(channel_gain(max_idx(N_idx), N_idx)));
            end
        end
        detect_real = 2*(real(SelectiveComb_data) > 0) - 1;
        detect_imag = 2*(imag(SelectiveComb_data) > 0) - 1;
        SelectiveComb_error(L_idx, N0_idx) = length(find(detect_real - x)) ...
                                           + length(find(detect_imag - y));
        SelectiveComb_BER = SelectiveComb_error / (2*N);
        %
        % Maxmimal Ratio Combining
        %
        MaxRatioComb_data(L_idx,:) = conj(channel_gain(L_idx,:)) ...
                                   .* rx_data(L_idx,:);
        if L_idx == 1
            detect_real = 2*(real(MaxRatioComb_data(1:L_idx,:)) > 0) - 1;
            detect_imag = 2*(imag(MaxRatioComb_data(1:L_idx,:)) > 0) - 1;
        else
            detect_real = 2*(real(sum(MaxRatioComb_data(1:L_idx,:)))>0) -1;
            detect_imag = 2*(imag(sum(MaxRatioComb_data(1:L_idx,:)))>0) -1;
        end
        MaxRatioComb_error(L_idx, N0_idx) = length(find(detect_real - x))  ...
                                          + length(find(detect_imag - y));
        MaxRatioComb_BER = MaxRatioComb_error / (2*N);
        %
        % Equal Gain Combining
        %
        EqualGainComb_data(L_idx,:) = rx_data(L_idx,:) ...
                                  .* exp(-1i*angle(channel_gain(L_idx,:)));
        if L_idx==1
            detect_real=2*((real(EqualGainComb_data(1:L_idx,:)) > 0)) - 1;
            detect_imag=2*((imag(EqualGainComb_data(1:L_idx,:)) > 0)) - 1;
        else
            detect_real=2*((real(sum(EqualGainComb_data(1:L_idx,:)))>0))-1;
            detect_imag=2*((imag(sum(EqualGainComb_data(1:L_idx,:)))>0))-1;
        end
        EqualGainComb_error(L_idx, N0_idx) = length( find(detect_real - x)) ...
                                           + length(find(detect_imag - y));
        EqualGainComb_BER = EqualGainComb_error / (2*N);
        %
        % Direct Combining
        %
        if L_idx == 1
            phase_equalize = exp(-1i*angle(channel_gain(L_idx,:)));
            detect_real = (2*(real(rx_data(1:L_idx,:) ...
                        .*phase_equalize) > 0)) - 1;
            detect_imag = (2*(imag(rx_data(1:L_idx,:) ...
                        .*phase_equalize) > 0)) - 1;
        else
            phase_equalize = exp(-1i*angle(sum(channel_gain(1:L_idx,:))));
            detect_real = (2*(real(sum(rx_data(1:L_idx,:) ...
                        .*phase_equalize)) > 0)) - 1;
            detect_imag = (2*(imag(sum(rx_data(1:L_idx,:) ...
                        .*phase_equalize)) > 0)) - 1;
        end
        DirectComb_error(L_idx,N0_idx) = length(find(detect_real - x)) ...
                                       + length(find(detect_imag - y));
        DirectComb_BER = DirectComb_error / (2*N);
     end 
end
% 
% plot and save, 4 diversity techniques
% 
figure(1);
% semilogy(SNR, SelectiveComb_BER,'-*'); % plot different L separately
semilogy(SNR, SelectiveComb_BER(1, :), '-*', ...
         SNR, SelectiveComb_BER(2, :), '-o', ...
         SNR, SelectiveComb_BER(3, :), '-x', ...
         SNR, SelectiveComb_BER(4, :), '-d');
title('Selective Combining for Ricean Fading Channel');
xlim([1 9]);
ylim([10^-5 1]);
xlabel('E_s/N_0(dB)', 'FontWeight', 'bold');
ylabel('QPSK Bit Error Rate P_b', 'FontWeight', 'bold');
legend('L=1', 'L=2', 'L=3', 'L=4', 'location', 'southwest');
grid on;
%
% save as Selective Combining for Ricean Fading Channel.png
fname = 'Selective Combining for Ricean Fading Channel.png';
print(fname, '-dpng'); 
%
figure(2);
% semilogy(SNR, MaxRatioComb_BER,'-o'); % plot different L separately
semilogy(SNR, MaxRatioComb_BER(1, :), '-*', ...
         SNR, MaxRatioComb_BER(2, :), '-o', ...
         SNR, MaxRatioComb_BER(3, :), '-x', ...
         SNR, MaxRatioComb_BER(4, :), '-d');
title('Maximal Ratio Combining for Ricean Fading Channel');
xlim([1 9]);
ylim([10^-5 1]);
xlabel('E_s/N_0(dB)', 'FontWeight', 'bold');
ylabel('QPSK Bit Error Rate P_b', 'FontWeight', 'bold');
legend('L=1', 'L=2', 'L=3', 'L=4', 'location', 'southwest');
grid on;
%
% save as Maximal Ratio Combining for Ricean Fading Channel.png
fname = 'Maximal Ratio Combining for Ricean Fading Channel.png';
print(fname, '-dpng'); 
%
figure(3);
% semilogy(SNR, EqualGainComb_BER,'-x'); % plot different L separately
semilogy(SNR, EqualGainComb_BER(1, :), '-*', ...
         SNR, EqualGainComb_BER(2, :), '-o', ...
         SNR, EqualGainComb_BER(3, :), '-x', ...
         SNR, EqualGainComb_BER(4, :), '-d');
title('Equal Gain Combining for Ricean Fading Channel');
xlim([1 9]);
ylim([10^-5 1]);
xlabel('E_s/N_0(dB)', 'FontWeight', 'bold');
ylabel('QPSK Bit Error Rate P_b', 'FontWeight', 'bold');
legend('L=1', 'L=2', 'L=3', 'L=4', 'location', 'southwest');
grid on;
%
% save as Equal Gain Combining for Ricean Fading Channel.png
fname = 'Equal Gain Combining for Ricean Fading Channel.png';
print(fname, '-dpng'); 
%
figure(4);
% semilogy(SNR, DirectComb_BER,'-d'); % plot different L separately
semilogy(SNR, DirectComb_BER(1, :), '-*', ...
         SNR, DirectComb_BER(2, :), '-o', ...
         SNR, DirectComb_BER(3, :), '-x', ...
         SNR, DirectComb_BER(4, :), '-d');
title('Direct Combining for Ricean Fading Channel');
xlim([1 9]);
ylim([10^-5 1]);
xlabel('E_s/N_0(dB)', 'FontWeight', 'bold');
ylabel('QPSK Bit Error Rate P_b', 'FontWeight', 'bold');
legend('L=1', 'L=2', 'L=3', 'L=4', 'location', 'southwest');
grid on;
%
% save as Direct Combining for Ricean Fading Channel.png
fname = 'Direct Combining for Ricean Fading Channel.png';
print(fname, '-dpng'); 
%
% 
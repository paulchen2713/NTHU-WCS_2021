% 
% Wireless Communication Systems HW4, 通訊所一年級 110064533 陳劭珩
% 
% 1.Assume that there are L(L=1,2,3,4) diversity branches of uncorrelated 
% Rayleigh fading signals. Each branch has the same average symbol energy-
% to¡Vnoise power ratio E_s/N_0, for E_s/N_0 = 1,3,5,7,9 dB.
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
% signal mapping
s00 = [ 1  1];
s01 = [ 1 -1];
s11 = [-1  1];
s10 = [-1 -1];
% 
% source input, store the real parts and imaginary parts of the complex 
% input signal as source_real, and source_imag, separately.
source_real = zeros(1, N);
source_imag = zeros(1, N);
r = zeros(1, L);       % received signal
BER = zeros(L, length(SNR), L);  % Bit Error Rate
% 
% L = 1~4 diversity branches
%
for L_idx = 1 : L
    for k = 1 : length(SNR)
        num_bit_error = zeros(L, 1);
        for i = 1 : N
            % Initialize the random source signal for different SNR cases
            temp = rand; % uniformly distributed in (0, 1)
            if temp < 0.25
                % if rand < 0.25, then set the source input as '00'
                source_real(i) = 0;
                source_imag(i) = 0;
            elseif 0.25 < temp && temp < 0.5
                % if 0.25 < rand < 0.5, then set the source input as '01'
                source_real(i) = 0;
                source_imag(i) = 1;
            elseif 0.5 < temp && temp < 0.75
                % if 0.5 < rand < 0.75, then set the source input as '11'
                source_real(i) = 1;
                source_imag(i) = 0;
            elseif 0.75 < temp
                % if 0.5 < rand < 0.75, then set the source input as '10'
                source_real(i) = 1;
                source_imag(i) = 1;
            end
            %
            %
            diversity=zeros(L,1);
            gain = 0;
            %
            % Generate the received signal
            for j = 1 : L_idx
                %
                % Rayleigh = the amplitube of two complex Gaussian combined
                rayleigh = sqrt(1/2) * (randn + 1i*randn); 
                % AWGN noise
                noise = sqrt(N0(k)/2) * (randn + 1i*randn);
                %
                % received signal = sorce signal * channel fading + noise
                if ((source_real(i)==0) && (source_imag(i)==0))
                    % r = s00 * Rayleigh fading + noise
                    r(j) = rayleigh * sqrt(Es/2) * s00 * [1; 1i] + noise;
                elseif ((source_real(i)==0) && (source_imag(i)==1))
                    % r = s01 * Rayleigh fading + noise
                    r(j) = rayleigh * sqrt(Es/2) * s01 * [1; 1i] + noise;
                elseif (source_real(i) == 1) && (source_imag(i) == 0)
                    % r = s10 * Rayleigh fading + noise
                    r(j) = rayleigh * sqrt(Es/2) * s10 * [1; 1i] + noise;
                elseif (source_real(i) == 1) && (source_imag(i) == 1)
                    % r = s11 * Rayleigh fading + noise
                    r(j) = rayleigh * sqrt(Es/2) * s11 * [1; 1i] + noise;
                end
                %
                % technique (1) Selective Combining
                if abs(conj(rayleigh) * r(j)) > abs(diversity(1))
                    diversity(1) = conj(rayleigh) * r(j);
                end
                %
                % technique (2) Maximal Ratio Combining
                diversity(2) = diversity(2) + conj(rayleigh) * r(j);
                %
                % technique (3) Equal Gain Combining
                temp_equal_gain = conj(rayleigh) * r(j) / abs(rayleigh);
                diversity(3) = diversity(3) + temp_equal_gain;
                %
                % technique (4) Direct Combining 
                diversity(4) = diversity(4) + r(j);
                gain = gain + rayleigh;
            end
            diversity(4) = conj(gain) * diversity(4) / abs(gain);
            %
            for j = 1 : L
                %
                % compute the correlation matrixs
                %
                c00 = real(diversity(j))*s00(1)*sqrt(Es/2) ...
                    + imag(diversity(j))*s00(2)*sqrt(Es/2);
                c01 = real(diversity(j))*s01(1)*sqrt(Es/2) ...
                    + imag(diversity(j))*s01(2)*sqrt(Es/2);
                c10 = real(diversity(j))*s10(1)*sqrt(Es/2) ...
                    + imag(diversity(j))*s10(2)*sqrt(Es/2);
                c11 = real(diversity(j))*s11(1)*sqrt(Es/2) ...
                    + imag(diversity(j))*s11(2)*sqrt(Es/2);
                %
                % decide the received signal based on the correlations
                %
                % constellation point (received_real, received_imag)
                max_corr = max([c00 c01 c10 c11]);
                if (c00 == max_corr)     % s00
                    received_real=0; 
                    received_imag=0;
                elseif (c01 == max_corr) % s01
                    received_real=0; 
                    received_imag=1;
                elseif (c10 == max_corr) % s10
                    received_real=1; 
                    received_imag=0;
                elseif (c11 == max_corr) % s11
                    received_real=1; 
                    received_imag=1;
                end
                %
                if (received_real ~= source_real(i))
                    num_bit_error(j) = num_bit_error(j) + 1;
                end
                if (received_imag ~= source_imag(i))
                    num_bit_error(j) = num_bit_error(j) + 1;  
                end
            end
        end
        BER(:, k, L_idx) = num_bit_error / (2*N);
    end
end
%
% plot and save, 4 diversity techniques
% 
for L_idx = 1 : L 
    figure; % plot different L separately
    semilogy(SNR, BER(L_idx,:,1), '-*', ...
             SNR, BER(L_idx,:,2), '-o', ...
             SNR, BER(L_idx,:,3), '-x', ...
             SNR, BER(L_idx,:,4), '-d');
    hold on;
    legend('L=1', 'L=2', 'L=3', 'L=4', 'location', 'southwest');
    xlabel('E_s/N_0(dB)', 'FontWeight', 'bold');
    ylabel('QPSK Bit Error Rate P_b', 'FontWeight', 'bold');
    xlim([1 9]);
    ylim([10^-5 1]);
    grid on;
    if L_idx == 1
        title('Selective Combining for Rayleigh Fading Channel');
        %
        % save as Selective Combining for Rayleigh Fading Channel.png
        fname = 'Selective Combining for Rayleigh Fading Channel.png';
        print(fname, '-dpng'); 
    elseif L_idx == 2
        title('Maximal Ratio Combining for Rayleigh Fading Channel');
        %
        % save as Maximal Ratio Combining for Rayleigh Fading Channel.png
        fname = 'Maximal Ratio Combining for Rayleigh Fading Channel.png';
        print(fname, '-dpng'); 
    elseif L_idx == 3
        title('Equal Gain Combining for Rayleigh Fading Channel');
        %
        % save as Equal Gain Combining for Rayleigh Fading Channel.png
        fname = 'Equal Gain Combining for Rayleigh Fading Channel.png';
        print(fname, '-dpng'); 
    elseif L_idx == 4
        title('Direct Combining for Rayleigh Fading Channel');
        %
        % save as Direct Combining for Rayleigh Fading Channel.png
        fname = 'Direct Combining for Rayleigh Fading Channel.png';
        print(fname, '-dpng'); 
    end
end
%

% 
% Wireless Communication Systems HW3, 通訊所一年級 110064533 陳劭珩
% 
% 1. Implement a Rayleigh fading channel simulator based on the Filtered
%    Gaussian Noise method
%    – Plot the channel output for fmT = 0.01, 0.1, 0.5 (t/T = 0~300)
%    – Plot the channel output autocorrelation for fmT = 0.01, 0.1, 0.5
%      (fmτ = 0~10)
% 
clear;
clc;
% 
fmT = [0.01, 0.1, 0.5]; % fmT = 0.01, 0.1, 0.5
t_over_T = 0:1:300;     % t/T = 0~300
% 
K = 1e5;                % epoch k = 1, 2, 3, ... , K. 10^5 samples
omega_p = 2;            % total power Ω_p
power = omega_p / 2;    % power = σ_g_I^2 = σ_g_Q^2 = Ω_p/2
shift_t = round(300 + (K-300)*rand);
% 
r = zeros(1, K);        % envelope r(t) = g_I(t) + j*g_Q(t)
g_I = zeros(1, K);      % g_I(t) = LPF(Gaussian noise)
g_Q = zeros(1, K);      % g_Q(t) = LPF(Gaussian noise)
envelope = zeros(1, K);
envelope_dB = zeros(1, K);
% 
for i = 1:3
    %
    % ζ^2 - 2ζ(2 - cos(πf_mT/2)) + 1 = 0
    % ζ = 2 - cos(πf_mT/2) - √(2 - cos(πf_mT/2)^2 - 1)
    %
    zeta = 2 - cos(pi*fmT(i)/2) - sqrt((2-cos(pi*fmT(i)/2))^2 - 1);
    %
    % σ_g_I^2 = power = (1 - ζ)/(1 + ζ) * σ^2, σ^2 variance of w_1,k w_2,k
    % σ = √((1 + ζ)/(1 - ζ) * power)
    %
    sigma = sqrt((1 + zeta)/(1 - zeta) * power);
    %
    % Gaussian noise source. w_1,k and w_2,k
    %
    w_I = randn(1, K) * sigma; % w_I = w_1,k
    w_Q = randn(1, K) * sigma; % w_Q = w_2,k
    %
    % pass the Gaussian noise w_1,k and w_2,k through the first-order LPF,
    % we can obtain th real and imaginary parts of the complex envelope
    %
    g_I(1) = (1 - zeta) * w_I(1);  
    g_Q(1) = (1 - zeta) * w_Q(1);
    %
    envelope(1) = sqrt(g_I(1)^2 + g_Q(1)^2);
    envelope_dB(1) = 10 * log10(envelope(1));
    %
    for k = 2:K
        %
        g_I(k) = zeta * g_I(k-1) + (1-zeta) * w_I(k-1);
        g_Q(k) = zeta * g_Q(k-1) + (1-zeta) * w_Q(k-1);
        %
        r(k) = g_I(k) + 1i*g_Q(k);
        envelope(k) = sqrt(g_I(k)^2 + g_Q(k)^2);
        envelope_dB(k) = 10 * log10(envelope(k));
        %
    end
    %
    % plot and save
    fig = figure(1);
    if i == 1
        fig;
        subplot(3, 2, 1);
        plot1_1 = plot(t_over_T, envelope_dB(0+shift_t : 300+shift_t));
        set(plot1_1, 'LineWidth', 1, 'LineStyle', '-', 'Color', 'r');
        xlabel('Time, t/T');
        ylabel('Envelope Level(dB)');
        title('Channel output when f_mT=0.01');
        %
        autocorrelation = autocorr(r(:), 10/fmT(i));
        subplot(3, 2, 2);
        time_delay = 0 : fmT(i) : 10;
        plot1_2 = plot(time_delay, autocorrelation);
        set(plot1_2, 'LineWidth', 1, 'LineStyle', '-', 'Color', 'r');
        xlabel('Time delay, f_mτ');
        ylabel('Autocorrelation, φgg(τ)');
        title('Channel output autocorrelation when f_mT = 0.01');
    elseif i == 2
        fig;
        subplot(3, 2, 3);
        plot2_1 = plot(t_over_T, envelope_dB(0+shift_t : 300+shift_t));
        set(plot2_1, 'LineWidth', 1, 'Color', [0.9763 0.8831 0.0338]);
        xlabel('Time, t/T');
        ylabel('Envelope Level(dB)');
        title('Channel output when f_mT=0.1');
        %
        autocorrelation = autocorr(r(:), 10/fmT(i));
        subplot(3, 2, 4);
        time_delay = 0 : fmT(i) : 10;
        plot2_2 = plot(time_delay, autocorrelation);
        set(plot2_2, 'LineWidth', 1, 'Color', [0.9763 0.8831 0.0338]);
        xlabel('Time delay, f_mτ');
        ylabel('Autocorrelation, φgg(τ)');
        title('Channel output autocorrelation when f_mT = 0.1');
    elseif i == 3
        fig;
        subplot(3, 2, 5);
        plot3_1 = plot(t_over_T, envelope_dB(0+shift_t : 300+shift_t));
        set(plot3_1, 'LineWidth', 1, 'LineStyle', '-', 'Color', 'g');
        xlabel('Time, t/T');
        ylabel('Envelope Level(dB)');
        title('Channel output when f_mT=0.5');
        %
        autocorrelation = autocorr(r(:), 10/fmT(i));
        subplot(3, 2, 6);
        time_delay = 0 : fmT(i) : 10;
        plot3_2 = plot(time_delay, autocorrelation);
        set(plot3_2, 'LineWidth', 1, 'LineStyle', '-', 'Color', 'g');
        xlabel('Time delay, f_mτ');
        ylabel('Autocorrelation, φgg(τ)');
        title('Channel output autocorrelation when f_mT = 0.5');
    end
    % save the plots and the parameters, if these files are not exist
    if isfile('HW3_1_110064533.png') == 0
        print(fig, 'HW3_1_110064533.png', '-dpng');
%       save('HW3_1_parameters.mat', 'w_I', 'w_Q', 'zeta', 'r', 'g_I', ...
%               'g_Q', 'envelope', 'envelope_dB', 'autocorrelation');
    end
end
%
%





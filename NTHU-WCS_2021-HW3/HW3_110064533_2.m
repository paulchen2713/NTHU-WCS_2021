% 
% Wireless Communication Systems HW3, 通訊所一年級 110064533 陳劭珩
% 
% 2. Implement a Rayleigh fading channel simulator based on the Sum of 
%    Sinusoids method
%    – Plot the channel output for M = 8, 16 (fmT = 0.01, 0.1, 0.5 and 
%      t/T = 0~300)
%    – Plot the channel output autocorrelation for M = 8, 16 (fmτ = 0~10)
% 
fmT = [0.01, 0.1, 0.5]; % fmT = 0.01, 0.1, 0.5
t_over_T = 0:1:300;     % t/T = 0~300
M = [8, 16];
N = 2 * (2*M + 1);
% 
K = 1e5;                % 10^5 samples
fm = fmT;
shift_t = round(300 + (K-300)*rand);
%
g = zeros(1, K);        % envelope r(t) = g_I(t) + j*g_Q(t)
g_I = zeros(1, K);      % g_I(t) = LPF(Gaussian noise)
g_Q = zeros(1, K);      % g_Q(t) = LPF(Gaussian noise)
envelope = zeros(1, K);
envelope_dB = zeros(1, K);
%
for j = 1:2
    for i = 1:3
        %
        n = 1:M(j);
        theta = 2*pi*n / N(j); % uniform incoming direction incident angle
        fn = fm(i)*cos(theta);
        alpha = 0;
        beta_n = pi*n / M(j);  % uniform phase
        %
        for t = 0:K-1
            %
            g_I(t+1) = sqrt(2) * (2*sum(cos(beta_n).*cos(2*pi*fn*t)) + ...
                sqrt(2)*cos(alpha).*cos(2*pi*fm(i)*t)); 
            g_Q(t+1) = sqrt(2) * (2*sum(sin(beta_n).*cos(2*pi*fn*t)) + ...
                sqrt(2)*sin(alpha).*cos(2*pi*fm(i)*t));
            %
            g(t+1) = g_I(t+1) + 1i*g_Q(t+1);
            envelope_dB(t+1) = 10 * log10(sqrt(g_I(t+1)^2 + g_Q(t+1)^2));
        end
        %
        % plot
        fig = figure(j);
        if i == 1
            fig;
            subplot(3, 2, 1);
            plot1_1 = plot(t_over_T, envelope_dB(0+shift_t : 300+shift_t));
            set(plot1_1, 'LineWidth', 1, 'LineStyle', '-', 'Color', 'r');
            xlabel('Time, t/T');
            ylabel('Envelope Level(dB)');
            title('Channel output when f_mT=0.01');
            %
            autocorrelation = autocorr(g(:), 10/fmT(i));
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
            autocorrelation = autocorr(g(:), 10/fmT(i));
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
            autocorrelation = autocorr(g(:), 10/fmT(i));
            subplot(3, 2, 6);
            time_delay = 0 : fmT(i) : 10;
            plot3_2 = plot(time_delay, autocorrelation);
            set(plot3_2, 'LineWidth', 1, 'LineStyle', '-', 'Color', 'g');
            xlabel('Time delay, f_mτ');
            ylabel('Autocorrelation, φgg(τ)');
            title('Channel output autocorrelation when f_mT = 0.5');
        end
    end
end


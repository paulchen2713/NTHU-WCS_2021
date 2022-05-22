%
% Wireless Communication Systems Term Project, 通訊所一年級 110064533 陳劭珩
%
% ref. Chap 7, Fundamentals of Wireless Communication 2005.
%      Chap 2, Wireless Communication Systems lecture handouts
% 
% Angular-domain ULA Radiation/Reception
% 
% * According to the angular domain model, evaluate the Radiation and 
%   Reception patterns of uniform linear arrays(ULA)
% 
% * Input parameters:
%   – The number of antennas N
%   – The normalized antenna separation Δ (normalized to wavelength λ_c)
%   – SIMO(reception) or MISO(radiation)
%   – The radiation or reception directions of the desired signals
%   – The radiation or reception direction of the interference signal
% 
% * Output results:
%   – The angular domain radiation/reception basis
%   – The correlation between different basis vectors
%   – The gain pattern of the ULA
%   – The gain of desired signal for different radiation/reception beams
%   – The signal to interference power ratio(SINR) for different beams
%   – The SINR of multiple input signals(multiple reception directions) 
%     with diversity combining(considerfading for signals and interference)
%   – (Any results that can present your work better)
%
close all;
clear;
clc; 
%
% Inputs
%
N = 16;              % number of antenna
delta = 0.5;         % normalized antenna separation, [0, 0.5)
sim_case = 2;        % simulation cases, 1 for SIMO, and 2 for MISO
num_desired_dir = 4; % number of reception directions of desired signal
L = N * delta;       % length of the antenna array
devide = 10;         % radian = pi / devide, devide = [1 180]
%
% select simulation cases, 1 for SIMO(reception), and 2 for MISO(radiation)
%
if mod(sim_case, 2) == 1
    %
    % SIMO(reception) simulation case
    %
    arr_desired_dir = zeros(num_desired_dir, 1);
    for i = 1 : num_desired_dir
        %
        % randomly initialize desired signal receive direction
        %
        temp_desired = randi([1 devide]);
        arr_desired_dir(i) = pi / temp_desired;
    end
    %
    % randomly initialize interference signal receive direction
    %
    temp_interference = randi([1 devide]);
    interference_dir = pi / temp_interference;
else
    %
    % MISO(radiation) simulation case
    %
    arr_desired_dir = zeros(num_desired_dir, 1);
    for i = 1 : num_desired_dir
        %
        % randomly initialize desired signal receive direction
        %
        temp_desired = randi([1 devide]);
        arr_desired_dir(i) = pi / temp_desired;
    end
    %
    % randomly initialize interference signal receive direction
    %
    temp_interference = randi([1 devide]);
    interference_dir = pi / temp_interference;
end
%
% Outputs
%
% The angular-domain radiation/reception basis
%
U = basis(0 , N, delta);
for i = 1 : (N-1)
    U = [ U basis(i/L, N, delta)];
end
%
% The correlation between different basis vectors
%
omega_r = linspace(-2, 2, 10000);
figure;
%
% fr_function defined in Eq(7.35)
%
% f_r(Ω_r) = (1/n_r) * exp(jπ*Δ_r*Ω_r(n_r - 1)) * (sin(π*L_r*Ω_r)) 
%           * (sin(π*L_r*Ω_r / n_r)), where L_r = n_r * Δ_r
%
plot(omega_r, abs(fr_function(omega_r, N, delta, L)));
xlabel('\Omega');
ylabel('| f(\Omega) |');
title('The correlation between different basis vectors');
%
% The gain pattern of the ULA
%
phi = linspace(0, 2*pi, 10000);
%
% defined in Eq.(7.43)
%
figure;
for i = 1 : num_desired_dir
    subplot(1, num_desired_dir, i);
    %
    % fr_function defined in Eq(7.35)
    %
    % f_r(Ω_r) = (1/n_r) * exp(jπ*Δ_r*Ω_r(n_r - 1)) * (sin(π*L_r*Ω_r)) 
    %           * (sin(π*L_r*Ω_r / n_r)), where L_r = n_r * Δ_r
    %
    polar(phi, abs(fr_function(cos(phi) ...
                   - cos(arr_desired_dir(i)), N, delta, L)));
    title(['Desired signal Direction: ', ...
            num2str(arr_desired_dir(i)*180/pi), '^o']);
end
%
% The gain of desired signal for using different radiation/reception beams 
%
% spatial-domain channel gain (considering multipath exists )
%
g_desired_array = zeros(N, num_desired_dir);
g_desired = zeros(N, 1);
for i = 1 : num_desired_dir
    g_desired_array(:, i) = sqrt(10000) * (randn() + 1i*randn()) * ...
        (sqrt(N) * basis(cos(arr_desired_dir(i)) , N, delta));
    g_desired = g_desired + g_desired_array(:, i);
end
%
% angular-domain channel matrix, defined in Chap2. WCS handoouts
%
if mod(sim_case,2) == 1
    G_a_desired = U' * g_desired;
else
    G_a_desired = g_desired' * U;
end
%
gain_desired = abs(G_a_desired);
figure; 
stem(gain_desired);
title('gain of desired signal for different radiation/reception beams');
set(gca,'XTick', 1:N);
ylabel('SINR(dB)');
%
% The signal-to-interference power ratio(SINR) for using different beams
%
% let x = 1, w = 0 ---> y = g
%
g_interference = sqrt(0.5) * (randn()+1i*randn()) * ...
    (sqrt(N) * basis(cos(interference_dir) , N, delta));
%
if mod(sim_case,2) == 1
    G_a_interference = U' * g_interference;
    x_a = 1;
else
    G_a_interference = g_interference' * U;
    x_a = conj(U) * 1;
end
%
y_a_desire = G_a_desired * x_a;
y_a_interference = G_a_interference * x_a;
%
SINR = (abs(y_a_desire).^2) ./ (abs(y_a_interference).^2);
figure; 
stem(10*log10(SINR));
set(gca, 'XTick', 1:N);
title('Signal-to-Interference power Ratio(SINR) for different beams');
ylabel('SINR(dB)');
%
% The SINR of multiple input signals (multiple reception directions) with
% diversity combining (considering fading for the signals and interference)
%
% "since MISO cannot peform combining, only SIMO have combining techniques"
%
if mod(sim_case, 2) == 1
    %
    % Maximum Ratio Combining(MRC)
    %
    r_a_desire = sum(y_a_desire .* conj(G_a_desired));
    r_a_interference = sum(y_a_interference .* conj(G_a_desired));
    %
    SINR_MRC = abs(r_a_desire)^2 / abs(r_a_interference)^2;
    %
    % Equal Gain Combining(EGC)
    %
    theta_desired = angle(G_a_desired);
    r_a_desire = sum(y_a_desire .* exp(-1*i*theta_desired));
    r_a_interference = sum(y_a_interference .* exp(-1*i*theta_desired));
    %
    SINR_EGC = abs(r_a_desire)^2 / abs(r_a_interference)^2;
    %
    % Selective Combining(SC)
    %
    g_desired_max_value = 0;
    g_desired_max_index = 1;
    for i = 1 : num_desired_dir
        if g_desired_max_value < abs(g_desired_array(:,i)).^2
            g_desired_max_value = abs(g_desired_array(:,i)).^2;
            g_desired_max_index = i;
        end
    end
    g_desired_max = g_desired_array(:,i);
    %
    if mod(sim_case, 2) == 1
        G_a_desired_max = U' * g_desired_max;
    else
        G_a_desired_max = g_desired_max' * U;
    end
    r_a_desire = G_a_desired_max * x_a;
    %
    SINR_SC = (abs(sum(r_a_desire))^2) ./ (abs(sum(y_a_interference))^2);
    %
    SINR_MRC_DB = 10 * log10(SINR_MRC);
    SINR_EGC_DB = 10 * log10(SINR_EGC);
    SINR_SC_DB = 10 * log10(SINR_SC);
    %
    figure;
    stem(1:3, [SINR_SC_DB, SINR_EGC_DB, SINR_MRC_DB]);
    title('The SINR with diversity combining');
    set(gca,'XTicklabel', {});
    set(gca,'XTick', 1:3);
    xlabel('(Selective Combining,             Equal Gain Combining,          Max Ratio Combining)');
    ylabel('SINR (dB)');
end
%
%
fprintf('-----------------Input-----------------\n');
fprintf('N = %d\n', N);
fprintf('delta = %f\n', delta);
if mod(sim_case, 2) == 1
    fprintf('case 1. SIMO with combining\n');
else
    fprintf('case 2. MISO without combining\n');
end
fprintf('number of reception directions of desired signal = %d\n', ...
                                                        num_desired_dir); 
fprintf('reception direction of desired signal(in degree) = ');
for i = 1 : num_desired_dir
    fprintf('%f ', arr_desired_dir(i)*180/pi);
end
fprintf('\n');
% disp(arr_desired_dir);
fprintf('reception direction of interference signal(in degree) = %f\n', ...
                                                  interference_dir*180/pi);
fprintf('L = %d\n', L);
%
fprintf('-----------------Output-----------------\n');
fprintf('First part\n');
disp(U);
fprintf('----------------------------------------\n');
fprintf('Secnod part\n');
fprintf('\t Plot in figure 1\n');
%
fprintf('----------------------------------------\n');
fprintf('Third part\n');
fprintf('\t Plot in figure 2\n');
%
fprintf('----------------------------------------\n');
fprintf('Forth part\n');
fprintf('\t Plot in figure 3\n');
%
fprintf('----------------------------------------\n');
fprintf('Fifth part\n');
fprintf('\t Plot in figure 4 \n');
%
fprintf('----------------------------------------\n');
fprintf('Sixth part\n');
%
if mod(sim_case, 2) == 1
%     fprintf('\t SINR with MRC = %f \n', num2str(SINR_MRC_DB));
%     fprintf('\t SINR with EGC = %f \n', num2str(SINR_EGC_DB));
%     fprintf('\t SINR with SC  = %f \n', num2str(SINR_SC_DB ));
    disp(['　　 SINR with MRC = ', num2str(SINR_MRC_DB)]);
    disp(['　　 SINR with EGC = ', num2str(SINR_EGC_DB)]);
    disp(['　　 SINR with SC  = ', num2str(SINR_SC_DB )]);
else
    fprintf('MISO cannot peform combining techniques. \n');
end
%
%
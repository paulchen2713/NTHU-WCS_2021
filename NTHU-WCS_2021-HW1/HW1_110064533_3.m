% 
% Wireless Communication Systems HW1, 通訊所一年級 110064533 陳劭珩
%
% 3.Assume that there are 600 channels equally shared by 1)one, 2)two, 
%   or 3)three operators by using the frequency reuse factor N = 5.
%   - Find the maximum offered traffic load per cell for the 3 cases with
%     the blocking rate equal to 1%, 3%, 5%, or 10%
%   - Which case (one, two, or three operators) is more efficient?
% 
clear;
clc;
%
total_channels = 600;         % the number of total channels
operators = [1 2 3];          % the number of operators
N = 5;                        % the frequency reuse factor N
channel_per_cell = [0 0 0];   % the number of channels per cell of 3 cases
rho_table3 = zeros(3, 4);     % the result of maximum offered traffic load ρ
step_size = 0.0001;
%
% initialize the number of channels per cell of the 3 cases
for i = 1 : 3
    channel_per_cell(i) = total_channels / N / operators(i);
    fprintf("%d operator shared %d channels\n", i, channel_per_cell(i));
end
% 
% Blocking Rate: 1%, 3%, 5%, 10%
% 
blocking_rate = [0.01 0.03 0.05 0.1];
%
% brute force
%
for rate_index = 1 : 4
    for m = 1 : 3 % the channel number index m
        rho = 0;             % the maximum offered traffic load ρ
        blocking_prob = 0;   % the blocking prob B(ρ, m)
        while blocking_prob < blocking_rate(1, rate_index)
            rho = rho + step_size*channel_per_cell(1, m);
            numerator = 1;   % ρ^m / m!
            denominator = 1; % sigma_k=0~m {ρ^k / k!}
            for k = 1 : channel_per_cell(1, m)
                numerator = numerator * (rho / k);
                denominator = denominator + numerator;
            end
            blocking_prob = numerator / denominator;
        end
        % 當逼近的blocking_prob超過題目規定的blocking rate時會跳出while loop
        % 要的是最接近但不超過的ρ 也就是超過條件的前一個 並將結果存至ρ table3
        rho_table3(m, rate_index) = rho - step_size*channel_per_cell(1, m);
    end
end
%
% save the results as a csv file, named 'rho_table3.txt'
csvwrite('rho_table3.txt', rho_table3);



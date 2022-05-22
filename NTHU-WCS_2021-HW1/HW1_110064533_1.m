% 
% Wireless Communication Systems HW1, 通訊所一年級 110064533 陳劭珩
%
% 1.calculate the total offered traffic load ρ 
% 
clear;
clc;
%
rho_table = zeros(220, 4); % the result of total offered traffic load ρ
step_size = 0.0001; 
% 
% Blocking Rate: 1%, 3%, 5%, 10%
% 
blocking_rate = [0.01 0.03 0.05 0.1];
%
% brute force
%
for rate_index = 1 : 4
    % 
    % Channel Nuber: m = 1~20
    %
    for channel_number = 1 : 20
        m = channel_number;  % the current channel number m
        rho = 0;             % the total offered traffic load ρ
        blocking_prob = 0;   % the blocking prob B(ρ, m)
        while blocking_prob < blocking_rate(1, rate_index)
            rho = rho + step_size*m;
            numerator = 1;   % ρ^m / m!
            denominator = 1; % sigma_k=0~m {ρ^k / k!}
            for k = 1 : m
                numerator = numerator * (rho / k);
                denominator = denominator + numerator;
            end
            blocking_prob = numerator / denominator;
        end
        % 當逼近的blocking_prob超過題目規定的blocking rate時會跳出while loop
        % 想要的是最接近但不超過的ρ 也就是超過條件的前一個ρ
        rho_table(m, rate_index) = rho - step_size*m; % 並將結果存至ρ table
    end
    % 
    % Channel Nuber: m = 200~220
    %
    for channel_number = 200 : 220
        m = channel_number;  % the current channel number m
        rho = 0;             % the total offered traffic load ρ
        blocking_prob = 0;   % the blocking prob B(ρ, m)
        while blocking_prob < blocking_rate(1, rate_index)
            rho = rho + step_size*m;
            numerator = 1;   % ρ^m / m!
            denominator = 1; % sigma_k=0~m {ρ^k / k!}
            for k = 1 : m
                numerator = numerator * (rho / k);
                denominator = denominator + numerator;
            end
            blocking_prob = numerator / denominator;
        end
        % 當逼近的blocking_prob超過題目規定的blocking rate時會跳出while loop
        % 想要的是最接近但不超過的ρ 所以減去1單位channel_number的step_size
        rho_table(m, rate_index) = rho - step_size*m; % 並將結果存至ρ table
    end
end
%
% save the result as a csv file, named 'rho_table.txt'
csvwrite('rho_table.txt', rho_table);


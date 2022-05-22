% 
% Wireless Communication Systems HW1, �q�T�Ҥ@�~�� 110064533 ���oҲ
%
% 1.calculate the total offered traffic load �l 
% 
clear;
clc;
%
rho_table = zeros(220, 4); % the result of total offered traffic load �l
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
        rho = 0;             % the total offered traffic load �l
        blocking_prob = 0;   % the blocking prob B(�l, m)
        while blocking_prob < blocking_rate(1, rate_index)
            rho = rho + step_size*m;
            numerator = 1;   % �l^m / m!
            denominator = 1; % sigma_k=0~m {�l^k / k!}
            for k = 1 : m
                numerator = numerator * (rho / k);
                denominator = denominator + numerator;
            end
            blocking_prob = numerator / denominator;
        end
        % ��G��blocking_prob�W�L�D�سW�w��blocking rate�ɷ|���Xwhile loop
        % �Q�n���O�̱�������W�L���l �]�N�O�W�L���󪺫e�@�ӣl
        rho_table(m, rate_index) = rho - step_size*m; % �ñN���G�s�ܣl table
    end
    % 
    % Channel Nuber: m = 200~220
    %
    for channel_number = 200 : 220
        m = channel_number;  % the current channel number m
        rho = 0;             % the total offered traffic load �l
        blocking_prob = 0;   % the blocking prob B(�l, m)
        while blocking_prob < blocking_rate(1, rate_index)
            rho = rho + step_size*m;
            numerator = 1;   % �l^m / m!
            denominator = 1; % sigma_k=0~m {�l^k / k!}
            for k = 1 : m
                numerator = numerator * (rho / k);
                denominator = denominator + numerator;
            end
            blocking_prob = numerator / denominator;
        end
        % ��G��blocking_prob�W�L�D�سW�w��blocking rate�ɷ|���Xwhile loop
        % �Q�n���O�̱�������W�L���l �ҥH��h1���channel_number��step_size
        rho_table(m, rate_index) = rho - step_size*m; % �ñN���G�s�ܣl table
    end
end
%
% save the result as a csv file, named 'rho_table.txt'
csvwrite('rho_table.txt', rho_table);


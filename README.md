# -
졸업작품
close all
clear
clc

Nt = 64;
Nr = 2;
Nu = 4;
Nus = Nr;
Nts = Nu * Nus;
Ray_number = 5;
lambda = 1;
d = lambda / 2;
k = 2 * pi / lambda;
Niter = 100;
Null_space_dim = Nt - (Ray_number * (Nu - 1));

SNRdB = 0 : 2 : 20;
L_SNR = length(SNRdB);
SNR = 10.^(SNRdB/10);
C = nchoosek(1:Null_space_dim,Nr);
search = Nt-Nr+1-Ray_number*(Nu-1);
% for i8 =1: search
%     C1(i8 , :) = [i8 i8+1 i8+2 i8+3];
% end   
for i8 =1: search
    C1(i8 , :) = [i8 i8+1];
end  

for i0 = 1 : Niter
    
    [H, Tx_UPA, Rx_UPA, alpha] = mmWave_channel_realization_USER_ULA(Ray_number, d, k, Nt, Nr, Nu);
    Tx_UPA_H = Tx_UPA';
    H_USER = zeros(Nr,Nt,Nu);
    Tx_UPA_USER = zeros( Ray_number,Nt,Nu);
    
    for user = 1 : Nu
        user_range_H = ((user - 1) * Nr) + 1 : user * Nr;
        user_range_Tx_UPA = ((user - 1) * Ray_number) + 1 : user * Ray_number;
        Tx_UPA_USER(:,:,user) = Tx_UPA_H(user_range_Tx_UPA, :);
        H_USER(:,:,user) = H(user_range_H, :);
        
        channel_gain = abs(alpha(:,user));
        [sorted_gain channel_index] = sort(channel_gain, 'descend');
        USER_channel_index(:,user) = channel_index;
    end
    
    for user = 1 : Nu
        Tx_UPA_temp = Tx_UPA_H;
        Tx_UPA_temp(((user - 1) * Ray_number) + 1 : user * Ray_number,:) = [];
        [U,S,V] = svd(Tx_UPA_temp);
         Null_space(:,((user - 1) * Null_space_dim) + 1 : user * Null_space_dim) = V(:,(Nt - Null_space_dim) + 1 :  Nt);
         Null_space_USER(:,1:Null_space_dim) = Null_space(:,((user - 1) * Null_space_dim) + 1 : user * Null_space_dim);
        
        % group search
        window_search = zeros(Nt,Nr,length(C1));
        for SNR_idx = 1:length(SNR)
            for j = 1 : length(C1)
                for i = 1 : Nr
                    window_search(:,i,j) = Null_space_USER(:, C1(j,i));
                end
                window_search_BD(j,SNR_idx) =log2(det(eye(Ray_number)+(((Tx_UPA_USER(:,:,user) * window_search(:,:,j)) * (Tx_UPA_USER(:,:,user) * window_search(:,:,j))') / (Nts / SNR(SNR_idx) * eye(Ray_number)))));
            end
            [max_window_Value, max_window_Index] = max(window_search_BD(:,SNR_idx));
            window_select(:,((user - 1) * Nr) + 1 : user * Nr,SNR_idx) =  Null_space_USER(:,C1(max_window_Index,:));
        end
       
        % full search
        full_search_gain = zeros(length(C),1);
        full_search = zeros(Nt,Nr,length(C));
        for SNR_idx = 1:length(SNR)
            for j = 1 : length(C)
                for i = 1 : Nr
                    full_search(:,i,j) = Null_space_USER(:, C(j,i));
                end
                full_search_BD(j,SNR_idx) =log2(det(eye(Nr)+(((H_USER(:,:,user) * full_search(:,:,j)) * (H_USER(:,:,user) * full_search(:,:,j))') / (Nts / SNR(SNR_idx) * eye(Nr)))));
            end
            [maxValue, maxIndex] = max(full_search_BD(:,SNR_idx));
            full_search_select(:,((user - 1) * Nr) + 1 : user * Nr,SNR_idx) = Null_space_USER(:, C(maxIndex,:));
        end
       
        % random select
        randomValue = randperm(Null_space_dim,Nr); % 1~Null_sapce_dim
        random_select(:,((user - 1) * Nr) + 1 : user * Nr) = Null_space_USER(:, randomValue);
    end
    
    for SNR_idx = 1:length(SNR)   
        
            % random select sum rate
            interference_random_select_BD_not_RF  = zeros(Nr,Nr * Nu);
            for user_i = 1: Nu 
                    user_range = ((user_i - 1) * Nr) + 1 : user_i * Nr;
                    random_select_temp = random_select;
                    random_select_temp( : , user_range) = [];
                    for user_sum_i = 1 : Nu - 1
                        user_sum_range = ((user_sum_i - 1) * Nr) + 1 : user_sum_i * Nr;
                        interference_random_select_BD_not_RF(:,user_range) = interference_random_select_BD_not_RF(:,user_range) + (H_USER(:,:,user_i) * random_select_temp(:,user_sum_range)) * (H_USER(:,:,user_i) * random_select_temp(:,user_sum_range))';
                    end
            end
            for user_i = 1: Nu
                user_range = ((user_i - 1) * Nr) + 1 : user_i * Nr;
                C_BD_not_RF_random_temp(user_i) =  log2(det(eye(Nr)+(((H_USER(:,:,user_i) * random_select(:,user_range)) * (H_USER(:,:,user_i) * random_select(:,user_range))') / (interference_random_select_BD_not_RF(:,user_range) + Nts / SNR(SNR_idx) * eye(Nr)))));
            end
            C_BD_random_not_RF(i0, SNR_idx) = sum(C_BD_not_RF_random_temp);
      
            % window search sum rate
            interference_window_select_BD  = zeros(Nr,Nr * Nu);% sum Interference
            for user_i = 1: Nu 
                    user_range = ((user_i - 1) * Nr) + 1 : user_i * Nr;
                    window_select_temp = window_select(:,:,SNR_idx);
                    window_select_temp( : , user_range) = [];
                    for user_sum_i = 1 : Nu - 1
                        user_sum_range = ((user_sum_i - 1) * Nr) + 1 : user_sum_i * Nr;
                        interference_window_select_BD(:,user_range) = interference_window_select_BD(:,user_range) + (H_USER(:,:,user_i) * window_select_temp(:,user_sum_range)) * (H_USER(:,:,user_i) * window_select_temp(:,user_sum_range))';
                    end
            end
            for user_i = 1: Nu
                user_range = ((user_i - 1) * Nr) + 1 : user_i * Nr;
                C_BD_window_temp(user_i) =  log2(det(eye(Nr)+(((H_USER(:,:,user_i) * window_select(:,user_range)) * (H_USER(:,:,user_i) * window_select(:,user_range))') / (interference_window_select_BD(:,user_range) + Nts / SNR(SNR_idx) * eye(Nr)))));
            end
            C_BD_window(i0, SNR_idx) = sum(C_BD_window_temp);
           
            % full search sum rate
            interference_full_search_select_BD_not_RF  = zeros(Nr,Nr * Nu);% sum Interference
            for user_i = 1: Nu 
                    user_range = ((user_i - 1) * Nr) + 1 : user_i * Nr;
                    full_search_select_temp = full_search_select(:,:,SNR_idx);
                    full_search_select_temp( : , user_range) = [];
                    for user_sum_i = 1 : Nu - 1
                        user_sum_range = ((user_sum_i - 1) * Nr) + 1 : user_sum_i * Nr;
                        interference_full_search_select_BD_not_RF(:,user_range) = interference_full_search_select_BD_not_RF(:,user_range) + (H_USER(:,:,user_i) * full_search_select_temp(:,user_sum_range)) * (H_USER(:,:,user_i) * full_search_select_temp(:,user_sum_range))';
                    end
            end
            for user_i = 1: Nu
                user_range = ((user_i - 1) * Nr) + 1 : user_i * Nr;
                C_BD_full_search_not_RF_temp(user_i) =  log2(det(eye(Nr)+(((H_USER(:,:,user_i) * full_search_select(:,user_range,SNR_idx)) * (H_USER(:,:,user_i) * full_search_select(:,user_range,SNR_idx))') / (interference_full_search_select_BD_not_RF(:,user_range) + Nts / SNR(SNR_idx) * eye(Nr)))));
            end
            C_BD_full_search_not_RF(i0, SNR_idx) = sum(C_BD_full_search_not_RF_temp);
    enda
    
    if mod(i0, 1) == 0
        fprintf('%f \n',i0 * 100 / Niter); % Process percentage
    end
    
end

full_not_RF_sum_rate = abs(sum(C_BD_full_search_not_RF,1)) ./ Niter;
window_not_RF_sum_rate = abs(sum(C_BD_window,1)) ./ Niter;
random_not_RF_sum_rate = abs(sum(C_BD_random_not_RF,1)) ./ Niter;

figure();
plot(SNRdB, full_not_RF_sum_rate, '-or');
hold on 
plot(SNRdB, window_not_RF_sum_rate, '-sb');
plot(SNRdB, random_not_RF_sum_rate, '-*k');
grid on;
legend(["AoD Matrix-based BD with Exhaustive Search","AoD Matrix-based BD with Group Search", "AoD Matrix-based BD without Null Space Vector Selection",], 'Location','northwest');
title_str = sprintf('[N_T=%d, N_R=%d, USER=%d, Ray number=%d]', Nt, Nr, Nu, Ray_number);
title(title_str);
xlabel('SNR [dB]');
ylabel('Sum Rate [bits/s/Hz]');

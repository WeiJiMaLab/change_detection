function get_subj_LL_cluster_par(model_num,is_color,run_num,subj_num)

fprintf('\nGetting subject %g LL for model %g ...',subj_num,model_num);

P_C_HAT_temp = [];

if (is_color == 0) || (is_color == 2)
    color_char = '';
else
    color_char = 'c';
end

if (model_num == 6) || (model_num == 7)
    for ii = 1:4
        load(['P_C_HAT/' color_char 'Bin_P_C_HAT_1_' num2str(model_num) '_' num2str(ii) '.mat']);
        P_C_HAT_temp = cat(1,P_C_HAT_temp,P_C_HAT);
    end
    P_C_HAT = P_C_HAT_temp;
else
    
    for ii = 1:4
        load(['P_C_HAT/' color_char 'Bin_P_C_HAT_1_' num2str(model_num) '_' num2str(ii) '.mat']);
        P_C_HAT_temp = cat(ndims(P_C_HAT),P_C_HAT_temp,P_C_HAT);
    end
    
    P_C_HAT = P_C_HAT_temp;
end


clear P_C_HAT_temp

P_C_HAT = permute(P_C_HAT,[ndims(P_C_HAT) 1:(ndims(P_C_HAT)-1)]);
P_C_HAT_YES_log = P_C_HAT;
P_C_HAT_YES_log(P_C_HAT_YES_log==0) = 1/(10*samples);
P_C_HAT_YES_log = log(P_C_HAT_YES_log);

P_C_HAT_NO_log = 1-P_C_HAT;
P_C_HAT_NO_log(P_C_HAT_NO_log==0) = 1/(10*samples);
P_C_HAT_NO_log = log(P_C_HAT_NO_log);

P_C_HAT_size = size(P_C_HAT);

run_prefix = num2str(run_num);

switch is_color
    case 0
        load ../data/subjCell;
        dataCell = subjCell;
        dataCell = revert_data(dataCell);

        prefix = '';
    case 1
        load ../data/subjColorCell;
        dataCell = subjColorCell;
        dataCell = revert_data(dataCell);

        prefix = 'Color_';
    case 2
        load(['./subjFakeCell' run_prefix]);
        dataCell = subjFakeCell;
        prefix = ['Fake_' run_prefix];
    case 3
        load(['./csubjFakeCell' run_prefix]);
        dataCell = subjFakeCell;
        prefix = ['cFake_' run_prefix];
    case 4
        load(['./sasubjFakeCell' run_prefix]);
        dataCell = subjFakeCell;
        prefix = ['saFake_' run_prefix];
    case 5
        load(['./csasubjFakeCell' run_prefix]);
        dataCell = subjFakeCell;
        prefix = ['csaFake_' run_prefix];
    case 6
        load ./tsubjFakeCell;
        dataCell = subjFakeCell;
        prefix = 'tFake_';
end

% LL = zeros(P_C_HAT_size(2:end));
% P_C_HAT_full = zeros([size(dataCell{1},1) size(LL)]);
num_delta = size(P_C_HAT,1)/4;

% for each subject
% for subj_num = 1:length(dataCell)
LL = zeros(P_C_HAT_size(2:end));
Data = dataCell{subj_num};
phi = pi/90*Data(:,56:(55+8));
theta = pi/90*Data(:,64:(63+8));

% get delta as an index
delta = ((.001*round(1000*(90/pi)*sum(acos(cos(theta-phi)),2))));
curr_N = Data(:,5);
SetSizes = unique(curr_N);
N_idx =zeros(size(curr_N));
for i = 1:length(curr_N)
    N_idx(i) = find(curr_N(i)==SetSizes);
end
curr_C_hat = Data(:,2);

% convert to indices
idx = (N_idx-1)*num_delta+delta+1;

% if it is set size non-parametric
if sum(model_num == [1 6 7])~= 0
    
    % for each trial
    for ii = 1:length(idx)
%         p_c_idx = idx(ii);
        if curr_C_hat(ii)
            LL(N_idx(ii),:,:,:,:) = LL(N_idx(ii),:,:,:,:) + shiftdim(P_C_HAT_YES_log(delta(ii)+1,N_idx(ii),:,:,:,:),1);
        else
            LL(N_idx(ii),:,:,:,:) = LL(N_idx(ii),:,:,:,:) + shiftdim(P_C_HAT_NO_log(delta(ii)+1,N_idx(ii),:,:,:,:),1);
        end
    end
else
    % for each trial
    for ii = 1:length(idx)
        p_c_idx = idx(ii);
        if curr_C_hat(ii)
            LL = LL + shiftdim(P_C_HAT_YES_log(p_c_idx,:,:,:,:,:));
        else
            LL = LL + shiftdim(P_C_HAT_NO_log(p_c_idx,:,:,:,:,:));
        end
    end
end

% save the LL
save(['LL/' prefix 'LL_' num2str(subj_num) '_' num2str(model_num)],'LL','runtime','samples','time_completed','samples');

% get max_LL vals
max_P_C_HAT = Data;
if model_num == 1 % NP ER
    temp_sum = sum(max(LL,[],2),1);
    [x I] = max(temp_sum(:));
    [B] = ind2sub(size(squeeze(temp_sum)),I);
    for n = 1:size(LL,1)
        [x I] = max(squeeze(LL(n,:,B)));
        A(n) = ind2sub(size(squeeze(LL(n,:,B))),I);
    end
    
    for i = 1:length(SetSizes)
        max_P_C_HAT(Data(:,5)==SetSizes(i),2)=P_C_HAT_full(Data(:,5)==SetSizes(i),i,A(i),B);
    end
elseif model_num == 6 % NP VR
    maxLLH = -Inf;
    for ii=1:size(LL,4)
        for jj=1:size(LL,3)
            for kk=1:4
                LLH_tmp = squeeze(LL(kk,:,jj,ii));
                maxLLH_tmp(kk) = max(LLH_tmp);
                j_idx(kk)  = find(LLH_tmp==max(LLH_tmp));
            end
            if sum(maxLLH_tmp)>maxLLH
                maxLLH = sum(maxLLH_tmp);
                A = (j_idx);
                B = (jj);
                C = (ii);
            end
        end
    end
    % for each stimulus
    size(P_C_HAT)
    for ii = 1:length(idx)
        max_P_C_HAT(ii,2)=P_C_HAT(delta(ii)+1,N_idx(ii),A(N_idx(ii)),B,C);
    end
elseif model_num == 7 % NPER w/ NP guessing
     
     % FIX PRIOR TO 0.5
    prior_mid_idx = ceil(size(LL,ndims(LL))/2);
     for ii = 1:size(LL,1)
        LL_temp = squeeze(LL(ii,:,:,prior_mid_idx));
        [x I] = max(LL_temp(:));
        [a b] = ind2sub(size(LL_temp),I);
        
        A(ii) = a;
        B(ii) = b;
        
     end
     C = prior_mid_idx;
    D = [];
     % for each stimulus
    for ii = 1:length(idx)
        max_P_C_HAT(ii,2)=P_C_HAT(delta(ii)+1,N_idx(ii),A(N_idx(ii)),B(N_idx(ii)),C);
    end
    
%     temp_sum = sum(max(max(LL,[],3),[],2),1);
%     [x I] = max(temp_sum(:));
%     [C] = ind2sub(size(squeeze(temp_sum)),I);
%     C = ceil(size(LL,4)/2);
%     for n = 1:size(LL,1)
%         tmp = squeeze(LL(n,:,:,C));
%         [x I] = max(tmp(:));
%         [A(n) B(n)] = ind2sub(size(squeeze(LL(n,:,:,C))),I);
%     end
%     for i = 1:length(SetSizes)
%         max_P_C_HAT(Data(:,5)==SetSizes(i),2)=P_C_HAT(Data(:,5)==SetSizes(i),i,A(i),B(i),C);
%     end
    
elseif model_num == 17 % SR
    
    [x I] = max(LL(:));
    [A B C D] = ind2sub(size(LL),I);
    for jj = 1:length(idx)
        max_P_C_HAT(jj,2) = P_C_HAT(idx(jj),A,B,C,D);
    end
    
else
    
    [x I] = max(LL(:));
    [A B C D] = ind2sub(size(LL),I);
    for jj = 1:length(idx)
        max_P_C_HAT(jj,2) = P_C_HAT(idx(jj),A,B,C,D);
    end
end

% save the max_P_C_HAT
save(['max_P_C_HAT/' prefix 'max_P_C_HAT_' num2str(subj_num) '_' num2str(model_num)],'max_P_C_HAT','runtime','samples','time_completed','samples');

fprintf('\nDone!\n')

% end

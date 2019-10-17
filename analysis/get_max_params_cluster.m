function [A B C D] = get_max_params_cluster(model_num,is_color,run_num)
% same as get_max_params, but accepts a run number. Used when generating
% fake data for analysis testing or apparent guessing rate analysis

run_prefix = num2str(run_num);

if is_color == 0
    prefix = '';
    load subjCell
    num_subj = length(subjCell);
elseif is_color == 1
    prefix = 'Color_';
    load subjColorCell
    num_subj = length(subjColorCell);
elseif is_color == 2
    prefix = ['Fake_' run_prefix];
    load(['subjFakeCell' run_prefix]);
    num_subj = length(subjFakeCell);
elseif is_color == 3
    prefix = ['cFake_' run_prefix];
    load(['csubjFakeCell' run_prefix]);
    num_subj = length(subjFakeCell);
elseif is_color == 4
    prefix = ['saFake_' run_prefix];
    load(['sasubjFakeCell' run_prefix]);
    num_subj = length(subjFakeCell);
elseif is_color == 5
    prefix = ['csaFake_' run_prefix];
    load(['csasubjFakeCell' run_prefix]);
    num_subj = length(subjFakeCell);
elseif is_color == 6
    prefix = 'tFake_';
    load tsubjFakeCell;
    num_subj = length(subjFakeCell);
end

for i = 1:num_subj
    
    load(['LL/' prefix 'LL_' num2str(i) '_' num2str(model_num)])
    
    % if NP
    if model_num == 6
        % look for max LLH using non-parametric J, shared prior, shared theta
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
                    A(i,:) = (j_idx);
                    B(i) = (jj);
                    C(i) = (ii);
                end
            end
        end
        D = [];
    elseif model_num == 1
        temp_sum = sum(max(LL,[],2),1);
        [x I] = max(temp_sum(:));
        [B(i)] = ind2sub(size(squeeze(temp_sum)),I);
        for n = 1:size(LL,1)
            [x I] = max(squeeze(LL(n,:,B(i))));
            A(i,n) = ind2sub(size(squeeze(LL(n,:,B(i)))),I);
        end
        C = [];
        D = [];
    elseif model_num == 7
        
%                 temp_sum = squeeze(sum(max(max(LL,[],3),[],2),1));
%         [x I] = max(temp_sum(:));
%         [C(i)] = ind2sub(size(squeeze(temp_sum)),I);
%         C(i) = ceil(size(LL,4)/2);
%         for n = 1:size(LL,1)
%             L_temp = squeeze(LL(n,:,:,C(i)));
%             [x I] = max(L_temp(:));
%             [A(i,n) B(i,n)] = ind2sub(size(squeeze(LL(n,:,:,C(i)))),I);
%         end
%         
                C(i) = ceil(size(LL,4)/2);
                for n = 1:size(LL,1)
                    L_temp = squeeze(LL(n,:,:,C(i)));
                    [x I] = max(L_temp(:));
                   [A(i,n) B(i,n)] = ind2sub(size(L_temp),I);
                end
                D = [];
    else
        [x I] = max(LL(:));
        [A(i) B(i) C(i) D(i)] = ind2sub(size(LL),I);
    end
end

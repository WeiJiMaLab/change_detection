% input a row vector of model numbers and get out BMC for each subject for
% that model in the rows

function [LL_model] = compute_BMC(model_vec,is_color)

% load parameter ranges
load parameter_ranges_new

% compute punishment terms (log(1/range))
mean_range = log(abs(max(meanvec) - min(meanvec)));
mean_rangeVR = log(abs(max(meanvecVR) - min(meanvecVR)));
meanvecEG = meanvec;
mean_rangeEG = log(abs(max(meanvecEG) - min(meanvecEG)));
mean_range_NPER = log(abs(max(meanvecNPER) - min(meanvecNPER)));
theta_range = log(abs(max(thetavec) - min(thetavec)));
k_range = log(abs(max(kvec) - min(kvec)));
lapse_range = log(abs(max(lapsevec) - min(lapsevec)));
prior_range = log(abs(max(priorvec) - min(priorvec)));
q_range = log(abs(max(qvec) - min(qvec)));
slotJ_range = log(abs(max(slotJvec)-min(slotJvec)));
alpha_range = log(abs(max(alphavec) - min(alphavec)));

if size(model_vec,1)~=1
    model_vec = model_vec';
end

% select correct kind of data
if is_color == 1
    prefix = 'Color_LL_';
    load ../data/subjColorCell;
    num_subj = length(subjColorCell);
elseif is_color == 2
    load ./subjFakeCell;
    prefix= 'Fake_LL_';
    num_subj = length(subjFakeCell);
elseif is_color == 4
    load ./sasubjFakeCell;
    prefix= 'saFake_LL_';
    num_subj = length(subjFakeCell);
else
    prefix = 'LL_';
    load ../data/subjCell;
    num_subj = length(subjCell);
end

LL_model = zeros(num_subj,length(model_vec));

% for each model
for model_idx = 1:length(model_vec)
    curr_model = model_vec(model_idx);
    LL_Cell = load_model_data(curr_model,num_subj,prefix);
    
    for j = 1:length(LL_Cell)
        
        if ~mod(j,10)
            fprintf('\n%g',j);
        end
        % get current LL
        curr_LL = LL_Cell{j};
        curr_LL_max = max(curr_LL(:));
        
        % exponentiate LL after subtracting max
        curr_LL_exp = exp(curr_LL - curr_LL_max);
        
        switch curr_model
            
            case 1
                % NP ER -- first need to build 5-dimensional matrix
                all_LLH = single(zeros(length(meanvec),length(meanvec),length(meanvec),length(meanvec),length(priorvec)));
                for ii=1:length(meanvec)
                    for jj=1:length(meanvec)
                        for kk=1:length(meanvec)
                            for ll=1:length(meanvec)
                                all_LLH(ii,jj,kk,ll,:) = curr_LL(1,ii,:) + curr_LL(2,jj,:) + curr_LL(3,kk,:) + curr_LL(4,ll,:);
                            end
                        end
                    end
                end
                LL_model(j,model_idx) = log(trapz(priorvec,trapz(meanvecNPER,trapz(meanvecNPER,trapz(meanvecNPER,trapz(meanvecNPER,exp(all_LLH - max(all_LLH(:))))))))) + ...
                    max(all_LLH(:)) - 4*mean_range_NPER-prior_range;
                fprintf('Subj %g completed\n',j);
            case 2 % ER
                LL_model(j,model_idx) = log(trapz(priorvec,trapz(meanvec,curr_LL_exp))) + ...
                    curr_LL_max - mean_range-prior_range;
            case 3 % VR
                LL_model(j,model_idx) = log(trapz(priorvec,trapz(thetavec,trapz(meanvecVR,curr_LL_exp)))) + ...
                    curr_LL_max - mean_rangeVR-prior_range-theta_range;
            case 4 % SA
                k_range_temp = 1:size(curr_LL,2);
                LL_model(j,model_idx) = log(trapz(priorvec,trapz(k_range_temp,trapz(slotJvec,curr_LL_exp)))) + ...
                    curr_LL_max - slotJ_range-prior_range-log(abs(k_range_temp(end)-1));
            case 5 % IL
                k_range_temp = 1:size(curr_LL,1);
                LL_model(j,model_idx) = log(trapz(qvec,trapz(lapsevec,trapz(k_range_temp,curr_LL_exp)))) + ...
                    curr_LL_max - lapse_range-q_range-log(abs(k_range_temp(end)-1));
            case 6
                % NP VR -- first need to build 6-dimensional matrix
                all_LLH = single(zeros(length(meanvecVR),length(meanvecVR),length(meanvecVR),length(meanvecVR),length(thetavec),length(priorvec)));
                for ii=1:length(meanvecVR)
                    ii
                    for jj=1:length(meanvecVR)
                        for kk=1:length(meanvecVR)
                            for ll=1:length(meanvecVR)
                                all_LLH(ii,jj,kk,ll,:,:) = curr_LL(1,ii,:,:) + curr_LL(2,jj,:,:) + curr_LL(3,kk,:,:) + curr_LL(4,ll,:,:);
                            end
                        end
                    end
                end
                [x I] = max(all_LLH(:));
                [a(j) b(j) c(j) d(j) f(j) g(j)] = ind2sub(size(all_LLH),I);
                LL_model(j,model_idx) = log(trapz(priorvec,trapz(thetavec,trapz(meanvecVR,trapz(meanvecVR,trapz(meanvecVR,trapz(meanvecVR,exp(all_LLH - max(all_LLH(:)))))))))) + ...
                    max(all_LLH(:)) - 4*mean_rangeVR-prior_range - theta_range;
                fprintf('Subj %g completed\n',j);
            case 7 % EG = ER with with guessing
                % ER -- first need to build 6-dimensional matrix
                all_LLH = single(zeros(length(meanvecEG),length(lapsevec),length(lapsevec),length(lapsevec),length(lapsevec),length(priorvec)));
                for ii=1:length(lapsevec)
                    for jj=1:length(lapsevec)
                        for kk=1:length(lapsevec)
                            for ll=1:length(lapsevec)
                                all_LLH(:,ii,jj,kk,ll,:) = squeeze(curr_LL(1,:,ii,:) + curr_LL(2,:,jj,:) + curr_LL(3,:,kk,:) + curr_LL(4,:,ll,:));
                            end
                        end
                    end
                end
                LL_model(j,model_idx) = log(trapz(priorvec,trapz(lapsevec,trapz(lapsevec,trapz(lapsevec,trapz(lapsevec,trapz(meanvecEG,exp(all_LLH - max(all_LLH(:)))))))))) + ...
                    max(all_LLH(:)) -mean_rangeEG - 4*lapse_range -prior_range;
            case 8 % SU = slots with uniform K
                k_range_temp = 1:size(curr_LL,2);
                LL_model(j,model_idx) = log(trapz(priorvec,trapz(k_range_temp,trapz(slotJvec,curr_LL_exp)))) + ...
                    curr_LL_max - slotJ_range-prior_range-log(abs(k_range_temp(end)-1));
            case 9 % ER + single K
                
                k_range_temp = 1:size(curr_LL,2);
                LL_model(j,model_idx) = log(trapz(priorvec,trapz(k_range_temp,trapz(meanvec,curr_LL_exp)))) + ...
                    curr_LL_max - mean_range-prior_range-log(abs(k_range_temp(end)-1));
                
            case 10 % ER plus uniform K
                k_range_temp = 1:size(curr_LL,2);
                LL_model(j,model_idx) = log(trapz(priorvec,trapz(k_range_temp,trapz(meanvec,curr_LL_exp)))) + ...
                    curr_LL_max - mean_range-prior_range-log(abs(k_range_temp(end)-1));
                
            case 12 % VR plus single K
                
                k_range_temp = 1:size(curr_LL,3);
                LL_model(j,model_idx) = log(trapz(priorvec,trapz(k_range_temp,trapz(thetavec,trapz(meanvecVR,curr_LL_exp))))) + ...
                    curr_LL_max - mean_rangeVR-prior_range-theta_range-log(abs(k_range_temp(end)-1));
            case 15 % Power Law VP
                LL_model(j,model_idx) = log(trapz(priorvec,trapz(alphavec,trapz(thetavec,trapz(meanvecVR,curr_LL_exp))))) + ...
                    curr_LL_max - mean_rangeVR-prior_range-theta_range-alpha_range;
            case 16
                LL_model(j,model_idx) = log(trapz(priorvec,trapz(alphavec,trapz(meanvec,curr_LL_exp)))) + ...
                    curr_LL_max - mean_range-prior_range-alpha_range;
            case 17
                LL_model(j,model_idx) = log(trapz(priorvec,trapz(kvec,trapz(meanvec,curr_LL_exp)))) + ...
                    curr_LL_max - mean_range-prior_range-k_range;
                
        end
    end
end

end

% load model data for a given subject and model and data type
function [LL_Cell] = load_model_data(model_num,num_subj,prefix)

LL_Cell = cell(1,num_subj);

for i = 1:num_subj
    try
        load(['./LL/' prefix num2str(i) '_' num2str(model_num)]);
    catch
        fprintf('\nCould not load model %g for subject %g!\n',model_num,i)
    end
    
    LL_Cell{i} = LL;
end

end

function run_model_a_on_subj_b_with_c_samples_new(model_num,subj_num,samples,is_color,set_size_idx)

load parameter_ranges_new

% choose which type of data to fit
switch is_color
    case 1
        load ../data/subjColorCell
        Data = revert_data(subjColorCell{subj_num});
        prefix = 'Color_';
    case 0
        load ../data/subjCell
        Data = revert_data(subjCell{subj_num});
        prefix = '';
    case 2
        load subjFakeCell
        Data = subjFakeCell{subj_num};
        prefix = 'Fake_';
    case 3
        load predictions_ort
        Data = data_fake(((set_size_idx-1)*91+1):(set_size_idx*91),:);
        prefix = 'Bin_';
    case 4
        load predictions_color
        Data = data_fake_color(((set_size_idx-1)*91+1):(set_size_idx*91),:);
        prefix = 'cBin_';
end

% run model on specified subject
fprintf('Generating predictions for model: %s on set size index: %i\n',model_names(model_num,:),set_size_idx)
t1 = tic;
switch model_num
    
    case 1
        [LL P_C_HAT max_P_C_HAT] = run_model_EN(Data,meanvecNPER,priorvec,samples);
    case 2
        [LL P_C_HAT max_P_C_HAT] = run_model_ER(Data,meanvec,priorvec,samples);
    case 3
        [LL P_C_HAT max_P_C_HAT] = run_model_VR(Data,meanvecVR,thetavec,priorvec,samples);
    case 4
        k_temp = 1:8;
        [LL P_C_HAT max_P_C_HAT] = run_model_SA(Data,slotJvec,k_temp,priorvec,samples);
    case 5
        k_temp = 1:8;
        [LL P_C_HAT max_P_C_HAT] = run_model_IL(Data,lapsevec,k_temp,qvec);
    case 6
        [LL P_C_HAT max_P_C_HAT] = run_model_VN(Data,meanvecVR,thetavec,priorvec,samples);
    case 7
        [LL P_C_HAT max_P_C_HAT] = run_model_EG(Data,meanvecNPER,lapsevec,priorvec,samples);
    case 8
        % SA with random uniform K
        [LL P_C_HAT max_P_C_HAT] = run_model_SU(Data,slotJvec,kvec,priorvec,samples);
    case 9
        % Item limit with ER
        k_temp = 1:8;
        [LL P_C_HAT max_P_C_HAT] = run_model_EI(Data,meanvec,k_temp,priorvec,samples);
    case 10
        % Uniform random item limit with ER
        [LL P_C_HAT max_P_C_HAT] = run_model_EU(Data,meanvec,kvec,priorvec,samples);
    case 11
%         [LL P_C_HAT max_P_C_HAT] = run_model_VR_constK(Data,meanvec,thetavec,priorvec,samples);        % OUT OF ORDER
    case 12
        k_temp = 1:8;
        [LL P_C_HAT max_P_C_HAT] = run_model_VI(Data,meanvecVR,thetavec,k_temp,priorvec,samples);
    case 15
        [LL P_C_HAT max_P_C_HAT] = run_model_VP_new(Data,meanvecVR,thetavec,alphavec,priorvec,samples);
    case 16
        [LL P_C_HAT max_P_C_HAT] = run_model_EP(Data,meanvec,alphavec,priorvec,samples);
    case 17
        [LL P_C_HAT max_P_C_HAT] = run_model_SR(Data,meanvec,kvec,priorvec,samples);
end

t2 = toc(t1);
time_completed = datestr(clock);
runtime = [num2str(t2/60) ' minutes'];

% save model results
if is_color <3
    save(['LL/' prefix 'LL_' num2str(subj_num) '_' num2str(model_num)],'LL','samples','time_completed','runtime','-v7.3');
end

if exist('P_C_HAT','var')
    save([ 'P_C_HAT/' prefix 'P_C_HAT_' num2str(subj_num) '_' num2str(model_num) '_' num2str(set_size_idx)],'P_C_HAT','samples','time_completed','runtime','-v7.3');
else
    fprintf('\nP_C_HAT not produced or saved!\n')
end

if is_color <3

if exist('max_P_C_HAT','var')
    save(['max_P_C_HAT/' prefix 'max_P_C_HAT_' num2str(subj_num) '_' num2str(model_num)],'max_P_C_HAT','time_completed','runtime','samples','-v7.3');
else
    fprintf('\nmax_P_C_HAT not produced or saved!\n')
end
end

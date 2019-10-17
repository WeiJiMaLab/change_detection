% helper function for compute_psych_curves_multiple
function [p_C_hat_mat HR FA] = compute_psych_curves(Subj_Data_Cell,num_delta_bins,make_plot_curr,Model_Data_Cell)

% do different things depending on the number of arguments in. 
switch nargin
    case 2 % don't plot anything
        make_plot = 0;
        make_plot_model = 0;
    case 3 % turn on plotting 
        make_plot = make_plot_curr;
        make_plot_model = 0;
    case 4 % plot subject data and model predictions
        
        % recursively calls itself
        [model_p_C_hat_mat model_HR model_FA] = compute_psych_curves(Model_Data_Cell,num_delta_bins,0);
        if make_plot_curr
            make_plot_model = 1;
            make_plot = 0;
        else
            make_plot_model = 0;
            make_plot = 0;
        end
end

% Get P_C_hat = 1 as a function of delta
% first bin will contain only delta == 0 points
delta_bin_vec = linspace(eps,pi/2+eps,num_delta_bins);

% Loop over subjects
num_subj = length(Subj_Data_Cell);

% intialize matrices for holding statistics
N_vec = unique(Subj_Data_Cell{1}(:,5));
p_C_hat_mat = zeros(num_subj,length(N_vec),length(delta_bin_vec));
HR = zeros(num_subj,length(N_vec));
FA= HR;

for subj_idx = 1:num_subj
    
    % get current subject data
    curr_data = Subj_Data_Cell{subj_idx};
    curr_data_delta = .5*sum(abs(circ_dist((pi/90)*curr_data(:,56:63),(pi/90)*curr_data(:,64:71))),2);
    
    % loop over set size
    for N_idx = 1:length(N_vec)
        
        % current set size
        curr_N = curr_data(:,5) == N_vec(N_idx);
        curr_N_data = curr_data(curr_N,:);
        curr_N_delta = curr_data_delta(curr_N);
        
        % first data point contains all no-change stimuli
        resp_C = mean(curr_N_data((curr_N_delta<eps),2),1);
        p_C_hat_mat(subj_idx,N_idx,1)= resp_C;
        
        % loop over deltas
        for delta_idx = 2:num_delta_bins
            
            resp_C = mean(curr_N_data((curr_N_delta>delta_bin_vec(delta_idx-1)) & ...
                (curr_N_delta<=delta_bin_vec(delta_idx)),2),1);
            
            p_C_hat_mat(subj_idx,N_idx,delta_idx)= resp_C;
            
        end
        
        % compute hit rate and false alarm rate
        HR(subj_idx,N_idx) = sum(curr_N_data((curr_N_data(:,1)~=0),2),1)/sum((curr_N_data(:,1)~=0),1);
        FA(subj_idx,N_idx) = sum(curr_N_data((curr_N_data(:,1)==0),2),1)/sum((curr_N_data(:,1)==0),1);
        
    end
    
end

% compute means and standard error
curr_mean = squeeze(mean(p_C_hat_mat,1));
curr_stderr = squeeze(std(p_C_hat_mat,[],1)/sqrt(num_subj));
FA_mean = mean(FA,1);
FA_stderr = std(FA,[],1)/sqrt(num_subj);

HR_mean = mean(HR,1);
HR_stderr = std(HR,[],1)/sqrt(num_subj);

% plot subject results
if make_plot
    close
    figure
    subplot(1,2,1)
    stdlines = [.9 .6 0;1 0 0;0 1 0;0 0 1];
    delta_bin_vec = [0 diff(delta_bin_vec)/2+delta_bin_vec(1:(end-1))];
    for i = 1:length(N_vec)
        errorbar(delta_bin_vec,curr_mean(i,:),curr_stderr(i,:),'Color',stdlines(i,:));
        hold on
    end
    legend([repmat('N=',length(N_vec),1) num2str(N_vec)],'Location','SouthEast')
    title('Probability report "Change"')
    xlabel('Magnitude of change in degrees')
    xlim([-.05 pi/2+.05])
    ylim([0 1])
    set(gca,'YTick',[0 .2 .4 .6 .8 1]);
    set(gca,'XTick',[0 pi/6 pi/3 pi/2]);
    set(gca,'XTickLabel',{'0','30','60','90'});
    set(gca,'YTickLabel',{'0','0.2','0.4','0.6','0.8','1'});
    
    subplot(1,2,2)
    
    HR_point_color = [0 0 0];
    FA_point_color = (1/255)*[148 138 84];
    
    hold on
    errorbar(N_vec,HR_mean,HR_stderr,'Color',HR_point_color);
    errorbar(N_vec,FA_mean,FA_stderr,'Color',FA_point_color);
    legend('Hit Rates','False Alarms')
    title('Hit and false alarm rates')
    
    xlabel('Set Size')
end

% plot model predictions as polygons on the original data
if make_plot_model
    figure
    stdcols = [1 .8 .5;1 .8 .8; .8 1 .8; .8 .8 1];
    subplot(1,2,1)
    hold on
    
    polyX = [delta_bin_vec delta_bin_vec(end:-1:1)];
    
    model_p_C_hat_mat_mean = squeeze(mean(model_p_C_hat_mat,1));
    model_p_C_hat_mat_stderr = squeeze(std(model_p_C_hat_mat,[],1)/sqrt(num_subj));
    
    if num_subj > 1
    for i = 1:length(N_vec)
        polyY(1:length(delta_bin_vec)) = model_p_C_hat_mat_mean(i,:)-model_p_C_hat_mat_stderr(i,:);
        polyY((2*length(delta_bin_vec)):-1:(length(delta_bin_vec)+1)) = model_p_C_hat_mat_mean(i,:)+model_p_C_hat_mat_stderr(i,:);
%         fill(polyX,polyY,stdcols(i,:),'LineStyle','none','facealpha',0.8);
        fill(polyX,polyY,stdcols(i,:),'LineStyle','none');
    end
    else
        for i = 1:length(N_vec)
        errorbar(delta_bin_vec,model_p_C_hat_mat_mean(i,:),model_p_C_hat_mat_stderr(i,:),'Color',.5*stdcols(i,:))
        end
    end
    
    
    stdlines = [.9 .6 0;1 0 0;0 1 0;0 0 1];
    hold on
    for i = 1:length(N_vec)
        errorbar(delta_bin_vec,curr_mean(i,:),curr_stderr(i,:),'o','Color',stdlines(i,:),'LineWidth',2);
    end
    legend([repmat('N=',length(N_vec),1) num2str(N_vec)],'Location','SouthEast')
    xlim([-.05 (pi/2+.05)])
    ylim([0 1])
    set(gca,'YTick',[.2 .4 .6 .8 1]);
    set(gca,'XTick',[0 pi/6 pi/3 pi/2]);
    set(gca,'XTickLabel',[]);
    set(gca,'YTickLabel',[]);
    title('Probability report "Change"')
    xlabel('Magnitude of change in radians')
    
    subplot(1,2,2)
    % HR and FA
    hold on
    
    polyX = [N_vec' N_vec(end:-1:1)'];
    
    HR_fill_color = 191/255*[1 1 1];
    HR_point_color = [0 0 0];
    FA_fill_color = (1/255)*[221 217 195];
    FA_point_color = (1/255)*[148 138 84];


    model_HR_mean = (mean(model_HR,1));
    model_HR_stderr = (std(model_HR,[],1)/sqrt(num_subj));
    model_FA_mean = (mean(model_FA,1));
    model_FA_stderr = (std(model_FA,[],1)/sqrt(num_subj));
    
    if num_subj >1
    polyY_HR(1:length(N_vec)) = model_HR_mean-model_HR_stderr;
    polyY_HR((2*length(N_vec)):-1:(length(N_vec)+1)) = model_HR_mean+model_HR_stderr;
    fill(polyX,polyY_HR,HR_fill_color,'LineStyle','none');
    polyY_FA(1:length(N_vec)) = model_FA_mean-model_FA_stderr;
    polyY_FA((2*length(N_vec)):-1:(length(N_vec)+1)) = model_FA_mean+model_FA_stderr;
    fill(polyX,polyY_FA,FA_fill_color,'LineStyle','none');
    else
        errorbar(N_vec,model_HR_mean,model_HR_stderr,'Color',.5*HR_fill_color)
        errorbar(N_vec,model_FA_mean,model_FA_stderr,'Color',.5*FA_fill_color)
    end
    
    
    errorbar(N_vec,HR_mean,HR_stderr,'o','Color',HR_point_color,'LineWidth',2);
    errorbar(N_vec,FA_mean,FA_stderr,'o','Color',FA_point_color,'LineWidth',2);
    xlim([(min(N_vec)-1) (max(N_vec)+1)])
    set(gca,'YTick',[.2 .4 .6 .8 1]);
    set(gca,'XTick',N_vec);
    set(gca,'XTickLabel',[]);
    set(gca,'YTickLabel',[]);
    ylim([0 1])
    legend('Hit Rates','False Alarms')
    title('Hit and false alarm rates')
    
    xlabel('Set Size')
    
    set(gcf,'OuterPosition',[0 0 1000 500])
end
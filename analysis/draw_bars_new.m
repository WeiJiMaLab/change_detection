% Compute BMC results and draw the bar plot

function [bars_mean,bars_std,bars_color_mean,bars_color_std] = draw_bars_new()

x_vec = 1:4;

% compute orientation BMC and subtract the VP results
bars = compute_BMC([5 4 17 16 15],0);
bars_temp = bars - repmat(bars(:,5),1,5);
bars_temp = bars_temp(:,1:4);

% do the same for color
bars_color = compute_BMC([5 4 17 16 15],1);
bars_color_temp = bars_color - repmat(bars_color(:,5),1,5);
bars_color_temp = bars_color_temp(:,1:4);

% compute mean and standard error for both
bars_mean = mean(bars_temp,1);
bars_std = std(bars_temp,[],1)/sqrt(size(bars,1));

bars_color_mean = mean(bars_color_temp,1);
bars_color_std = std(bars_color_temp,[],1)/sqrt(size(bars_color,1));

% plot the results for orientation and color
subplot(1,2,1)
bar(x_vec,bars_mean)
hold on
errorbar(x_vec,bars_mean,bars_std,'.')
title('Orientation BMC (relative to PL VP)')
set(gca,'XTickLabel',{'IL','SA','1/N SR','PL EP'}')
% uncomment to plot bar values at the end of the bars
% text(x_vec(1:(end-1)) - .5,bars_mean(1:(end-1))+ sign(bars_mean(1:(end-1)))*10,num2str(bars_mean(1:(end-1))'))
hold off

subplot(1,2,2)
bar(x_vec,bars_color_mean)
hold on
errorbar(x_vec,bars_color_mean,bars_color_std,'.')
title('Color BMC (relative to PL VP)')
set(gca,'XTickLabel',{'IL','SA','1/N SR','PL EP'}')
% uncomment to plot bar values at the end of the bars
% text(x_vec(1:(end-1)) - .5,bars_color_mean(1:(end-1))+ sign(bars_color_mean(1:(end-1)))*10,num2str(bars_color_mean(1:(end-1))'))
set(gcf,'Name','BMC')

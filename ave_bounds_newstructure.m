function [result_ave,upper,lower] = ave_bounds_newstructure(result)
%ave_bounds takes in the desired quantity and produces the upper/lower
%bound lengths for use with error bar plotting

result_ave = mean(result(1,:),'omitnan');
result_std = std(result(1,:),'omitnan');
uncert = result(2,:);
weight = 1./uncert.^2;
if sum(result(2,:)) == 0
    result_stdmean = result_std/sqrt(length(result(:,2)));
else
    % result_stdmean = sqrt(sum(result(2,:).^2)/length(result(2,:))); %This is pooled statistics
    result_ave = sum(result(1,:).*weight)/sum(weight);
    result_stdmean = sqrt(1/sum(weight));
end
result_stdmean = result_stdmean*tinv(0.975,length(result)); %Multiplying to get 95%
upper = result_stdmean;
lower = result_stdmean;
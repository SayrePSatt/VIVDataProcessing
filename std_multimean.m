function pooledStat = std_multimean(numSamp,testAve,testStd)
%Calculates the average of averages and mean standard deviation
aveofAve = mean(testAve);
meanStd = sqrt(sum(testStd.^2)/(sum(numSamp)-numel(numSamp)));
pooledStat = [aveofAve; meanStd];
end
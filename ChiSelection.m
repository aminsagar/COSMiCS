%% Create initial estimations dataset 
firstcurve = horzcat(valoresS, Intensities(:,TwoEstimates(1)), Errors(:,TwoEstimates(1)));
for i = 1:length(Intensities(1,:))
curvetocompare = horzcat(valoresS,Intensities(:,i));
chiforcomp1(i,1) = i;
chiforcomp1(i,2) = compare2curves(firstcurve, curvetocompare);
end
[~,idx] = sort(chiforcomp1(:,2));
sortedchi1values = chiforcomp1(idx,:);
secondcurve = horzcat(valoresS, Intensities(:,TwoEstimates(2)), Errors(:,TwoEstimates(2)));
for i = 1:length(Intensities(1,:))
curvetocompare = horzcat(valoresS, Intensities(:,i));
chiforcomp2(i,1) = i;
chiforcomp2(i,2) = compare2curves(secondcurve, curvetocompare);
end
[~,idx] = sort(chiforcomp2(:,2));
sortedchi2values = chiforcomp2(idx,:);

for i = 1:length(sortedchi1values)
curvenumber = sortedchi1values(i,1);
[~,idx]=min(abs(curvenumber-sortedchi2values(:,1)));
chirankavglist(i,1) = curvenumber;
chirankavglist(i,2) = (i+idx)/2;
end
[~,idx] = sort(chirankavglist(:,1));
sortedchirankavg = chirankavglist(idx,:);
n = 1;
chiforcompavg = (chiforcomp1 + chiforcomp2)/2;
% [M,I] = max (chiforcompavg);
% chicompcomb(:,1) = chiforcomp1;
% chicompcomb(:,2) = chiforcomp2;
% chistddev(:,1) = std(chicompcomb,0,2);
[M,I] = max (sortedchirankavg(:,2));
ThirdCurveNumber = I;

PureComponents = horzcat(TwoEstimates(1),TwoEstimates(2),I);
figure()
subplot(2,1,1)
hold on
scatter (chirankavglist(:,1), chirankavglist(:,2), chirankavglist(:,2), chirankavglist(:,2), 'filled');
[M,I] = max(chirankavglist(:,2));
scatter(chirankavglist(I,1), M, M,'black','LineWidth',1.5);
xlabel ('Curve Number', 'FontSize', 18)
ylabel ('Average Rank', 'FontSize', 18)
title('Average Chi Rank', 'FontSize',18);
subplot(2,1,2)


for i = 1:length(PureComponents)
    PureComponents(i)
    semilogy(Intensities(:,PureComponents(i)), 'LineWidth',2);
    ylabel ('Intensity', 'FontSize', 18)
    xlabel ('q', 'FontSize', 18)
    hold on
    title('Selected Estimates', 'FontSize',18);
end
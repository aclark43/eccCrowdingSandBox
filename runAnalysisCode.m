%% Load Data
load('X:\Ashley\EccentriciyCrowdingData\eccCrowdingData.mat')


%% Plot All Thresholds in logMAR
figure;
for ii = 1:length(subjectsAll)
    for c = 1:2
        for eNum = 1:4
            idxS1 = find([dataT(:).subject] == ii & ...
                [dataT(:).condition] == c & ...
                [dataT(:).ecc] == eNum);
            plot(eccVals(eNum), dataT(idxS1).logMAR,'o',...
                'MarkerFaceColor',[.7 .7 .7],'MarkerEdgeColor',[.7 .7 .7]);
            hold on
        end
    end
end
idxS1 = find([dataT(:).condition] == 1 & ...
    [dataT(:).ecc] == 1);
idxS2 = find([dataT(:).condition] == 1 & ...
    [dataT(:).ecc] == 2);
idxS3 = find([dataT(:).condition] == 1 & ...
    [dataT(:).ecc] == 3);
idxS4 = find([dataT(:).condition] == 1 & ...
    [dataT(:).ecc] == 4);

hold on;
Forlegs(1) = errorbar([0 10 15 25],...
    [mean([dataT(idxS1).logMAR]), mean([dataT(idxS2).logMAR]),...
    mean([dataT(idxS3).logMAR]),mean([dataT(idxS4).logMAR])],...
    [sem([dataT(idxS1).logMAR]), sem([dataT(idxS2).logMAR]),...
    sem([dataT(idxS3).logMAR]),sem([dataT(idxS4).logMAR])],'o',...
    'MarkerFaceColor','b','MarkerEdgeColor','b');

md1 = fitlm([0 10 15 25],...
    [mean([dataT(idxS1).logMAR]), mean([dataT(idxS2).logMAR]),...
    mean([dataT(idxS3).logMAR]),mean([dataT(idxS4).logMAR])]);
x = 0:25; % Defines the domain as [-15,25] with a refinement of 0.25
m  = md1.Coefficients.Estimate(2);  % Specify your slope
x1 = 0; % Specify your starting x
y1 = mean([dataT(idxS1).logMAR]);  % Specify your starting y
y = m*(x - x1) + y1;
hPlot = plot(x,y,'--b'); 
% text(5,-0.05,sprintf('b = %.2f',m))
hold on

[~,p1] = ttest([dataT(idxS1).logMAR],[dataT(idxS2).logMAR])
[~,p2] = ttest([dataT(idxS1).logMAR],[dataT(idxS3).logMAR])
[~,p3] = ttest([dataT(idxS1).logMAR],[dataT(idxS4).logMAR])

groups = {[0 10],[0 15],[0 25]};
H=sigstar(groups,[p1 p2 p3]);   


idxS1 = find([dataT(:).condition] == 2 & ...
    [dataT(:).ecc] == 1);
idxS2 = find([dataT(:).condition] == 2 & ...
    [dataT(:).ecc] == 2);
idxS3 = find([dataT(:).condition] == 2 & ...
    [dataT(:).ecc] == 3);
idxS4 = find([dataT(:).condition] == 2 & ...
    [dataT(:).ecc] == 4);

Forlegs(2) = errorbar([0 10 15 25],...
    [mean([dataT(idxS1).logMAR]), mean([dataT(idxS2).logMAR]),...
    mean([dataT(idxS3).logMAR]),mean([dataT(idxS4).logMAR])],...
    [sem([dataT(idxS1).logMAR]), sem([dataT(idxS2).logMAR]),...
    sem([dataT(idxS3).logMAR]),sem([dataT(idxS4).logMAR])],'o',...
    'MarkerFaceColor','r','MarkerEdgeColor','r');

md1 = fitlm([0 10 15 25],...
    [mean([dataT(idxS1).logMAR]), mean([dataT(idxS2).logMAR]),...
    mean([dataT(idxS3).logMAR]),mean([dataT(idxS4).logMAR])]);
x = 0:25; % Defines the domain as [-15,25] with a refinement of 0.25
m  = md1.Coefficients.Estimate(2);  % Specify your slope
x1 = 0; % Specify your starting x
y1 = mean([dataT(idxS1).logMAR]);  % Specify your starting y
y = m*(x - x1) + y1;
hPlot = plot(x,y,'--r');
% text(5,0.1,sprintf('b = %.2f',m))

ylabel('Acuity (logMAR)');
xlabel('Eccentricity (arcmin)');

xlim([-5 30])
xticks([0 10 15 25])
legend(Forlegs,{'Uncrowded','Crowded'},'Location','northwest')
[~,p1] = ttest([dataT(idxS1).logMAR],[dataT(idxS2).logMAR])
[~,p2] = ttest([dataT(idxS1).logMAR],[dataT(idxS3).logMAR])
[~,p3] = ttest([dataT(idxS1).logMAR],[dataT(idxS4).logMAR])

groups = {[0 10],[0 15],[0 25]};
H=sigstar(groups,[p1 p2 p3]);   
% legend off

saveas(gcf, '../../../CrowdingTopology/Documents/Overleaf_CrowdingTopology/figures/acrossTheFovea/CvsU.png');
saveas(gcf, '../../../CrowdingTopology/Documents/Overleaf_CrowdingTopology/figures/acrossTheFovea/CvsU.epsc');


%% Create HeatMaps 2, 5
figure('Position',[2000 0 800 8000]);
counter = 1;
for ii = 1:length(subjectsAll)
    for e = 1:4
        xAll = [];
        yAll = [];
        currentPath = subjectCond{ii...
            }.Uncrowded.Unstabilized.(eccNames{e}).threshInfo.em.allTraces;
        for i = 1:length(currentPath.x)
            xAll = [currentPath.x{i} xAll];
            yAll = [currentPath.y{i} yAll];
        end   
             
        %         figure;
        subplot(6,4,counter)
        if ii == 1 && e == 1
            heatMapIm{ii,e} = generateHeatMapSimple( ...
                xAll, ...
                yAll, ...
                'Bins', 40,...
                'StimulusSize', 0,...
                'AxisValue', 15,...
                'Uncrowded', 0,...
                'Borders', 1);
        else heatMapIm{ii,e} = generateHeatMapSimple( ...
                xAll, ...
                yAll, ...
                'Bins', 40,...
                'StimulusSize', 0,...
                'AxisValue', 15,...
                'Uncrowded', 0,...
                'Borders', 0);
        end
        percentWithin5(ii,e) = length(find(xAll < 5 & xAll > -5 & ...
            yAll < 5 & yAll > -5))/length(xAll)*100;
        
        
        
        axis square
        line([0 0],[-15 15],'Color','k')
        line([-15 15],[0,0],'Color','k')
        if ii == 1
            title(sprintf(' %i Eccentricity',eccVals(e)))
        end
        if e == 1
            ylabel(ii)
        end
        counter = counter + 1;
    end
end
% saveas(gcf, '../../../Grants/MPRO1_2022/figures/Ashley10arcHM.png');
% saveas(gcf, '../../../Grants/MPRO1_2022/figures/Ashley10arcHM.epsc');
suptitle('Uncrowded');
saveas(gcf, '../../../CrowdingTopology/Documents/Overleaf_CrowdingTopology/figures/acrossTheFovea/EMmaps.png');
saveas(gcf, '../../../CrowdingTopology/Documents/Overleaf_CrowdingTopology/figures/acrossTheFovea/EMmaps.epsc');

figure('Position',[2000 0 800 8000]);
counter = 1;
for ii = 1:length(subjectsAll)
    for e = 1:4
        xAll = [];
        yAll = [];
        if ii == 3 && e == 1 ||...
                ii == 4 && e == 1
            namesAll = subjectCond{ii}.Crowded.Unstabilized.Center0ecc.threshInfo.em.ecc_0;
            namesPass = fieldnames(namesAll);
            for i = 1:length(namesPass)
                for j = 1:length(namesAll.(namesPass{i}).position)
                    xAll = [namesAll.(namesPass{i}).position(j).x xAll];
                    yAll = [namesAll.(namesPass{i}).position(j).y yAll];
                end
            end
        else currentPath = subjectCond{ii...
                }.Crowded.Unstabilized.(eccNames{e}).threshInfo.em.allTraces;
            
            for i = 1:length(currentPath.x)
                xAll = [currentPath.x{i} xAll];
                yAll = [currentPath.y{i} yAll];
            end
        end
        %         figure;
        subplot(6,4,counter)
        if ii == 1 && e == 1
            heatMapIm2{ii,e} = generateHeatMapSimple( ...
                xAll, ...
                yAll, ...
                'Bins', 40,...
                'StimulusSize', 0,...
                'AxisValue', 15,...
                'Uncrowded', 0,...
                'Borders', 1);
        else heatMapIm2{ii,e} = generateHeatMapSimple( ...
                xAll, ...
                yAll, ...
                'Bins', 40,...
                'StimulusSize', 0,...
                'AxisValue', 15,...
                'Uncrowded', 0,...
                'Borders', 0);
        end
        axis square
        line([0 0],[-15 15],'Color','k')
        line([-15 15],[0,0],'Color','k')
        if ii == 1
            title(sprintf(' %i Eccentricity',eccVals(e)))
        end
        if e == 1
            ylabel(ii)
        end
        counter = counter + 1;
    end
end
% saveas(gcf, '../../../Grants/MPRO1_2022/figures/Ashley10arcHM.png');
% saveas(gcf, '../../../Grants/MPRO1_2022/figures/Ashley10arcHM.epsc');
suptitle('Crowded');
saveas(gcf, '../../../CrowdingTopology/Documents/Overleaf_CrowdingTopology/figures/acrossTheFovea/EMmapsCrow.png');
saveas(gcf, '../../../CrowdingTopology/Documents/Overleaf_CrowdingTopology/figures/acrossTheFovea/EMmapsCrow.epsc');


%%% Create Combined Map
figure;
counterU = 1;
counterC = 5;
for e = 1:4
    newMatU = (...
        heatMapIm{1,e}+...
        heatMapIm{2,e}+...
        heatMapIm{3,e}+...
        heatMapIm{4,e}+...
        heatMapIm{5,e}+...
        heatMapIm{6,e})/6;
    %     figure;
    subplot(2,4,counterU)
    hold on
    temp = pcolor(linspace(-15, 15, size(newMatU, 1)),...
        linspace(-15, 15, size(newMatU, 1)),...
        newMatU');
    hold on
    load('./MyColormaps.mat');
    mycmap(1,:) = [1 1 1];
    set(gcf, 'Colormap', mycmap)
    axis([-15 15 -15 15])
    shading interp;
    line([-15 15],[0 0],'Color','k');
    line([0,0],[-15 15],'Color','k');
    axis square
    if e == 1
        ylabel({'Uncrowded','Y (arcmin)'});
        xlabel('X (arcmin)');
    else
        axis off
    end
    title(sprintf('Eccentricity %i',eccVals(e)))
    
    newMatC = (...
        heatMapIm2{1,e}+...
        heatMapIm2{2,e}+...
        heatMapIm2{3,e}+...
        heatMapIm2{4,e}+...
        heatMapIm2{5,e}+...
        heatMapIm2{6,e})/6;
    subplot(2,4,counterC)
    hold on
    temp = pcolor(linspace(-15, 15, size(newMatC, 1)),...
        linspace(-15, 15, size(newMatC, 1)),...
        newMatC');
    hold on
    load('./MyColormaps.mat');
    mycmap(1,:) = [1 1 1];
    set(gcf, 'Colormap', mycmap)
    axis([-15 15 -15 15])
    shading interp;
    line([-15 15],[0 0],'Color','k');
    line([0,0],[-15 15],'Color','k');
    axis square
    counterU = counterU + 1;
    counterC = counterC + 1;
     if e == 1
        ylabel({'Crowded','Y (arcmin)'});
        xlabel('X (arcmin)');

     else
         axis off
     end
end
saveas(gcf, '../../../CrowdingTopology/Documents/Overleaf_CrowdingTopology/figures/acrossTheFovea/EMmapsComb.png');
saveas(gcf, '../../../CrowdingTopology/Documents/Overleaf_CrowdingTopology/figures/acrossTheFovea/EMmapsComb.epsc');

%% Percent Differences
idxU1 = find([dataT(:).condition] == 1 & ...
    [dataT(:).ecc] == 1);
idxC3 = find([dataT(:).condition] == 2 & ...
    [dataT(:).ecc] == 3);
idxC4 = find([dataT(:).condition] == 2 & ...
    [dataT(:).ecc] == 4);

% Define initial and final performance values
initial_performance = mean([dataT(idxU1).thresh]); % Replace with the actual initial performance value
final_performance = mean([dataT(idxC4).thresh]);    % Replace with the actual final performance value

% Calculate percent drop in performance
percent_drop = ((initial_performance - final_performance) / initial_performance) * 100;

%% ANOVA

idxC1 = find([dataT(:).condition] == 2 & ...
    [dataT(:).ecc] == 1);
idxC2 = find([dataT(:).condition] == 2 & ...
    [dataT(:).ecc] == 2);
idxC3 = find([dataT(:).condition] == 2 & ...
    [dataT(:).ecc] == 3);
idxC4 = find([dataT(:).condition] == 2 & ...
    [dataT(:).ecc] == 4);
idxU1 = find([dataT(:).condition] == 1 & ...
    [dataT(:).ecc] == 1);
% Example data (replace with your actual data)
% Example data for 5 groups (replace with your actual data)
group1_data = [dataT(idxU1).thresh];
group2_data = [dataT(idxC1).thresh];
group3_data = [dataT(idxC2).thresh];
group4_data = [dataT(idxC3).thresh];
group5_data = [dataT(idxC4).thresh];

% Combine data into a single column vector
data = [group1_data; group2_data; group3_data; group4_data; group5_data];

% Create a grouping variable
group = [ones(6, 1); 2 * ones(6, 1); 3 * ones(6, 1); 4 * ones(6, 1); 5 * ones(6, 1)];

% Perform one-way ANOVA
[p, tbl, stats] = anova1(data(:), group, 'on');

% Perform Tukey-Kramer post hoc test
c = multcompare(stats);

% Display results
disp('Tukey-Kramer Post Hoc Test Results:');
disp(c);

% Identify significant group differences
significant_groups = find(c(:, end) < 0.05);

% Display significant group differences
disp('Significant Group Differences:');
for i = 1:size(significant_groups, 1)
    group1 = c(significant_groups(i), 1);
    group2 = c(significant_groups(i), 2);
    fprintf('Groups %d and %d are significantly different.\n', group1, group2);
end


%% Paired ttest
idxC1 = find([dataT(:).condition] == 2 & ...
    [dataT(:).ecc] == 1);
idxC2 = find([dataT(:).condition] == 2 & ...
    [dataT(:).ecc] == 2);
idxC3 = find([dataT(:).condition] == 2 & ...
    [dataT(:).ecc] == 3);
idxC4 = find([dataT(:).condition] == 2 & ...
    [dataT(:).ecc] == 4);
idxU1 = find([dataT(:).condition] == 1 & ...
    [dataT(:).ecc] == 1);
idxU2 = find([dataT(:).condition] == 1 & ...
    [dataT(:).ecc] == 2);
idxU3 = find([dataT(:).condition] == 1 & ...
    [dataT(:).ecc] == 3);
idxU4 = find([dataT(:).condition] == 1 & ...
    [dataT(:).ecc] == 4);

% Example data for 10 conditions with paired observations (replace with your actual data)
condition1 = [dataT(idxC1).thresh];
condition2 = [dataT(idxC2).thresh];
condition3 = [dataT(idxC3).thresh];
condition4 = [dataT(idxC4).thresh];
condition5 = [dataT(idxU1).thresh];
condition6 = [dataT(idxU2).thresh];
condition7 = [dataT(idxU3).thresh];
condition8 = [dataT(idxU4).thresh];


% Combine data into a cell array
conditions_data = {condition1, condition2, condition3, condition4, ...
                   condition5, condition6, condition7, condition8};

% Number of conditions
num_conditions = numel(conditions_data);

% Perform paired t-tests for all pairwise comparisons
p_values = zeros(num_conditions);

for i = 1:num_conditions
    for j = i+1:num_conditions
        [~, p_values(i, j)] = ttest(conditions_data{i}, conditions_data{j});
    end
end

% Display results
fprintf('Paired t-Test P-Values:\n');
for i = 1:num_conditions
    for j = i+1:num_conditions
        fprintf('Condition %d vs. Condition %d: %.4f\n', i, j, p_values(i, j));
    end
end



%% Crowded vs Uncrowded
figure('units','normalized','outerposition',[0.2 0.05 .4 .6]) %plots all absolute thresholds for each condition
counter = 1;
counterLeg = 1;
clear data
for ii = 1:length(subjectsAll)
    clear oneRow
    hold on
    for cIdx = 1:length(conditions)
        for sIdx = 1:length(stabilization)
            clear allValsSubj
            cond = conditions{cIdx};
            stab = stabilization{sIdx};
            for numEcc = 1:length(ecc)
                    oneRow(numEcc) = subjectCond{ii}.(cond).(stab).(eccNames{numEcc}).threshInfo.thresh;
                    boots(numEcc) = confInter95(subjectCond{ii}.(cond).(stab).(eccNames{numEcc}).threshInfo.threshB);
                    oneRowEcc(numEcc) = eccVals(numEcc);
                    numTrials(ii,numEcc) = sum(subjectCond{ii}.(cond).(stab).(eccNames{numEcc}).threshInfo.em.valid);
            end
            data(counter).condition = cIdx;
            data(counter).subject = ii;
            data(counter).thresholds = oneRow;
            data(counter).eccentricities = oneRowEcc;
            data(counter).numTrials = numTrials;
            data(counter).threshFreq = 60./((oneRow));
            data(counter).boots = boots;
            hold on
            %             yyaxis left
            if bitget(counter,1) %odd
%                 leg(counterLeg) = plot(eccVals(~isnan(oneRow))+1,oneRow(~isnan(oneRow)),'-o','Color',c1(ii,:),...
%                     'MarkerFace',c1(ii,:),'MarkerSize', 10, 'LineWidth',3,'MarkerEdgeColor','k');
                leg(counterLeg) = errorbar(eccVals(~isnan(oneRow))+1,oneRow(~isnan(oneRow)),...
                    boots(~isnan(oneRow)),boots(~isnan(oneRow)),'-o','Color',c1(ii,:),...
                    'MarkerFace',c1(ii,:),'MarkerSize', 10, 'LineWidth',3,'MarkerEdgeColor','k');
                counterLeg = counterLeg + 1;
            else %even
                plot(eccVals(~isnan(oneRow))+1,oneRow(~isnan(oneRow)),'-s','Color',c2(ii,:),...
                    'MarkerFace',c2(ii,:),'MarkerSize', 10, 'LineWidth',3,'MarkerEdgeColor','k');
            end
            
            hold on
           
            counter = counter + 1;
        end
        forBouma{cIdx}(ii,:) = oneRow;
    end
    
end
hold on
ylabel('Threshold (arcmin)');

set(gca,'xtick',eccVals,'xticklabel',string(eccVals), 'FontSize',14)
xlabel('Eccentricity (arcmin)')
xtickangle(60)
hold on
errorbar([0 10 15 25],mean(forBouma{1}),std(forBouma{1}),'-o',...
    'MarkerSize',10,'Color','k','LineWidth',4);
errorbar([0 10 15 25],mean(forBouma{2}),std(forBouma{2}),'-s',...
    'MarkerSize',10,'Color','k','LineWidth',4);

temp = ([mean(forBouma{1})])/2*20
temp2 = ([mean(forBouma{2})])/2*20
temp3 = [data(:).thresholds]
ecc0Th = temp3(1:4:length(temp3))'
ecc10th= temp3(2:4:length(temp3))'
ecc15th= temp3(3:4:length(temp3))'
ecc25th= temp3(4:4:length(temp3))'


[~,pU,bC,rU] = LinRegression([0 10 15 25],...
    mean(forBouma{1}),...
    0,NaN,1,0);
[~,pU,bC,rU] = LinRegression([0 10 15 25],...
    mean(forBouma{1}),...
    0,NaN,1,0);



[h,p0] = ttest(forBouma{1}(:,1)', forBouma{2}(:,1)');
[h,p10] = ttest(forBouma{1}(:,2)', forBouma{2}(:,2)');
[h,p15] = ttest(forBouma{1}(:,3)', forBouma{2}(:,3)');
[h,p25] = ttest(forBouma{1}(:,4)', forBouma{2}(:,4)');

ylim([1 3.5])
legend(leg, subjectsAll,...
    'Location','northwest');
title('All Thresholds per Subject');
saveas(gcf, '../../../CrowdingTopology/Documents/Overleaf_CrowdingTopology/figures/acrossTheFovea/AllThreshCrUncr.png');
saveas(gcf, '../../../CrowdingTopology/Documents/Overleaf_CrowdingTopology/figures/acrossTheFovea/AllThreshCrUncr.epsc');

%% EM and Slope for Acuity
% % % for ii = 1:length(subjectsAll)
% % % %     if ii == 3
% % % %         continue;
% % % %     end
% % %     mdl1 = fitlm([eccVals],forBouma(ii,:));
% % %     coefsTemp(ii) = mdl1.Coefficients{'x1','SE'};
% % %     emInfo(ii) = subjectCond{ii}.Uncrowded.Unstabilized.Center0ecc.threshInfo.em.allTraces.dCoefDsq;
% % % end
% % % figure;
% % % [~,p,bC,r] = LinRegression(coefsTemp,...
% % %     emInfo,...
% % %     0,NaN,1,0);
%% Percent change in VA
temp = mean(forBouma{1});
x1 = temp(1);
temp2 = mean(forBouma{2});
x2 = temp2(1);
% percentChange = (abs(x1 - x2)/((x1+x2)/2))*100;
percentChange = (abs(x1 - x2)/x1)*100;

x3 = temp2(2);
% percentChange2 = (abs(x3 - x2)/((x3+x2)/2))*100; INSTEAD x1-x3/x1, then
% DIFFERENCE percentChange2 = percentChange1
percentChange2 = (abs(x1 - x3)/x1)*100;
deffPercChange = percentChange2-percentChange


[~,p,bC,r] = LinRegression(eccVals, temp2,...
    0,NaN,1,0);

%% Bouma Law with Ecc Data too!
%%% plot again in terms of critical spacing (* 1.4)
% subjectsAll = {'Z002DDPI','AshleyDDPI','Z046DDPI','Z084DDPI','Z091','Z138'};
load ('MATFiles/peripheralCrowdingData'); %'peripheralData'
subjMatch = [0 5 2 6];
% 'Zoe', 'Z091', 'AshleyDDPI','Z138'
% create data set for MP
counter = length(data)+1;
for ii = 2:length(subjMatch)
    for i = 1:2
        data(counter).subject = subjMatch(ii);
        data(counter).condition = i;
        data(counter).thresholds = peripheralData{i}(ii,5:end);
        data(counter).eccentricities = peripheralData{3}(5:end);
        data(counter).numTrials = peripheralData{5}(ii,5:end);
        data(counter).threshFreq = 60./((peripheralData{i}(ii,5:end)));
        data(counter).boots = peripheralData{6}{i}(ii,5:end);
        counter = counter + 1;
    end
end
% save('FovealPeripheralAcuityCrowding.mat','data');
figure;
subplot(1,2,1)
    for ii = 1:length(subjectsAll)
    plot(eccVals, forBouma{2}(ii,:)* (1.4),'-o','LineWidth',3,...
        'Color',c1(ii,:),'MarkerFace',c1(ii,:)) %plotting in terms of critical spacing
    hold on
    end
    coefficients1 = polyfit(eccVals,mean(forBouma{2})*(1.4), 1);
%     xFit = linspace(min(eccVals), max(eccVals), 1000);
xFit = linspace(min(eccVals), 240, 1000);
    % xFitTemp - 
    yFit = polyval(coefficients1 , xFit);
    % m = (yFit(end)-yFit(1))/(xFit(end)-xFit(1)) ; %same as coefficents(1)
    temp = plot(xFit, yFit,'--','Color','k',...
        'MarkerSize',12,'MarkerFaceColor','k','LineWidth',3);
    hold on
    xticks([0 10 15 25])
    % plot(eccVals,mean(forBouma{2})*(1.4),'o','MarkerSize',12,'MarkerFaceColor','k');
    % [~,p,bCB,r] = LinRegression(eccVals, mean(forBouma{2})*(1.4),...
    %     0,NaN,1,0);
    axis square
    xlim([-5 30])
    ylim([2 5])
    xlabel('Eccentricity (arcmin)')
    ylabel('Critical Spacing (arcmin)')
    sigma = eccVals;
    w = mean(forBouma{2});
    BoumaLaw = (0.5 * sigma);
    hold on
    % temp(2) = plot(sigma, BoumaLaw,'--','Color','m',...
    %     'MarkerSize',12,'MarkerFaceColor','m','LineWidth',3);
    % hold on
    % temp(3) = plot(sigma, BoumaLaw+w,'--','Color','b',...
    %     'MarkerSize',12,'MarkerFaceColor','b','LineWidth',3);
    % legend(temp,{sprintf('Best Fit Line, b = %.2f',coefficients(1)),...
    %     'Bouma Law','Bouma + w'},'Location','northwest')



colorForPeriph = [1 0 0;c1(5,:);c1(2,:);c1(6,:)];
subplot(1,2,2)
% t = tiledlayout(2,1);
% ax1 = axes(t);
% figure;
    for ii = 1:size(peripheralData{1})
    plot(peripheralData{3}(5:8), peripheralData{2}(ii,5:8)*(1.4),'-o','LineWidth',3,...
        'Color',colorForPeriph(ii,:),'MarkerFace',colorForPeriph(ii,:)) %plotting in terms of critical spacing
    hold on
    end
    coefficients = polyfit(peripheralData{3}(5:8),...
        mean(peripheralData{2}(:,5:8))*(1.4), 1);
%     xFit2 = linspace(min(peripheralData{3}(5:8)), max(peripheralData{3}(5:8)), 1000);
xFit2 = linspace(30, 240, 1000);

yFit2 = polyval(coefficients , xFit2);
%     m = (yFit(end)-yFit(1))/(xFit(end)-xFit(1)) ; %same as coefficents(1)
    temp2 = plot(xFit2, yFit2,'--','Color','k',...
        'MarkerSize',12,'MarkerFaceColor','k','LineWidth',3);
%     hold on
%     plot(peripheralData{3}(5:8),mean(peripheralData{2}(:,5:8))*(1.4),...
%         'o','MarkerSize',12,'MarkerFaceColor','k');
    axis square
    xlim([40 380])
    ylim([4 80])
    xlabel('Eccentricity (arcmin)')
    ylabel('Critical Spacing (arcmin)')
    xticks([60 120 240 360])
    sigma = peripheralData{3}(5:8);
    w = mean(peripheralData{2}(:,5:8));
    BoumaLaw = (0.5 * sigma);
    hold on
%     temp2(2) = plot(sigma, BoumaLaw,'--','Color','m',...
%         'MarkerSize',12,'MarkerFaceColor','m','LineWidth',3);

%     legend(temp2,{sprintf('Best Fit Line, b = %.2f',coefficients(1)),...
%         'Boumas Law'},'Location','northwest')

    percentChange = (abs(0.15 - 0.05)/0.15)*100

% ax2 = axes(t);
% plot(ax2,peripheralData{3}(5:8)/60,mean(peripheralData{2}(:,5:8))*(1.4),...
%     'o','MarkerSize',12,'MarkerFaceColor','k')
% ax2.XAxisLocation = 'top';
% ax2.YAxisLocation = 'right';
% ax2.Color = 'none';
% ax1.Box = 'off';
% ax2.Box = 'off';

% Foveal and Periph Together
t = figure;
% ax1 = axes(t);
% plot(xFit, yFit,'--','Color','k',...
%     'MarkerSize',12,'MarkerFaceColor','k','LineWidth',3);
% xticks([0 30 60 120 240])
% xticklabels([0 0.5 1 2 4])
% xlim([0 250])
% ylim([0 40]);
% set(gca,'FontSize',16) 
% axis square
% xlabel('Eccentricity (deg)')
% ax1.XAxisLocation = 'top';
% 
% ax2 = axes(t);
temp3(1) = plot(xFit, yFit,'--','Color','k',...
    'MarkerSize',12,'MarkerFaceColor','k','LineWidth',3);
ax2.XAxisLocation = 'bottom';
hold on
temp3(2) = plot(xFit2, yFit2,'--','Color','b',...
    'MarkerSize',12,'MarkerFaceColor','b','LineWidth',3);
xticks([0 30 60 120 240])
ylim([0 40]);

xlabel('Eccentricity (arcmin)')
ylabel('Critical Spacing (arcmin)');
set(gca,'FontSize',16) 
xlim([0 250])
% legend(temp3,{'b = 0.05','b = 0.15'},'Location','northwest');
axis square
% temp3(1) = plot(xFit, yFit,'--','Color','k',...
%     'MarkerSize',12,'MarkerFaceColor','k','LineWidth',3);
% ax2.XAxisLocation = 'top';
% plot(ax2,x2,y2,'-k')


% saveas(gcf, '../../../CrowdingTopology/Documents/Overleaf_CrowdingTopology/figures/acrossTheFovea/PvsF.png');
% saveas(gcf, '../../../CrowdingTopology/Documents/Overleaf_CrowdingTopology/figures/acrossTheFovea/PvsF.epsc');

%%% plot Boumas law
%%% "critical spacing for identification of small letters is roughly half
%%% the eccentricity"; stim of 1 = cs of 1.4
sigma = eccVals;
w = mean(forBouma{1});
BoumaLaw = (0.5 * sigma);
hold on
leg3(1) = plot(sigma, BoumaLaw,'--','Color','r',...
    'MarkerSize',12,'MarkerFaceColor','k','LineWidth',3);
hold on
leg3(2) = plot(sigma, BoumaLaw+w,'--','Color','b',...
    'MarkerSize',12,'MarkerFaceColor','k','LineWidth',3);
xlim([0 30])
ylim([-1 10])
ylabel('Critical Spacing (arcmin)')
xlabel('Ecc Arcmin')
title('Critical Spacing x1')
% title('Critical Spacing x3')

legend(leg3,{'(0.5 * sigma)','(0.5 * sigma)+Width'});

% saveas(gcf, '../../../CrowdingTopology/Documents/Overleaf_CrowdingTopology/figures/acrossTheFovea/FoveavsPeriph.png');
% saveas(gcf, '../../../CrowdingTopology/Documents/Overleaf_CrowdingTopology/figures/acrossTheFovea/FoveavsPeriph.epsc');
% 
% 
% saveas(gcf, '../../../CrowdingTopology/Documents/Overleaf_CrowdingTopology/figures/acrossTheFovea/BoumaLawFovea.png');
% saveas(gcf, '../../../CrowdingTopology/Documents/Overleaf_CrowdingTopology/figures/acrossTheFovea/BoumaLawFoveax3.png');

%% Comparing with stabilized in AO machine

for ii = 1:length(subjectsAll)
    for c = 1:2
    idx1 = find([dataT.subject] == ii & ...
        [dataT.ecc] == 1 & ...
        [dataT.condition] == c);
    stimWidth = (dataT(idx1).acuity/20)*2;
    edgeSpacing0ecc(c,ii) = ((stimWidth)*1.4)-stimWidth;
    
      idx2 = find([dataT.subject] == ii & ...
        [dataT.ecc] == 3 & ...
        [dataT.condition] == c);
    stimWidth = (dataT(idx2).acuity/20)*2;
    edgeSpacing15ecc(c,ii) = ((stimWidth)*1.4)-stimWidth;
    end
end

temp = load('C:\Users\Ruccilab\Downloads\combinedCrowdingDF.mat');

figure;
for ii = 1:length(subjectsAll)
%     plot([0 15],[edgeSpacing0ecc(1,ii) edgeSpacing15ecc(1,ii)],'--o',...
%         'Color','b');
    hold on
   forlegs(1) =  plot([0 15],[edgeSpacing0ecc(2,ii) edgeSpacing15ecc(2,ii)],'--d',...
        'Color','r');

end
hold on

allSubjectsAO = fieldnames(temp.crowdingDF);
for ii = 1:length(fieldnames(temp.crowdingDF))
   forlegs(2) =  plot([0 15],[temp.crowdingDF.(allSubjectsAO{ii}).AtPRL.meanThreshE2E...
        temp.crowdingDF.(allSubjectsAO{ii}).Ecc15.meanThreshE2E],'--o',...
        'Color','g');
    hold on
end

ylabel('Edge-Edge Threshold')
xlabel('Eccentricity (0 = PRL)')
xlim([-5 20])
xticks([0 15])

legend(forlegs,{'Crowded DDPI','Crowded Stabilized AO'})
%% Plot just periph in degrees

changeEccToDeg = 60;
changeEccToDegANDEdgeToEdge = 60;

figure;
    for ii = 1:size(peripheralData{1})
    plot(peripheralData{3}(5:8)/changeEccToDeg, ...
        peripheralData{2}(ii,5:8)/changeEccToDegANDEdgeToEdge,'-o','LineWidth',3,...
        'Color',colorForPeriph(ii,:),'MarkerFace',colorForPeriph(ii,:)) %plotting in terms of critical spacing
    hold on
    end
    coefficients = polyfit(peripheralData{3}(5:8)/changeEccToDeg,...
        mean(peripheralData{2}(:,5:8))/changeEccToDegANDEdgeToEdge, 1);
%     xFit2 = linspace(min(peripheralData{3}(5:8)), max(peripheralData{3}(5:8)), 1000);
xFit2 = [];yFit2 = [];
xFit2 = linspace(1, 10, 1000);
yFit2 = polyval(coefficients , xFit2);
%     m = (yFit(end)-yFit(1))/(xFit(end)-xFit(1)) ; %same as coefficents(1)
    temp2 = plot(xFit2, yFit2,'--','Color','k',...
        'MarkerSize',12,'MarkerFaceColor','k','LineWidth',3);
ylabel('Horizontal C-C Spacing');
xlabel('Eccentricity (deg)');

hold on
coefficients = polyfit([1 8],[.5 4], 1);
yFit3 = polyval(coefficients , xFit2);
temp2 = plot(xFit2, yFit3,'--','Color','r',...
        'MarkerSize',12,'MarkerFaceColor','k','LineWidth',3);
    

%% Plot DC for each ecc
for ii = 1:length(subjectsAll)
    for e = 1:4
        dsqU(ii,e) = subjectCond{ii}.Uncrowded.Unstabilized.(eccNames{e}).threshInfo.em.allTraces.dCoefDsq;
    end
end
for ii = 1:length(subjectsAll)
    for e = 1:4
        dsqC(ii,e) = subjectCond{ii}.Crowded.Unstabilized.(eccNames{e}).threshInfo.em.allTraces.dCoefDsq; 
    end
end
figure
for ii = 1:length(subjectsAll)
    plot(eccVals,dsqU(ii,:),'--','b');
    hold on
end
for ii = 1:length(subjectsAll)
    plot(eccVals,dsqC(ii,:),'--','r');
    hold on
end

%% Normalized All to 0,0 Uncrowded
figure('units','normalized','outerposition',[0.2 0.05 .4 .6]) %plots all absolute thresholds for each condition
counter = 1;
counterLeg = 1;
for ii = 1:length(subjectsAll)
    clear oneRow
    hold on
    for cIdx = 1:length(conditions)
        for sIdx = 1:length(stabilization)
            clear allValsSubj
            cond = conditions{cIdx};
            stab = stabilization{sIdx};
            for numEcc = 1:length(ecc)
                    oneRow(numEcc) = subjectCond{ii}.(cond).(stab).(eccNames{numEcc}).threshInfo.thresh - ...
                        subjectCond{ii}.Uncrowded.(stab).(eccNames{1}).threshInfo.thresh;
                    oneRowEcc(numEcc) = eccVals(numEcc);
                    numTrials(ii,numEcc) = sum(subjectCond{ii}.(cond).(stab).(eccNames{numEcc}).threshInfo.em.valid);
            end
            
            hold on
            %             yyaxis left
            if bitget(counter,1) %odd
                leg(counterLeg) = plot(eccVals(~isnan(oneRow))+1,oneRow(~isnan(oneRow)),'-o','Color',c1(ii,:),...
                    'MarkerFace',c1(ii,:),'MarkerSize', 10, 'LineWidth',3);
                counterLeg = counterLeg + 1;
            else %even
                plot(eccVals(~isnan(oneRow))+1,oneRow(~isnan(oneRow)),'--o','Color',c2(ii,:),...
                    'MarkerFace',c2(ii,:),'MarkerSize', 10, 'LineWidth',3);
            end
            
            hold on
           
            counter = counter + 1;
        end
    end
    forLaterOneRow(ii,:) = oneRow;
end
hold on
ylabel('Threshold (arcmin)');
line([0 30],[0 0], 'Color','k','LineStyle','--');
set(gca,'xtick',eccVals,'xticklabel',string(eccVals), 'FontSize',14)
xlabel('Eccentricity (arcmin)')
xtickangle(60)
legend(leg, subjectsAll,...
    'Location','northwest');
title('Normalized to 0 Uncrowded at 0');

saveas(gcf, '../../../CrowdingTopology/Documents/Overleaf_CrowdingTopology/figures/acrossTheFovea/NormalizedCr.png');

%% Difference between Crowded and Uncrowded
figure('units','normalized','outerposition',[0.2 0.05 .4 .6]) %plots all absolute thresholds for each condition
counter = 1;
counterLeg = 1;
for ii = 1:length(subjectsAll)
    clear oneRow
    hold on
    for sIdx = 1:length(stabilization)
        clear allValsSubj
        cond = conditions{cIdx};
        stab = stabilization{sIdx};
        for numEcc = 1:length(ecc)
            oneRow(numEcc) = subjectCond{ii}.Crowded.(stab).(eccNames{numEcc}).threshInfo.thresh - ...
                subjectCond{ii}.Uncrowded.(stab).(eccNames{numEcc}).threshInfo.thresh;
            oneRowEcc(numEcc) = eccVals(numEcc);
        end
        forLaterOneRow(ii,:) = oneRow;
        hold on
        leg2(counterLeg) = plot(eccVals(~isnan(oneRow)),oneRow(~isnan(oneRow)),'o','Color',[.7 .7 .7],...
            'MarkerFace',[.7 .7 .7],'MarkerSize', 10, 'LineWidth',3);
        counterLeg = counterLeg + 1;
        hold on
        counter = counter + 1;
    end
end
hold on
ylabel({'Crowded-Uncrowded','\Delta Acuity Threshold (arcmin)'});

set(gca,'xtick',eccVals,'xticklabel',string(eccVals), 'FontSize',14)
xlabel('Eccentricity (arcmin)')
% xtickangle(60)

c = brewermap(12,'Paired');
hold on
coefficients = polyfit([0 10 15 25],mean(forLaterOneRow), 1);
xFit = linspace(0, 25, 1000);
yFit = polyval(coefficients , xFit);
m = (yFit(end)-yFit(1))/(xFit(end)-xFit(1)) ; %same as coefficents(1)

temp = plot(xFit, yFit,'--','Color','k',...
    'MarkerSize',12,'MarkerFaceColor','k','LineWidth',3);
hold on
errorbar([0 10 15 25],mean(forLaterOneRow),sem(forLaterOneRow),'o',...
    'Color','black','MarkerSize',12,'MarkerFaceColor','k');
ylim([0 2])
xlim([-1 26])

% yyaxis right
% plot([0 10 15 25],mean(forLaterOneRow)*20,'o','Color','k');
% ylim([5 35])
% ylabel('Snellen Acuity Equivalent')
% set(gca, 'YDir','reverse')
% legend([leg2 temp], [subjectsAll {sprintf('Best Fit Line, m=%.2f', m)}],...
%     'Location','northwest');

saveas(gcf, '../../../CrowdingTopology/Documents/Overleaf_CrowdingTopology/figures/acrossTheFovea/DeltaDiff.png');
saveas(gcf, '../../../CrowdingTopology/Documents/Overleaf_CrowdingTopology/figures/acrossTheFovea/DeltaDiff.epsc');

for ii = 1:length(subjectsAll)
    IsolatedThresh25(1,ii) = subjectCond{ii ...
        }.Uncrowded.Unstabilized.Side25ecc.threshInfo.thresh;
    IsolatedThresh0(1,ii) = subjectCond{ii ...
        }.Uncrowded.Unstabilized.Center0ecc.threshInfo.thresh;
end

mean(forLaterOneRow(:,1)./IsolatedThresh0')

mean(forLaterOneRow(:,4)./IsolatedThresh25')

figure;
[~,pU,bC,rU] = LinRegression([0 10 15 25],...
    mean(forLaterOneRow),...
    0,NaN,1,0);
%% Critical Spacing Specifically
figure('units','normalized','outerposition',[0.2 0.05 .4 .6]) %plots all absolute thresholds for each condition
counter = 1;
counterLeg = 1;
for ii = 1:length(subjectsAll)
    clear oneRow
    hold on
    for sIdx = 1:length(stabilization)
        clear allValsSubj
        cond = conditions{cIdx};
        stab = stabilization{sIdx};
        for numEcc = 1:length(ecc)
            oneRow(numEcc) = subjectCond{ii}.Crowded.(stab).(eccNames{numEcc}).threshInfo.thresh * 1.4;
            oneRowEcc(numEcc) = eccVals(numEcc);
        end
        forLaterOneRow(ii,:) = oneRow;
        hold on
        leg2(counterLeg) = plot(eccVals(~isnan(oneRow)),oneRow(~isnan(oneRow)),'o','Color', [.7 .7 .7],...
            'MarkerFace', [.7 .7 .7],'MarkerSize', 10, 'LineWidth',3);
        counterLeg = counterLeg + 1;
        hold on
        counter = counter + 1;
    end
end
hold on
ylabel('Critical Spacing (arcmin)');

set(gca,'xtick',eccVals,'xticklabel',string(eccVals), 'FontSize',14)
xlabel('Eccentricity (arcmin)')
% xtickangle(60)

c = brewermap(12,'Paired');
% title('Difference of Crowded - Uncrowded');
hold on
coefficients = polyfit([0 10 15 25],mean(forLaterOneRow), 1);
xFit = linspace(0, 25, 1000);
yFit = polyval(coefficients , xFit);
m = (yFit(end)-yFit(1))/(xFit(end)-xFit(1)) ; %same as coefficents(1)

temp = plot(xFit, yFit,'--','Color','k',...
    'MarkerSize',12,'MarkerFaceColor','k','LineWidth',3);
hold on
errorbar([0 10 15 25],mean(forLaterOneRow),sem(forLaterOneRow),'o','Color','black',...
    'MarkerSize',12,'MarkerFaceColor','k');

% legend([leg2 temp], [subjectsAll {sprintf('Best Fit Line, m=%.2f', m)}],...
%     'Location','northwest');
ylim([2.25 4.75])
xlim([-1 26])
saveas(gcf, '../../../CrowdingTopology/Documents/Overleaf_CrowdingTopology/figures/acrossTheFovea/CritSpac.png');
saveas(gcf, '../../../CrowdingTopology/Documents/Overleaf_CrowdingTopology/figures/acrossTheFovea/CritSpac.epsc');

figure;
[~,pU,bC,rU] = LinRegression([0 10 15 25],...
    mean(forLaterOneRow),...
    0,NaN,1,0);

%% Plot for Stabilized

% figure;
counter = 1;
for ii = find(subWithStab)
    for e = 1:4
%         for c = 1:2
            threshUStabCmp(counter,e) = log10(subjectCond{ii}.Crowded.Unstabilized.(...
                eccNames{e}).threshInfo.thresh/2);
            threshStabCmp(counter,e) = log10(subjectCond{ii}.Crowded.Stabilized.(...
                eccNames{e}).threshInfo.thresh/2);
            threshUStabCmpU(counter,e) = log10(subjectCond{ii}.Uncrowded.Unstabilized.(...
                eccNames{e}).threshInfo.thresh/2);
            threshStabCmpU(counter,e) = log10(subjectCond{ii}.Uncrowded.Stabilized.(...
                eccNames{e}).threshInfo.thresh/2);
            
%         end
    end
    counter = counter + 1;
end

figure;
plot(eccVals,threshStabCmp(1,:),'o','MarkerFaceColor',[.7 .7 .7],'MarkerEdgeColor',[.7 .7 .7])
hold on
plot(eccVals,threshStabCmp(2,:),'o','MarkerFaceColor',[.7 .7 .7],'MarkerEdgeColor',[.7 .7 .7])

plot(eccVals,threshStabCmpU(1,:),'o','MarkerFaceColor',[.7 .7 .7],'MarkerEdgeColor',[.7 .7 .7])
plot(eccVals,threshStabCmpU(2,:),'o','MarkerFaceColor',[.7 .7 .7],'MarkerEdgeColor',[.7 .7 .7])

clear forlegs
forlegs(1) = errorbar(eccVals,mean(threshStabCmp(:,:)),sem(threshStabCmp(:,:)),'o',...
    'MarkerFaceColor','r','MarkerEdgeColor','r','Color','r');
forlegs(2) =errorbar(eccVals,mean(threshStabCmpU(:,:)),sem(threshStabCmp(:,:)),'o',...
    'MarkerFaceColor','b','MarkerEdgeColor','b','Color','b');

md1 = fitlm([0 10 15 25],...
    mean(threshStabCmp));
x = 0:25; % Defines the domain as [-15,25] with a refinement of 0.25
m  = md1.Coefficients.Estimate(2);  % Specify your slope
x1 = 0; % Specify your starting x
y1 = mean(threshStabCmp(:,1));  % Specify your starting y
y = m*(x - x1) + y1;
hPlot = plot(x,y,'--r');

md1 = fitlm([0 10 15 25],...
    mean(threshStabCmpU));
x = 0:25; % Defines the domain as [-15,25] with a refinement of 0.25
m  = md1.Coefficients.Estimate(2);  % Specify your slope
x1 = 0; % Specify your starting x
y1 = mean(threshStabCmpU(:,1));  % Specify your starting y
y = m*(x - x1) + y1;
hPlot = plot(x,y,'--b');
legend(forlegs,{'Crowded','Uncrowded'},'Location','northwest');

xlabel('Eccentricity (arcmin)');
ylabel('Acuity (logMAR)');
xticks([0 10 15 25])
xlim([-2 28])
title('Stabilized');
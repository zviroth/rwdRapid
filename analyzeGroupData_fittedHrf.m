dataFolder = '~/data/rwdFmri/';
load([dataFolder 'rwdRapidData_fittedHrf.mat'], 'dataFolder', 'subFolders', 'samplerate', 'roiNames', ...
    'frames', 'junkedFrames', 'TR', 'trialsPerRun', 'trialLength', 'nVolumes',...
    'numContrasts', 'numFreqs', 'numSubs', 'numRois', 'onlyCorrect', ...
    'concatInfo', 'contrastBetasWithTask', 'sfBetasWithTask', 'sfContrastBetasWithTask', ...
    'taskResp','taskBaseline', 'nullBetas', 'contrastBetas','sfBetas','sfContrastBetas', ...
    'trialTC','meanDeconvNull','meanDeconvNullConstant','meanDeconvFreqConstant','nullDeconvBetas',...
    'meanDeconvContrastConstant', 'meanDeconvContrast', ...
    'meanDeconvFreqConstant', 'meanDeconvFreq',...
    'meanDeconvContrastFreq', 'meanDeconvContrastFreqConstant');
rwdString = {'high','low'};
plotColors = {[0 0 1], [1 0 0], [0 1 0], [0.5 1 0.2], [1 0.2 0.5] };
plotStyles = {'-','--',':','-.','-','--',':','-.'};
% subFolders = {'007420190221/'};%nice pupil response
% subFolders = {'007520190318/'};%almost no pupil difference
% subFolders = {'008320190320/'};%nice pupil response
% subFolders = {'009020190325/'};%no pupil size difference
% subFolders = {'001920190403/'};%very bad eye tracking
%%
numRois=2;
clear stdResp
for iSub = 1:numSubs
    subString{iSub} = subFolders{iSub}(1:4);
end

i=0;

%average fMRI response for high and low
i=i+1;
figure(i) 
clf
rows=numRois;
cols=numSubs;
subSumTrialTC=zeros(numRois,2,trialLength);
for iSub = 1:numSubs
    for iRoi=1:numRois
        subplot(rows,cols,iSub+(iRoi-1)*cols)
        for rwd=1:2
            temp = trialTC{iSub,iRoi}(rwd,:)';
            subSumTrialTC(iRoi,rwd,:) = squeeze(subSumTrialTC(iRoi,rwd,:)) + temp;
            plot(temp, 'Color',plotColors{rwd},'linewidth', 1); hold on
        end
        xlabel('time (TR)');
        title([subString{iSub} ' ' roiNames{iRoi}])
        if iSub==1; legend('avg high', 'avg low'); end
        stdResp(iSub,iRoi,:) = std(trialTC{iSub,iRoi}');
    end
end
subSumTrialTC = subSumTrialTC./numSubs;

%deconvolved fMRI response for null and stim trials, separately for high
%and low reward runs
i=i+1;
figure(i) %fMRI
clf
rows=numRois;
cols=numSubs;
subSumDeconvNull=zeros(numRois,2,2,trialLength);
for iSub = 1:numSubs
    for iRoi=1:numRois
        subplot(rows,cols,iSub+(iRoi-1)*cols)
        for rwd=1:2
            for c=1:2%stim,null
                %add constant to deconvolved response. constant differs
                %between reward type
%                 temp = squeeze(meanDeconvNull(iSub,iRoi, rwd,c,:) + meanDeconvNullConstant(iSub,iRoi, rwd));
                temp = squeeze(meanDeconvNull(iSub,iRoi, rwd,c,:));
                subSumDeconvNull(iRoi,rwd,c,:) = squeeze(subSumDeconvNull(iRoi,rwd,c,:)) + temp;%summing across subjects
                plot(temp, plotStyles{c}, 'Color',plotColors{rwd},'linewidth', 1); hold on
            end
%             plot(squeeze(meanDeconvNull{iSub}(iRoi, rwd,:,:))', plotStyles{rwd}, 'linewidth', 1); hold on
        end
        xlabel('time (TR)');
        title([subString{iSub} ' ' roiNames{iRoi}])
        if iSub==1 & iRoi==1; legend('high stim','high null','low stim', 'low null'); end
    end
end
subSumDeconvNull = subSumDeconvNull./numSubs;

%null/sttim betas
i=i+1;
figure(i)
clf
rows=numRois;
cols=numSubs;
subSumNullBetas=zeros(numRois,2,2);
for iSub = 1:numSubs
    for iRoi=1:numRois
        subplot(rows,cols,iSub+(iRoi-1)*cols)
        for rwd=1:2
            temp = nanmean(nullBetas{iSub,iRoi,rwd}(1:end-1,:),2);
%             temp = temp+repmat(nanmean(nullBetas{iSub,iRoi,rwd}(end,:)),2,1);%adding the constant beta ??
            subSumNullBetas(iRoi,rwd,:) = squeeze(subSumNullBetas(iRoi,rwd,:)) + temp;
            plot(temp, 'Color',plotColors{rwd},'linewidth', 1); hold on
        end
%         xlabel('contrast');
        xticks([ 1 2]);
        xticklabels({'stim', 'null'});
        title([subString{iSub} ' ' roiNames{iRoi}])
        if iSub==1 & iRoi==1; legend('high','low'); end
    end
end
subSumNullBetas = subSumNullBetas./numSubs;


% contrast timecourse without removing task response
i=i+1;
figure(i) %fMRI
clf
rows=numRois;
cols=numSubs;
subSumDeconvContrast=zeros(numRois,2,numContrasts,trialLength);
for iSub = 1:numSubs
    for iRoi=1:numRois
        subplot(rows,cols,iSub+(iRoi-1)*cols)
        for rwd=1:2
            for c=1:numContrasts% low contrast, high contrast
                temp = squeeze(meanDeconvContrast(iSub,iRoi, rwd,c,:) + meanDeconvContrastConstant(iSub,iRoi, rwd));%adding constant beta
                temp = squeeze(meanDeconvContrast(iSub,iRoi, rwd,c,:));
                subSumDeconvContrast(iRoi,rwd,c,:) = squeeze(subSumDeconvContrast(iRoi,rwd,c,:)) + temp;%summing across subjects
                plot(temp, plotStyles{c}, 'Color',plotColors{rwd},'linewidth', 1); hold on
            end
        end
        xlabel('time (TR)');
        if iSub==1 & iRoi==1; legend('high rwd low contrast','high rwd high contrast','low rwd low contrast', 'low rwd high contrast'); end
    end
end
subSumDeconvContrast = subSumDeconvContrast./numSubs;

%contrast betas without removing task response
i=i+1;
figure(i)
clf
rows=numRois;
cols=numSubs;
subSumContrastBetasWithTask=zeros(numRois,2,numContrasts);
for iSub = 1:numSubs
    for iRoi=1:numRois
        subplot(rows,cols,iSub+(iRoi-1)*cols)
        for rwd=1:2
            temp = nanmean(contrastBetasWithTask{iSub,iRoi,rwd}(1:end-1,:),2);%contrastBetasWithTask{iRoi,rwdTypeNum}(:,i)
%             temp = temp+repmat(nanmean(contrastBetasWithTask{iSub,iRoi,rwd}(end,:)),numContrasts,1);%adding the constant beta ??
            subSumContrastBetasWithTask(iRoi,rwd,:) = squeeze(subSumContrastBetasWithTask(iRoi,rwd,:)) + temp;
            plot(temp, 'Color',plotColors{rwd},'linewidth', 1); hold on
        end
        xlabel('contrast');
        xticks([ 1 2]);
        xticklabels({'low', 'high'});
        title([subString{iSub} ' ' roiNames{iRoi}])
        if iSub==1 & iRoi==1; legend('high rwd','low rwd'); end
    end
end
subSumContrastBetasWithTask = subSumContrastBetasWithTask./numSubs;

%contrast betas after removing task response
i=i+1;
figure(i)
clf
rows=numRois;
cols=numSubs;
subSumContrastBetas=zeros(numRois,2,numContrasts);
for iSub = 1:numSubs
    for iRoi=1:numRois
        subplot(rows,cols,iSub+(iRoi-1)*cols)
        for rwd=1:2
            temp = nanmean(contrastBetas{iSub,iRoi,rwd}(1:end-1,:),2);
%             temp = temp+repmat(nanmean(contrastBetas{iSub,iRoi,rwd}(end,:)),numContrasts,1);%adding the constant beta ??
            subSumContrastBetas(iRoi,rwd,:) = squeeze(subSumContrastBetas(iRoi,rwd,:)) + temp;
            plot(temp, 'Color',plotColors{rwd},'linewidth', 1); hold on
        end
        xlabel('time (TR)');
        title([subString{iSub} ' ' roiNames{iRoi}])
        if iSub==1 & iRoi==1; legend('high','low'); end
    end
end
subSumContrastBetas = subSumContrastBetas./numSubs;


% spatial frequency timecourse without removing task response
i=i+1;
figure(i) %fMRI
clf
rows=numRois;
cols=numSubs;
subSumDeconvFreq=zeros(numRois,2,numFreqs,trialLength);
for iSub = 1:numSubs
    for iRoi=1:numRois
        subplot(rows,cols,iSub+(iRoi-1)*cols)
        for rwd=1:2
            for c=1:numFreqs% low contrast, high contrast
%                 temp = squeeze(meanDeconvFreq(iSub,iRoi, rwd,c,:) + meanDeconvFreqConstant(iSub,iRoi, rwd));%add constant beta
                temp = squeeze(meanDeconvFreq(iSub,iRoi, rwd,c,:));
                subSumDeconvFreq(iRoi,rwd,c,:) = squeeze(subSumDeconvFreq(iRoi,rwd,c,:)) + temp;%summing across subjects
                plot(temp, 'Color',plotColors{rwd},'linewidth', 1); hold on
%                 plot(squeeze(meanDeconvFreq(iSub,iRoi, rwd,c,:)), plotStyles{c}, 'Color',plotColors{rwd},'linewidth', 1); hold on
            end
        end
        if iSub==1 & iRoi==1; legend('high','low'); end
    end
end
subSumDeconvFreq = subSumDeconvFreq./numSubs;



%spatial frequency betas without removing task response
i=i+1;
figure(i)
clf
rows=numRois;
cols=numSubs;
subSumSfBetasWithTask=zeros(numRois,2,numFreqs);
for iSub = 1:numSubs
    for iRoi=1:numRois
        subplot(rows,cols,iSub+(iRoi-1)*cols)
        for rwd=1:2
            temp = nanmean(sfBetasWithTask{iSub,iRoi,rwd}(1:end-1,:),2);
%             temp = temp+repmat(nanmean(sfBetasWithTask{iSub,iRoi,rwd}(end,:)),numFreqs,1);%adding the constant beta ??
            subSumSfBetasWithTask(iRoi,rwd,:) = squeeze(subSumSfBetasWithTask(iRoi,rwd,:)) + temp;
            plot(temp, 'Color',plotColors{rwd},'linewidth', 1); hold on
        end
        xlabel('time (TR)');
        title([subString{iSub} ' ' roiNames{iRoi}])
        if iSub==1 & iRoi==1; legend('high','low'); end
    end
end
subSumSfBetasWithTask = subSumSfBetasWithTask./numSubs;

% %spatial frequency betas after removing task response
% i=i+1;
% figure(i)
% clf
% rows=numRois;
% cols=numSubs;
% subSumSfBetas=zeros(numRois,2,numFreqs);
% for iSub = 1:numSubs
%     for iRoi=1:numRois
%         subplot(rows,cols,iSub+(iRoi-1)*cols)
%         for rwd=1:2
%             temp = nanmean(sfBetas{iSub,iRoi,rwd}(1:end-1,:),2);
% %             temp = temp+repmat(nanmean(sfBetas{iSub,iRoi,rwd}(end,:)),numFreqs,1);%adding the constant beta ??
%             subSumSfBetas(iRoi,rwd,:) = squeeze(subSumSfBetas(iRoi,rwd,:)) + temp;
%             plot(temp, 'Color',plotColors{rwd},'linewidth', 1); hold on
%         end
%         xlabel('time (TR)');
%         title([subString{iSub} ' ' roiNames{iRoi}])
%         if iSub==1 & iRoi==1; legend('high','low'); end
%     end
% end
% subSumSfBetas = subSumSfBetas./numSubs;

% spatial frequency & contrast timecourse without removing task response
i=i+1;
figure(i) 
clf
rows=numRois;
cols=numSubs;
subSumDeconvContrastFreq=zeros(numRois,2,numFreqs*numContrasts,trialLength);
for iSub = 1:numSubs
    for iRoi=1:numRois
        subplot(rows,cols,iSub+(iRoi-1)*cols)
        for rwd=1:2
            for c=1:numFreqs*numContrasts
                temp = squeeze(meanDeconvContrastFreq(iSub,iRoi, rwd,c,:));
                subSumDeconvContrastFreq(iRoi,rwd,c,:) = squeeze(subSumDeconvContrastFreq(iRoi,rwd,c,:)) + temp;%summing across subjects
                plot(temp, 'Color',plotColors{rwd},'linewidth', 1); hold on
            end
        end
        if iSub==1 & iRoi==1; legend('high','low'); end
    end
end
subSumDeconvContrastFreq = subSumDeconvContrastFreq./numSubs;


%spatial frequency & contrast betas without removing task response
i=i+1;
figure(i)
clf
rows=numRois;
cols=numSubs;
subSumSfContrastBetasWithTask=zeros(numRois,2,numFreqs*numContrasts);
subRdm=zeros(numSubs,numRois,2,nchoosek(numContrasts*numFreqs,2));
for iSub = 1:numSubs
    for iRoi=1:numRois
        subplot(rows,cols,iSub+(iRoi-1)*cols)
        for rwd=1:2
            temp = nanmean(sfContrastBetasWithTask{iSub,iRoi,rwd}(1:end-1,:),2);
%             temp = temp+repmat(nanmean(sfContrastBetasWithTask{iSub,iRoi,rwd}(end,:)),numFreqs,1);%adding the constant beta ??
            subSumSfContrastBetasWithTask(iRoi,rwd,:) = squeeze(subSumSfContrastBetasWithTask(iRoi,rwd,:)) + temp;
            plot([temp(1:5); NaN; temp(6:end)], 'Color',plotColors{rwd},'linewidth', 1); hold on
            
            goodVox = ~isnan(sum(sfContrastBetasWithTask{iSub,iRoi,rwd}(1:end-1,:)));
             corMat = corr(sfContrastBetasWithTask{iSub,iRoi,rwd}(1:end-1,goodVox)');
             subRdm(iSub,iRoi,rwd,:) = squareform(1-corMat);
        end
        xlabel('time (TR)');
        title([subString{iSub} ' ' roiNames{iRoi}])
        if iSub==1 & iRoi==1; legend('high','low'); end
    end
end
subSumSfContrastBetasWithTask = subSumSfContrastBetasWithTask./numSubs;

%spatial frequency & contrast betas after removing task response
i=i+1;
figure(i)
clf
rows=numRois;
cols=numSubs;
subSumSfContrastBetas=zeros(numRois,2,numFreqs*numContrasts);
for iSub = 1:numSubs
    for iRoi=1:numRois
        subplot(rows,cols,iSub+(iRoi-1)*cols)
        for rwd=1:2
            temp = nanmean(sfContrastBetas{iSub,iRoi,rwd}(1:end-1,:),2);
%             temp = temp+repmat(nanmean(sfContrastBetasWithTask{iSub,iRoi,rwd}(end,:)),numFreqs,1);%adding the constant beta ??
            subSumSfContrastBetas(iRoi,rwd,:) = squeeze(subSumSfContrastBetas(iRoi,rwd,:)) + temp;
            plot([temp(1:5); NaN; temp(6:end)], 'Color',plotColors{rwd},'linewidth', 1); hold on
        end
        xlabel('time (TR)');
        title([subString{iSub} ' ' roiNames{iRoi}])
        if iSub==1 & iRoi==1; legend('high','low'); end
    end
end
subSumSfContrastBetas = subSumSfContrastBetas./numSubs;

%% plot RDMs
roiRdm = squeeze(mean(subRdm));%roiRdm(iRoi,rwd,:);
i=i+1;
figure(i)
clf
rows=2;
cols=numRois;
for rwd=1:2
    for iRoi=1:numRois
        subplot(rows,cols,iRoi+ (rwd-1)*cols)
        imagesc(squareform(squeeze(roiRdm(iRoi,rwd,:))));

        title([roiNames{iRoi} ' ' rwdString{rwd}]);
    end
end


%%
%plot averages
i=i+1;
figure(i)
clf
rows=numRois;
cols=10;
j=1;
for iRoi=1:numRois
    subplot(rows,cols,j+(iRoi-1)*cols)
    for rwd=1:2
            plot(squeeze(subSumTrialTC(iRoi,rwd,:)), 'Color',plotColors{rwd},'linewidth', 1); hold on
    end
    title(roiNames{iRoi})
    xlabel('time (TR)');
    if iRoi==1; legend('high','low'); end
end
%null and stim
j=j+1;
for iRoi=1:numRois
    subplot(rows,cols,j+(iRoi-1)*cols)
    for rwd=1:2
        for c=1:2
            plot(squeeze(subSumDeconvNull(iRoi,rwd,c,:)), plotStyles{c},'Color',plotColors{rwd},'linewidth', 1); hold on
        end
    end
    xlabel('time (TR)');
    title(roiNames{iRoi})
    if iRoi==1; legend('high stim','high null','low stim', 'low null'); end
end

% null minus stim
j=j+1;
for iRoi=1:numRois
    subplot(rows,cols,j+(iRoi-1)*cols)
    for rwd=1:2
        plot(squeeze(subSumDeconvNull(iRoi,rwd,1,:)) - squeeze(subSumDeconvNull(iRoi,rwd,2,:)), 'Color',plotColors{rwd},'linewidth', 1); hold on
    end
    xlabel('time (TR)');
    title(roiNames{iRoi})
    if iRoi==1; legend('high stim-null','low stim-null'); end
end

%subSumNullBetas
j=j+1;
for iRoi=1:numRois
    subplot(rows,cols,j+(iRoi-1)*cols)
    for rwd=1:2
            plot(squeeze(subSumNullBetas(iRoi,rwd,:)), plotStyles{1},'Color',plotColors{rwd},'linewidth', 1); hold on
%             plot(squeeze(subSumContrastBetas(iRoi,rwd,:)), plotStyles{2},'Color',plotColors{rwd},'linewidth', 1); hold on
    end
    title(roiNames{iRoi})
%     xlabel('contrast');
    xticks([ 1 2]);
    xticklabels({'stim','null'});
    if iRoi==1; legend('high rwd','low rwd'); end
end

%contrast timecourse
j=j+1;
for iRoi=1:numRois
    subplot(rows,cols,j+(iRoi-1)*cols)
    for rwd=1:2
        for c=1:numContrasts
            plot(squeeze(subSumDeconvContrast(iRoi,rwd,c,:)), plotStyles{c},'Color',plotColors{rwd},'linewidth', 1); hold on
        end
    end
    xlabel('time (TR)');
    title(roiNames{iRoi})
    if iRoi==1; legend('high rwd low contrast','high rwd high contrast','low rwd low contrast', 'low rwd high contrast'); end
end
% high contrast minus low contrast
j=j+1;
for iRoi=1:numRois
    subplot(rows,cols,j+(iRoi-1)*cols)
    for rwd=1:2
        plot(squeeze(subSumDeconvContrast(iRoi,rwd,2,:)) - squeeze(subSumDeconvContrast(iRoi,rwd,1,:)), 'Color',plotColors{rwd},'linewidth', 1); hold on
    end
    xlabel('time (TR)');
    title(roiNames{iRoi})
    if iRoi==1; legend('high high-low contrast','low high-low contrast'); end
end
%contrast betas
j=j+1;
for iRoi=1:numRois
    subplot(rows,cols,j+(iRoi-1)*cols)
    for rwd=1:2
            plot(squeeze(subSumContrastBetasWithTask(iRoi,rwd,:)), plotStyles{1},'Color',plotColors{rwd},'linewidth', 1); hold on
%             plot(squeeze(subSumContrastBetas(iRoi,rwd,:)), plotStyles{2},'Color',plotColors{rwd},'linewidth', 1); hold on
    end
    title(roiNames{iRoi})
    xlabel('contrast');
    xticks([ 1 2]);
    xticklabels({'low','high'});
    if iRoi==1; legend('high rwd','low rwd'); end
end
%frequency timecourse
j=j+1;
for iRoi=1:numRois
    subplot(rows,cols,j+(iRoi-1)*cols)
    for rwd=1:2
        for c=1:numFreqs
            plot(squeeze(subSumDeconvFreq(iRoi,rwd,c,:)), 'Color',plotColors{rwd},'linewidth', 1); hold on
        end
    end
    xlabel('time (TR)');
    title(roiNames{iRoi})
%     if iRoi==1; legend('high rwd low contrast','high rwd high contrast','low rwd low contrast', 'low rwd high contrast'); end
end
% freq betas
j=j+1;
for iRoi=1:numRois
    subplot(rows,cols,j+(iRoi-1)*cols)
    for rwd=1:2
            plot(squeeze(subSumSfBetasWithTask(iRoi,rwd,:)), plotStyles{1},'Color',plotColors{rwd},'linewidth', 1); hold on
%             plot(squeeze(subSumSfBetas(iRoi,rwd,:)), plotStyles{2},'Color',plotColors{rwd},'linewidth', 1); hold on
    end
    title(roiNames{iRoi})
    xlabel('frequency');
%         xticks([ 1:5]);
%         xticklabels({'high','low'});
    if iRoi==1; legend('high rwd','low rwd'); end
end


%sf + contrast betas
j=j+1;
for iRoi=1:numRois
    subplot(rows,cols,j+(iRoi-1)*cols)
    for rwd=1:2
            plot([squeeze(subSumSfContrastBetasWithTask(iRoi,rwd,1:5)); NaN; squeeze(subSumSfContrastBetasWithTask(iRoi,rwd,6:end))], plotStyles{1},'Color',plotColors{rwd},'linewidth', 1); hold on
%             plot([squeeze(subSumSfContrastBetas(iRoi,rwd,1:5)); NaN; squeeze(subSumSfContrastBetas(iRoi,rwd,6:end))], plotStyles{2},'Color',plotColors{rwd},'linewidth', 1); hold on
    end
    title(roiNames{iRoi})
    if iRoi==1; legend('high rwd','low rwd'); end
end


% %sf + contrast with task removed
% j=j+1;
% for iRoi=1:numRois
%     subplot(rows,cols,j+(iRoi-1)*cols)
%     for rwd=1:2
% %             plot([squeeze(subSumSfContrastBetasWithTask(iRoi,rwd,1:5)); NaN; squeeze(subSumSfContrastBetasWithTask(iRoi,rwd,6:end))], plotStyles{1},'Color',plotColors{rwd},'linewidth', 1); hold on
%             plot([squeeze(subSumSfContrastBetas(iRoi,rwd,1:5)); NaN; squeeze(subSumSfContrastBetas(iRoi,rwd,6:end))], plotStyles{1},'Color',plotColors{rwd},'linewidth', 1); hold on
%     end
%     title(roiNames{iRoi})
%     if iRoi==1; legend('high rwd','low rwd'); end
% end
%%

%surface plot of average frequency time courses
i=i+1;
figure(i)
clf
rows = numRois;
cols=2;
X=repmat(1:trialLength, numFreqs,1);
Y=repmat(1:numFreqs, trialLength,1)';

for iRoi=1:numRois
    for rwd=1:2
        subplot(rows,cols,rwd+(iRoi-1)*cols)
        Z=squeeze(subSumDeconvFreq(iRoi,rwd,:,:));
        surf(X,Y,Z);
        title([roiNames{iRoi} ' ' rwdString{rwd} ' rwd']);
        xlabel('time (TR)');
        ylabel('frequency');
    end
end

%surface plot of average high-low contrast frequency time courses
i=i+1;
figure(i)
clf
rows = numRois;
cols=2;
X=repmat(1:trialLength, numFreqs,1);
Y=repmat(1:numFreqs, trialLength,1)';
for iRoi=1:numRois
    for rwd=1:2
        subplot(rows,cols,rwd+(iRoi-1)*cols)
        Z=squeeze(subSumDeconvContrastFreq(iRoi,rwd,numFreqs+1:2*numFreqs,:) - subSumDeconvContrastFreq(iRoi,rwd,1:numFreqs,:));
        surf(X,Y,Z);
        title([roiNames{iRoi} ' ' rwdString{rwd} ' rwd: high-low contrast' ]);
        xlabel('time (TR)');
        ylabel('frequency');
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = plotRwdStimSegments(meanPupil, allTimes)

itime = 1;
startT=1;
endT = startT+allTimes(itime);
plot(startT:endT-1,meanPupil(startT:endT-1),'lineWidth',4,'color',[0 1 1]);
itime = itime+1;
startT=endT;
endT = startT+allTimes(itime);
plot(startT:endT-1,meanPupil(startT:endT-1),'lineWidth',4,'color',[0.1 0.7 0.6]);
itime = itime+1;
startT=endT;
endT = startT+allTimes(itime);
plot(startT:endT-1,meanPupil(startT:endT-1),'lineWidth',4,'color',[0 1 1]);
itime = itime+1;
startT=endT;
endT = startT+allTimes(itime);
plot(startT:endT-1,meanPupil(startT:endT-1),'lineWidth',4,'color',[1 1 0.0]);

end
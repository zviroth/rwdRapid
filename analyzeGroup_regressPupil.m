nullTrialsPredictor=0;%either only null trials or all trials
regressGlobalMean = 1;
dataFolder = '~/data/rwdFmri/';
if nullTrialsPredictor
    nullTrialsString = 'onlyNulls';
else
    nullTrialsString = 'allTrials';
end
globalMeanString = '';
if regressGlobalMean
    globalMeanString = '_globalRegressed';
end
load([dataFolder 'rwdRapidData_regressPupil_' nullTrialsString globalMeanString '.mat'], 'dataFolder', 'subFolders', 'samplerate', 'roiNames', ...
    'frames', 'junkedFrames', 'TR', 'trialsPerRun', 'trialLength', 'nVolumes',...
    'numContrasts', 'numFreqs', 'numSubs', 'numRois', 'onlyCorrect', ...
    'concatInfo', 'nullBetas', 'contrastBetas','sfBetas','sfContrastBetas', ...
    'deconvLength',...
    'contrastBetas', 'sfBetas', 'sfContrastBetas', 'nullBetas', ...
    'nullDeconv', 'contrastDeconv', 'freqDeconv', 'contrastFreqDeconv', ...
    'binMeanTseries','binTrialTC','numBinVoxels','binHrf',...
    'nullDeconvStd', 'contrastDeconvStd', 'freqDeconvStd', 'contrastFreqDeconvStd',...
    'nbins','binBorders','trialTCstd',...
    'nullDeconvFftAmp','nullDeconvFftPh','contrastDeconvFftAmp','contrastDeconvFftPh',...
    'freqDeconvFftAmp','freqDeconvFftPh','contrastFreqDeconvFftAmp','contrastFreqDeconvFftPh','stimulus',...
    'globalMean','regressBetasGlobal','runRwd',...
    'physioSampleRate','sacRwd','binnedSac', 'smoothSac', ...
    'rwdPupil','meanPupil','taskTimes',...
    'sacRun','pupilRun', 'numRuns','startTimes','endTimes','eyetrackerTime',...
    'nullTrials','meanNullPupil','meanStimPupil',...
    'deconvTaskResp','stdPupil','nullPupilTrials');
rwdString = {'H','L'};
plotColors = { [1 0 0], [0 0 1], [0 1 0], [0.5 1 0.2], [1 0.2 0.5] };
plotStyles = {'.','-','--',':','-.','-','--',':','-.'};
linewidth = 2;
markersize=20;
errbarwidth = 1;
%%

for iSub = 1:numSubs
    subString{iSub} = subFolders{iSub}(1:4);
end

i=0;
rows=8;
cols=nbins;


i=i+1;
figure(i) 
clf

row=0;

%ASSUME ONLY 1 ROI
iRoi=1;
%average fMRI response for high and low
row=row+1;

meanBinTC = mean(binTrialTC,1);%average over subjects
for ibin=1:nbins
    subplot(rows,cols,ibin+(row-1)*cols)
    for rwd=1:2
        plot(squeeze(meanBinTC(:,:,ibin,rwd,:))','Color',plotColors{rwd},'linewidth',linewidth);
        hold on
    end
    xlabel('time (TR)');
    title([roiNames{iRoi} ' bin ' num2str(ibin)])
    if ibin==1; legend({'high rwd','low rwd'},'location','northwest'); end
end

%std of average fMRI response
row=row+1;
for ibin=1:nbins
    subplot(rows,cols,ibin+(row-1)*cols)
        plot(squeeze(mean(trialTCstd(:,ibin,:)))','Color','k','linewidth',linewidth);
        hold on
    xticks([ 1 2]);
    xticklabels({'high rwd', 'low rwd'});
    title([roiNames{iRoi} ' bin ' num2str(ibin) ' std'])
%     if ibin==1; legend({'high rwd','low rwd'},'location','northwest'); end
end


%deconvolved fMRI response for null/stim
row=row+1;
meanBinNullTC = mean(nullDeconv,1);%nullDeconv(iSub,iRoi,ibin,rwd,:,:)
for ibin=1:nbins
    subplot(rows,cols,ibin+(row-1)*cols)
    for rwd=1:2
        for c=1:size(meanBinNullTC,5)
        
            plot(squeeze(meanBinNullTC(:,:,ibin,rwd,c,:))',plotStyles{c+1},'Color',plotColors{rwd},'linewidth',linewidth);
            hold on
        end
    end
    xlabel('time (TR)');
    title([roiNames{iRoi} ' bin ' num2str(ibin)])
    if ibin==1; legend({'stim','null'},'location','northwest'); end
end

%std of deconvolved fMRI response for null/stim
row=row+1;
for ibin=1:nbins
    subplot(rows,cols,ibin+(row-1)*cols)
    for rwd=1:2
        plot(squeeze(mean(nullDeconvStd(:,ibin,rwd,:)))',plotStyles{2},'Color',plotColors{rwd},'linewidth',linewidth);
        hold on
    end
    xticks([ 1 2]);
    xticklabels({'stim', 'null'});
    title([roiNames{iRoi} ' bin ' num2str(ibin) ' std'])
    if ibin==1; legend({'high rwd','low rwd'},'location','northwest'); end
end

%betas for null/stim
row=row+1;
meanNullBetas = mean(nullBetas,1);%nullBetas(iSub,iRoi,ibin,rwd,:)
for ibin=1:nbins
    subplot(rows,cols,ibin+(row-1)*cols)
    for rwd=1:2
        plot(squeeze(meanNullBetas(:,:,ibin,rwd,:))',plotStyles{2},'Color',plotColors{rwd},'linewidth',linewidth);
        hold on
    end
    xticks([ 1 2]);
    xticklabels({'stim', 'null'});
    title([roiNames{iRoi} ' bin ' num2str(ibin) ' betas'])
    if ibin==1; legend({'high rwd','low rwd'},'location','northwest'); end
end


%deconvolved fMRI response for high/low contrast
row=row+1;
meanBinContrastTC = mean(contrastDeconv,1);%contrastDeconv(iSub,iRoi,ibin,rwd,:,:)
for ibin=1:nbins
    subplot(rows,cols,ibin+(row-1)*cols)
    for rwd=1:2
        for c=1:size(meanBinContrastTC,5)
            plot(squeeze(meanBinContrastTC(:,:,ibin,rwd,c,:))',plotStyles{1+c},'Color',plotColors{rwd},'linewidth',linewidth);
            hold on
        end
    end
    xlabel('time (TR)');
    title([roiNames{iRoi} ' bin ' num2str(ibin)])
    if ibin==1; legend({'low contrast','high contrast'},'location','northwest'); end
end


%std of deconvolved fMRI response for contrast
row=row+1;
for ibin=1:nbins
    subplot(rows,cols,ibin+(row-1)*cols)
    for rwd=1:2
        plot(squeeze(mean(contrastDeconvStd(:,ibin,rwd,:)))',plotStyles{2},'Color',plotColors{rwd},'linewidth',linewidth);
        hold on
    end
    xticks([ 1 2]);
    xticklabels({'low contrast', 'high contrast'});
    title([roiNames{iRoi} ' bin ' num2str(ibin) ' std'])
    if ibin==1; legend({'high rwd','low rwd'},'location','northwest'); end
end

%betas for contrast
row=row+1;
meanContrastBetas = mean(contrastBetas,1);%nullBetas(iSub,iRoi,ibin,rwd,:)
for ibin=1:nbins
    subplot(rows,cols,ibin+(row-1)*cols)
    for rwd=1:2
        plot(squeeze(meanContrastBetas(:,:,ibin,rwd,:))',plotStyles{2},'Color',plotColors{rwd},'linewidth',linewidth);
        hold on
    end
    xticks([ 1 2]);
    xticklabels({'low contrast', 'high contrast'});
    title([roiNames{iRoi} ' bin ' num2str(ibin) ' betas'])
    if ibin==1; legend({'high rwd','low rwd'},'location','northwest'); end
end




%% SF
i=i+1;
figure(i) 
clf
rows=6;

row=0;
%deconvolved fMRI response for SF
row=row+1;
meanFreqBinTC = mean(freqDeconv,1);%average over subjects
for ibin=1:nbins
    subplot(rows,cols,ibin+(row-1)*cols)
    for rwd=1:2
        for c=1:size(meanFreqBinTC,5)
            plot(squeeze(meanFreqBinTC(:,:,ibin,rwd,c,:))',plotStyles{c+1},'Color',plotColors{rwd},'linewidth',linewidth);
            hold on
        end
    end
    xlabel('time (TR)');
    if ibin==1; ylabel('BOLD (std)'); end
    title([num2str(binBorders(ibin))  '\circ - ' num2str(binBorders(ibin+1)) '\circ']);
%     title([roiNames{iRoi} ' bin ' num2str(ibin)])
%     if ibin==1; legend({'low contrast','high contrast'},'location','northwest'); end
end

%std of deconvolved fMRI response for SF
row=row+1;
ymin = min(freqDeconvStd(:));
ymax = max(freqDeconvStd(:));
for ibin=1:nbins
    subplot(rows,cols,ibin+(row-1)*cols)
    for rwd=1:2
%         plot(1:numFreqs,squeeze(mean(freqDeconvStd(:,ibin,rwd,1:end-1)))',plotStyles{1},'Color',plotColors{rwd},'linewidth',linewidth,'markersize',markersize);%not plotting for "task"
        errorbar(1:numFreqs, squeeze(mean(freqDeconvStd(:,ibin,rwd,:)))',squeeze(std(freqDeconvStd(:,ibin,rwd,:)))'./sqrt(numSubs),plotStyles{1},'Color',plotColors{rwd},'linewidth',errbarwidth,'markersize',markersize);
        hold on
    end
    xticks([]);
%     title([roiNames{iRoi} ' bin ' num2str(ibin) ' std'])
    if ibin==1; legend({'high rwd','low rwd'},'location','northwest'); ylabel('STD'); end
    ylim([ymin ymax]);
end

%std of deconvolved fMRI response for SF
row=row+1;
ymin = min(freqDeconvStd(:));
ymax = max(freqDeconvStd(:));
for ibin=1:nbins
    subplot(rows,cols,ibin+(row-1)*cols)
    for rwd=1:2
%         plot(1:numFreqs,squeeze(mean(freqDeconvStd(:,ibin,rwd,1:end-1)))',plotStyles{1},'Color',plotColors{rwd},'linewidth',linewidth,'markersize',markersize);%not plotting for "task"
        errorbar(1:numFreqs, squeeze(mean(freqDeconvStd(:,ibin,rwd,:)))',squeeze(std(freqDeconvStd(:,ibin,rwd,:)))'./sqrt(numSubs),plotStyles{1},'Color',plotColors{rwd},'linewidth',errbarwidth,'markersize',markersize);
        hold on
    end
    xticks([]);
%     title([roiNames{iRoi} ' bin ' num2str(ibin) ' std'])
    if ibin==1; legend({'high rwd','low rwd'},'location','northwest'); ylabel('STD'); end
%     ylim([ymin ymax]);
end

% 
%betas for SF
row=row+1;
meanSfBetas = mean(sfBetas,1);%nullBetas(iSub,iRoi,ibin,rwd,:)
stdSfBetas = std(sfBetas);
ymin = min(meanSfBetas(:));
ymax = max(meanSfBetas(:));
for ibin=1:nbins
    subplot(rows,cols,ibin+(row-1)*cols)
    for rwd=1:2
%         plot(1:numFreqs,squeeze(meanSfBetas(:,:,ibin,rwd,1:end-1))',plotStyles{1},'Color',plotColors{rwd},'linewidth',linewidth,'markersize',markersize);%not plotting for "task"
        errorbar(1:numFreqs,squeeze(meanSfBetas(:,:,ibin,rwd,:))',squeeze(stdSfBetas(:,:,ibin,rwd,:))'./sqrt(numSubs),plotStyles{1},'Color',plotColors{rwd},'linewidth',errbarwidth,'markersize',markersize);%not plotting for "task"
        hold on
    end
    xticks([]);
%     title([roiNames{iRoi} ' bin ' num2str(ibin) ' betas'])
    if ibin==1; legend({'high rwd','low rwd'},'location','northwest'); ylabel('betas'); end
%     ylim([ymin ymax]);
end

%FFT amp for SF
row=row+1;
meanFreqFftAmp = mean(freqDeconvFftAmp,1);%nullBetas(iSub,iRoi,ibin,rwd,:)
stdFreqFftAmp = std(freqDeconvFftAmp);
ymin = min(meanFreqFftAmp(:));
ymax = max(meanFreqFftAmp(:));
for ibin=1:nbins
    subplot(rows,cols,ibin+(row-1)*cols)
    for rwd=1:2
%         plot(1:numFreqs,squeeze(meanFreqFftAmp(:,ibin,rwd,1:end-1))',plotStyles{1},'Color',plotColors{rwd},'linewidth',linewidth,'markersize',markersize);%not plotting for "task"
        errorbar(1:numFreqs,squeeze(meanFreqFftAmp(:,ibin,rwd,:))',squeeze(stdFreqFftAmp(:,ibin,rwd,:))'./sqrt(numSubs),plotStyles{1},'Color',plotColors{rwd},'linewidth',errbarwidth,'markersize',markersize);%not plotting for "task"
        hold on
    end
    xticks(1:numFreqs);
    xticklabels(num2str(stimulus.freqs','%0.1f'));
%     xlabel('frequency (cpd)');
%     title([roiNames{iRoi} ' bin ' num2str(ibin) ' FFT'])
    if ibin==1; legend({'high rwd','low rwd'},'location','northwest'); ylabel('fft amplitude'); end
%     ylim([ymin ymax]);
end

%average MAX SF
row=row+1;
ymin = min(freqDeconv(:));
ymax = max(freqDeconv(:));
for ibin=1:nbins
    subplot(rows,cols,ibin+(row-1)*cols)
    for rwd=1:2
%         plot(1:numFreqs,mean(max(squeeze(freqDeconv(:,iRoi,ibin,rwd,1:end-1,:)),[],3)),plotStyles{1},'Color',plotColors{rwd},'linewidth',linewidth,'markersize',markersize);%not plotting for "task"
        errorbar(1:numFreqs,mean(max(squeeze(freqDeconv(:,iRoi,ibin,rwd,:,:)),[],3)),std(max(squeeze(freqDeconv(:,iRoi,ibin,rwd,:,:)),[],3))./sqrt(numSubs),plotStyles{1},'Color',plotColors{rwd},'linewidth',errbarwidth,'markersize',markersize);%not plotting for "task"
        hold on
    end
%     title([roiNames{iRoi} ' bin ' num2str(ibin) ' max'])
    if ibin==1; legend({'high rwd','low rwd'},'location','northwest'); ylabel('max'); end
    xticks(1:numFreqs);
    xticklabels(num2str(stimulus.freqs','%0.1f'));
    xlabel('frequency (cpd)');
%     ylim([ymin ymax]);
end

%% Fit SF tuning curves
options = optimset('MaxFunEvals',1000, 'MaxIter',1000);%used by fminsearch
initParams = [3, 1, 1, 0];
for ibin=1:nbins
    for rwd=1:2
        respStd = squeeze(mean(freqDeconvStd(:,ibin,rwd,:)));%responses
        myfunc = @(x) gaussianRms(x, respStd, (1:numFreqs)');
        [x fval exflag] = fminsearch(myfunc, initParams, options);
        x(3) = abs(x(3));
        sfStdTuning(ibin,rwd,:) = x;
        sfStdExitFlag(ibin,rwd) = exflag;
        sfStdRms(ibin,rwd) = fval;%fit RMS 
        
        respBetas = squeeze(meanSfBetas(:,:,ibin,rwd,:));%responses
        myfunc = @(x) gaussianRms(x, respBetas, (1:numFreqs)');
        [x fval exflag] = fminsearch(myfunc, initParams, options);
        x(3) = abs(x(3));
        sfBetasTuning(ibin,rwd,:) = x;
        sfBetasExitFlag(ibin,rwd) = exflag;
        sfBetasRms(ibin,rwd) = fval;%fit RMS 
        
        respFft = squeeze(meanFreqFftAmp(:,ibin,rwd,:));%responses
        myfunc = @(x) gaussianRms(x, respFft, (1:numFreqs)');
        [x fval exflag] = fminsearch(myfunc, initParams, options);
        x(3) = abs(x(3));
        sfFftTuning(ibin,rwd,:) = x;
        sfFftExitFlag(ibin,rwd) = exflag;
        sfFftRms(ibin,rwd) = fval;%fit RMS 
        
        respMax = mean(max(squeeze(freqDeconv(:,iRoi,ibin,rwd,:,:)),[],3))';%responses
        myfunc = @(x) gaussianRms(x, respMax, (1:numFreqs)');
        [x fval exflag] = fminsearch(myfunc, initParams, options);
        x(3) = abs(x(3));
        sfMaxTuning(ibin,rwd,:) = x;
        sfMaxExitFlag(ibin,rwd) = exflag;
        sfMaxRms(ibin,rwd) = fval;%fit RMS 
    end
end


%% SF
% i=i+1;
% figure(i) 
% clf
% rows=4;

interpFreqs = 0:0.1:numFreqs+1;
row=3;
for ibin=1:nbins
    subplot(rows,cols,(row-1)*cols+ibin)
    for rwd=1:2
        x= sfStdTuning(ibin,rwd,:);
        z = x(4) + x(3)*exp(-((interpFreqs-x(1)).^2/2/x(2)^2));
        plot(interpFreqs,z,plotStyles{2},'Color',plotColors{rwd},'linewidth',linewidth);
        hold on
    end
%     if ibin==1; legend({'high rwd','low rwd'},'location','northwest'); ylabel('STD'); end
end
row=row+1;
for ibin=1:nbins
    subplot(rows,cols,(row-1)*cols+ibin)
    for rwd=1:2
        x= sfBetasTuning(ibin,rwd,:);
        z = x(4) + x(3)*exp(-((interpFreqs-x(1)).^2/2/x(2)^2));
        plot(interpFreqs,z,plotStyles{2},'Color',plotColors{rwd},'linewidth',linewidth);
        hold on
    end
%     if ibin==1; ylabel('betas'); end
end
row=row+1;
for ibin=1:nbins
    subplot(rows,cols,(row-1)*cols+ibin)
    for rwd=1:2
        x= sfFftTuning(ibin,rwd,:);
        z = x(4) + x(3)*exp(-((interpFreqs-x(1)).^2/2/x(2)^2));
        plot(interpFreqs,z,plotStyles{2},'Color',plotColors{rwd},'linewidth',linewidth);
        hold on
    end
%     if ibin==1; ylabel('fft amplitude'); end
end
row=row+1;
for ibin=1:nbins
    subplot(rows,cols,(row-1)*cols+ibin)
    for rwd=1:2
        x= sfMaxTuning(ibin,rwd,:);
        z = x(4) + x(3)*exp(-((interpFreqs-x(1)).^2/2/x(2)^2));
        plot(interpFreqs,z,plotStyles{2},'Color',plotColors{rwd},'linewidth',linewidth);
        hold on
    end
%     if ibin==1; ylabel('max'); end
end





%%
%Spatial frequency STD per subject
i=i+1;
figure(i) 
clf
rows=nbins;
cols=numSubs;
row=0;
%average fMRI response for SF
row=row+1;
for iSub=1:numSubs
    for ibin=1:nbins
        subplot(rows,cols,iSub+(ibin-1)*cols)
        for rwd=1:2
            plot(squeeze(freqDeconvStd(iSub,ibin,rwd,:)),plotStyles{2},'Color',plotColors{rwd},'linewidth',linewidth);%not plotting for "task"
            hold on
        end
%         title([roiNames{iRoi} ' bin ' num2str(ibin) ' std'])
        if ibin==1 & iSub==1; title('std'); legend({'high rwd','low rwd'},'location','northwest'); end
        if ibin>1
            xticks([]);           
        end
    end
end

%%
%Spatial frequency BETAS per subject
i=i+1;
figure(i) 
clf
rows=nbins;
cols=numSubs;
row=0;
%average fMRI response for SF
row=row+1;
for iSub=1:numSubs
    for ibin=1:nbins
        subplot(rows,cols,iSub+(ibin-1)*cols)
        for rwd=1:2
            plot(squeeze(sfBetas(iSub,iRoi,ibin,rwd,:)),plotStyles{2},'Color',plotColors{rwd},'linewidth',linewidth);%not plotting for "task"
            hold on
        end
%         title([roiNames{iRoi} ' bin ' num2str(ibin) ' std'])
        if ibin==1 & iSub==1; title('betas'); legend({'high rwd','low rwd'},'location','northwest'); end
        if ibin>1
            xticks([]);           
        end
    end
end

%%
%Spatial frequency MAX per subject
i=i+1;
figure(i) 
clf
rows=nbins;
cols=numSubs;
row=0;
%average fMRI response for SF
row=row+1;
for iSub=1:numSubs
    for ibin=1:nbins
        subplot(rows,cols,iSub+(ibin-1)*cols)
        for rwd=1:2
            plot(max(squeeze(freqDeconv(iSub,iRoi,ibin,rwd,:,:))'),plotStyles{2},'Color',plotColors{rwd},'linewidth',linewidth);%not plotting for "task"
            hold on
        end
%         title([roiNames{iRoi} ' bin ' num2str(ibin) ' std'])
        if ibin==1 & iSub==1; title('max'); legend({'high rwd','low rwd'},'location','northwest'); end
        if ibin>1
            xticks([]);           
        end
    end
end









%% SF + contrast




i=i+1;
figure(i) 
clf
rows=5;
cols=nbins;
row=0;

%deconvolved fMRI response for SF+contrast
row=row+1;
meanContrastFreqDeconv = mean(contrastFreqDeconv,1);%average over subjects
for ibin=1:nbins
    subplot(rows,cols,ibin+(row-1)*cols)
    for rwd=1:2
        for c=1:size(meanContrastFreqDeconv,5)
            plot(squeeze(meanContrastFreqDeconv(:,:,ibin,rwd,c,:))',plotStyles{2},'Color',plotColors{rwd},'linewidth',linewidth);

            hold on
        end
    end
    xlabel('time (TR)');
    title([num2str(binBorders(ibin))  '\circ - ' num2str(binBorders(ibin+1)) '\circ']);
%     if ibin==1; legend({'low contrast','high contrast'},'location','northwest'); end
end

%std of deconvolved fMRI response for SF+contrast
ymin = min(contrastFreqDeconvStd(:));
ymax = max(contrastFreqDeconvStd(:));
row=row+1;
for ibin=1:nbins
    subplot(rows,cols,ibin+(row-1)*cols)
    for rwd=1:2
%         plot(squeeze(mean(contrastFreqDeconvStd(:,ibin,rwd,1:end-1)))',plotStyles{1},'Color',plotColors{rwd},'linewidth',linewidth);%not plotting for "task"
        errorbar(1:numContrasts*numFreqs, squeeze(mean(contrastFreqDeconvStd(:,ibin,rwd,:)))',squeeze(std(contrastFreqDeconvStd(:,ibin,rwd,:)))'./sqrt(numSubs),plotStyles{1},'Color',plotColors{rwd},'linewidth',errbarwidth,'markersize',markersize);%not plotting for "task"
        hold on
    end
    ylim([ymin ymax]);
    if ibin==1; legend({'high rwd','low rwd'},'location','northwest'); ylabel('std'); end
end
%std of deconvolved fMRI response for SF+contrast
row=row+1;
for ibin=1:nbins
    subplot(rows,cols,ibin+(row-1)*cols)
    for rwd=1:2
%         plot(squeeze(mean(contrastFreqDeconvStd(:,ibin,rwd,1:end-1)))',plotStyles{1},'Color',plotColors{rwd},'linewidth',linewidth);%not plotting for "task"
        errorbar(1:numContrasts*numFreqs, squeeze(mean(contrastFreqDeconvStd(:,ibin,rwd,:)))',squeeze(std(contrastFreqDeconvStd(:,ibin,rwd,:)))'./sqrt(numSubs),plotStyles{1},'Color',plotColors{rwd},'linewidth',errbarwidth,'markersize',markersize);%not plotting for "task"
        hold on
    end
    if ibin==1; legend({'high rwd','low rwd'},'location','northwest'); ylabel('std'); end
end

% 
%betas for SF+contrast
row=row+1;
meanSfContrastBetas = mean(sfContrastBetas,1);%nullBetas(iSub,iRoi,ibin,rwd,:)
stdSfContrastBetas = std(sfContrastBetas);%nullBetas(iSub,iRoi,ibin,rwd,:)
for ibin=1:nbins
    subplot(rows,cols,ibin+(row-1)*cols)
    for rwd=1:2
%         plot(squeeze(meanSfContrastBetas(:,:,ibin,rwd,1:end-1))',plotStyles{1},'Color',plotColors{rwd},'linewidth',linewidth);%not plotting for "task"
        errorbar(1:numContrasts*numFreqs, squeeze(meanSfContrastBetas(:,:,ibin,rwd,:))',squeeze(stdSfContrastBetas(:,:,ibin,rwd,:))'./sqrt(numSubs),plotStyles{1},'Color',plotColors{rwd},'linewidth',errbarwidth,'markersize',markersize);%not plotting for "task"

        hold on
    end
%     title([roiNames{iRoi} ' bin ' num2str(ibin) ' betas'])
    if ibin==1; legend({'high rwd','low rwd'},'location','northwest'); ylabel('betas'); end
end

%average MAX SF+contrast
row=row+1;
for ibin=1:nbins
    subplot(rows,cols,ibin+(row-1)*cols)
    for rwd=1:2
        plot(mean(max(squeeze(contrastFreqDeconv(:,iRoi,ibin,rwd,1:end-1,:)),[],3)),plotStyles{1},'Color',plotColors{rwd},'linewidth',linewidth);%not plotting for "task"
        errorbar(1:numContrasts*numFreqs, mean(max(squeeze(contrastFreqDeconv(:,iRoi,ibin,rwd,:,:)),[],3)),std(max(squeeze(contrastFreqDeconv(:,iRoi,ibin,rwd,:,:)),[],3))./sqrt(numSubs),plotStyles{1},'Color',plotColors{rwd},'linewidth',errbarwidth,'markersize',markersize);%not plotting for "task"

        hold on
    end
    if ibin==1; legend({'high rwd','low rwd'},'location','northwest'); ylabel('max'); end
    xticks(1:2:numContrasts*numFreqs);
    temp = num2str(stimulus.freqs','%0.1f');
    dbltemp = [ temp; temp];
    xticklabels(dbltemp(1:2:end,:));
    xlabel('frequency (cpd)');
end



%% Fit SF+contrast tuning curves
options = optimset('MaxFunEvals',1000, 'MaxIter',1000);%used by fminsearch
initParams = [3, 1, 1, 0];
for ibin=1:nbins
    for rwd=1:2
        for c=1:numContrasts
            respStd = squeeze(mean(contrastFreqDeconvStd(:,ibin,rwd,(c-1)*numFreqs+1:c*numFreqs)));%responses
            myfunc = @(x) gaussianRms(x, respStd, (1:numFreqs)');
            [x fval exflag] = fminsearch(myfunc, initParams, options);
            x(3) = abs(x(3));
            sfContrastStdTuning(ibin,rwd,c,:) = x;
            sfContrastStdExitFlag(ibin,rwd,c) = exflag;
            sfContrastStdRms(ibin,rwd,c) = fval;%fit RMS
            
            respBetas = squeeze(meanSfContrastBetas(:,:,ibin,rwd,(c-1)*numFreqs+1:c*numFreqs));%responses
            myfunc = @(x) gaussianRms(x, respBetas, (1:numFreqs)');
            [x fval exflag] = fminsearch(myfunc, initParams, options);
            x(3) = abs(x(3));
            sfContrastBetasTuning(ibin,rwd,c,:) = x;
            sfContrastBetasExitFlag(ibin,rwd) = exflag;
            sfContrastBetasRms(ibin,rwd) = fval;%fit RMS
            
            respMax = mean(max(squeeze(contrastFreqDeconv(:,iRoi,ibin,rwd,(c-1)*numFreqs+1:c*numFreqs,:)),[],3))';%responses
            myfunc = @(x) gaussianRms(x, respMax, (1:numFreqs)');
            [x fval exflag] = fminsearch(myfunc, initParams, options);
            x(3) = abs(x(3));
            sfContrastMaxTuning(ibin,rwd,c,:) = x;
            sfContrastMaxExitFlag(ibin,rwd) = exflag;
            sfContrastMaxRms(ibin,rwd) = fval;%fit RMS
        end
    end
end



%%
interpFreqs = 0:0.1:numFreqs+1;
row=3;
for ibin=1:nbins
    subplot(rows,cols,(row-1)*cols+ibin)
    for rwd=1:2
        for c=1:numContrasts
            x= sfContrastStdTuning(ibin,rwd,c,:);
            z = x(4) + x(3)*exp(-((interpFreqs-x(1)).^2/2/x(2)^2));
            plot((c-1)*(interpFreqs(end)-1)+interpFreqs,z,plotStyles{2},'Color',plotColors{rwd},'linewidth',linewidth);
            hold on
        end
    end
%     if ibin==1; legend({'high rwd','low rwd'},'location','northwest'); ylabel('STD'); end
end

row=row+1;
for ibin=1:nbins
    subplot(rows,cols,(row-1)*cols+ibin)
    for rwd=1:2
        for c=1:numContrasts
            x= sfContrastBetasTuning(ibin,rwd,c,:);
            z = x(4) + x(3)*exp(-((interpFreqs-x(1)).^2/2/x(2)^2));
            plot((c-1)*(interpFreqs(end)-1)+interpFreqs,z,plotStyles{2},'Color',plotColors{rwd},'linewidth',linewidth);
            hold on
        end
    end
%     if ibin==1; ylabel('betas'); end
end
row=row+1;
for ibin=1:nbins
    subplot(rows,cols,(row-1)*cols+ibin)
    for rwd=1:2
        for c=1:numContrasts
        x= sfContrastMaxTuning(ibin,rwd,c,:);
        z = x(4) + x(3)*exp(-((interpFreqs-x(1)).^2/2/x(2)^2));
        plot((c-1)*(interpFreqs(end)-1)+interpFreqs,z,plotStyles{2},'Color',plotColors{rwd},'linewidth',linewidth);
        hold on
        end
    end
%     if ibin==1; ylabel('fft amplitude'); end
end


%% graph of high-low rwd effects, as function of eccentricity
i=i+1;
figure(i) 
clf
rows=2;
cols=4;
j=0;


%average task response STD
j=j+1;
subplot(rows,cols,j)
meanTrialStd = squeeze(mean(trialTCstd));
plot(meanTrialStd,'linewidth',linewidth);
legend(rwdString);
title('mean response STD');
%H minus L average task response STD
subplot(rows,cols,j+cols)
plot(meanTrialStd(:,1)-meanTrialStd(:,2),'k','linewidth',linewidth);
title('H - L');
%average stim-null response STD
j=j+1;
subplot(rows,cols,j)
meanNullStd = squeeze(mean(nullDeconvStd));
diffNullStd = squeeze(meanNullStd(:,:,1) - meanNullStd(:,:,2));
plot(diffNullStd,'linewidth',linewidth);%stim-null
title('stim - null STD');
legend(rwdString);
%H minus L stim-null STD
subplot(rows,cols,j+cols)
plot(diffNullStd(:,1) - diffNullStd(:,2),'k','linewidth',linewidth);%stim-null
title('H - L');

%average high/low contrast response STD
j=j+1;
subplot(rows,cols,j)
meanContrastStd = squeeze(mean(contrastDeconvStd));
diffContrastStd = squeeze(meanContrastStd(:,:,1) - meanContrastStd(:,:,2));

plot(diffContrastStd,'linewidth',linewidth);%high contrast - low contrast
title('high - low contrast STD');
legend(rwdString);
%H-L high/low contrast response STD
subplot(rows,cols,j+cols)
plot(diffContrastStd(:,1) - diffContrastStd(:,2),'k','linewidth',linewidth);%high contrast - low contrast
title('H - L');
%%
squeeze(sum(numBinVoxels))



%%
%%%%%
function [c] = gaussianRms(x, orig, xi)
% mean = x(1);
% sigma = x(2);
% amp = x(3);
% baseline = x(4);
x(3) = abs(x(3));%so that we don't have inverted tuning curves
zi = x(4) + x(3)*exp(-((xi-x(1)).^2/2/x(2)^2));

c = sqrt(sum((orig-zi).^2));

end
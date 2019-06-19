dataFolder = '~/data/rwdFmri/';
load([dataFolder 'rwdRapidData_combinedGLM.mat'], 'dataFolder', 'subFolders', 'samplerate', 'roiNames', ...
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
    'freqDeconvFftAmp','freqDeconvFftPh','contrastFreqDeconvFftAmp','contrastFreqDeconvFftPh');
rwdString = {'H','L'};
plotColors = {[0 0 1], [1 0 0], [0 1 0], [0.5 1 0.2], [1 0.2 0.5] };
plotStyles = {'-','--',':','-.','-','--',':','-.'};

%%

for iSub = 1:numSubs
    subString{iSub} = subFolders{iSub}(1:4);
end

i=0;
rows=8;
cols=nbins;
linewidth = 2;


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
        
            plot(squeeze(meanBinNullTC(:,:,ibin,rwd,c,:))',plotStyles{c},'Color',plotColors{rwd},'linewidth',linewidth);
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
        plot(squeeze(mean(nullDeconvStd(:,ibin,rwd,:)))',plotStyles{1},'Color',plotColors{rwd},'linewidth',linewidth);
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
        plot(squeeze(meanNullBetas(:,:,ibin,rwd,:))',plotStyles{1},'Color',plotColors{rwd},'linewidth',linewidth);
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
            plot(squeeze(meanBinContrastTC(:,:,ibin,rwd,c,:))',plotStyles{c},'Color',plotColors{rwd},'linewidth',linewidth);
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
        plot(squeeze(mean(contrastDeconvStd(:,ibin,rwd,:)))',plotStyles{1},'Color',plotColors{rwd},'linewidth',linewidth);
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
        plot(squeeze(meanContrastBetas(:,:,ibin,rwd,:))',plotStyles{1},'Color',plotColors{rwd},'linewidth',linewidth);
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
            plot(squeeze(meanFreqBinTC(:,:,ibin,rwd,c,:))',plotStyles{c},'Color',plotColors{rwd},'linewidth',linewidth);
            hold on
        end
    end
    xlabel('time (TR)');
    title([roiNames{iRoi} ' bin ' num2str(ibin)])
%     if ibin==1; legend({'low contrast','high contrast'},'location','northwest'); end
end

%std of deconvolved fMRI response for SF
row=row+1;
for ibin=1:nbins
    subplot(rows,cols,ibin+(row-1)*cols)
    for rwd=1:2
        plot(squeeze(mean(freqDeconvStd(:,ibin,rwd,1:end-1)))',plotStyles{1},'Color',plotColors{rwd},'linewidth',linewidth);%not plotting for "task"
        hold on
    end
    title([roiNames{iRoi} ' bin ' num2str(ibin) ' std'])
    if ibin==1; legend({'high rwd','low rwd'},'location','northwest'); end
end
% 
%betas for SF
row=row+1;
meanSfBetas = mean(sfBetas,1);%nullBetas(iSub,iRoi,ibin,rwd,:)
for ibin=1:nbins
    subplot(rows,cols,ibin+(row-1)*cols)
    for rwd=1:2
        plot(squeeze(meanSfBetas(:,:,ibin,rwd,1:end-1))',plotStyles{1},'Color',plotColors{rwd},'linewidth',linewidth);%not plotting for "task"
        hold on
    end
    title([roiNames{iRoi} ' bin ' num2str(ibin) ' betas'])
    if ibin==1; legend({'high rwd','low rwd'},'location','northwest'); end
end

%FFT amp for SF
row=row+1;
meanFreqFftAmp = mean(freqDeconvFftAmp,1);%nullBetas(iSub,iRoi,ibin,rwd,:)
for ibin=1:nbins
    subplot(rows,cols,ibin+(row-1)*cols)
    for rwd=1:2
        plot(squeeze(meanFreqFftAmp(:,ibin,rwd,1:end-1))',plotStyles{1},'Color',plotColors{rwd},'linewidth',linewidth);%not plotting for "task"
        hold on
    end
    title([roiNames{iRoi} ' bin ' num2str(ibin) ' FFT'])
    if ibin==1; legend({'high rwd','low rwd'},'location','northwest'); end
end

%average MAX SF
row=row+1;
for ibin=1:nbins
    subplot(rows,cols,ibin+(row-1)*cols)
    for rwd=1:2
        plot(mean(max(squeeze(freqDeconv(:,iRoi,ibin,rwd,1:end-1,:)),[],3)),plotStyles{1},'Color',plotColors{rwd},'linewidth',linewidth);%not plotting for "task"
%         plot(mean(squeeze(freqDeconv(:,iRoi,ibin,rwd,1:end-1,5))),plotStyles{1},'Color',plotColors{rwd},'linewidth',linewidth);%not plotting for "task"

        hold on
    end
    title([roiNames{iRoi} ' bin ' num2str(ibin) ' max'])
    if ibin==1; legend({'high rwd','low rwd'},'location','northwest'); end
end

%MAX of average SF timecourse
row=row+1;
for ibin=1:nbins
    subplot(rows,cols,ibin+(row-1)*cols)
    for rwd=1:2
        plot(max(squeeze(meanFreqBinTC(:,iRoi,ibin,rwd,1:end-1,:)),[],2),plotStyles{1},'Color',plotColors{rwd},'linewidth',linewidth);%not plotting for "task"
%         plot(mean(squeeze(freqDeconv(:,iRoi,ibin,rwd,1:end-1,5))),plotStyles{1},'Color',plotColors{rwd},'linewidth',linewidth);%not plotting for "task"

        hold on
    end
    title([roiNames{iRoi} ' bin ' num2str(ibin) ' max'])
    if ibin==1; legend({'high rwd','low rwd'},'location','northwest'); end
end



%%
%% SF
i=i+1;
figure(i) 
clf
rows=5;

row=0;

%deconvolved fMRI response for SF
row=row+1;
meanContrastFreqDeconv = mean(contrastFreqDeconv,1);%average over subjects

for ibin=1:nbins
    subplot(rows,cols,ibin+(row-1)*cols)
    for rwd=1:2
        for c=1:size(meanContrastFreqDeconv,5)
            plot(squeeze(meanContrastFreqDeconv(:,:,ibin,rwd,c,:))',plotStyles{1},'Color',plotColors{rwd},'linewidth',linewidth);
            hold on
        end
    end
    xlabel('time (TR)');
    title([roiNames{iRoi} ' bin ' num2str(ibin)])
%     if ibin==1; legend({'low contrast','high contrast'},'location','northwest'); end
end

%std of deconvolved fMRI response for SF
row=row+1;
for ibin=1:nbins
    subplot(rows,cols,ibin+(row-1)*cols)
    for rwd=1:2
        plot(squeeze(mean(contrastFreqDeconvStd(:,ibin,rwd,1:end-1)))',plotStyles{1},'Color',plotColors{rwd},'linewidth',linewidth);%not plotting for "task"
        hold on
    end
    title([roiNames{iRoi} ' bin ' num2str(ibin) ' std'])
    if ibin==1; legend({'high rwd','low rwd'},'location','northwest'); end
end
% 
%betas for SF
row=row+1;
meanSfContrastBetas = mean(sfContrastBetas,1);%nullBetas(iSub,iRoi,ibin,rwd,:)
for ibin=1:nbins
    subplot(rows,cols,ibin+(row-1)*cols)
    for rwd=1:2
        plot(squeeze(meanSfContrastBetas(:,:,ibin,rwd,1:end-1))',plotStyles{1},'Color',plotColors{rwd},'linewidth',linewidth);%not plotting for "task"
        hold on
    end
    title([roiNames{iRoi} ' bin ' num2str(ibin) ' betas'])
    if ibin==1; legend({'high rwd','low rwd'},'location','northwest'); end
end

%average MAX SF
row=row+1;
for ibin=1:nbins
    subplot(rows,cols,ibin+(row-1)*cols)
    for rwd=1:2
        plot(mean(max(squeeze(contrastFreqDeconv(:,iRoi,ibin,rwd,1:end-1,:)),[],3)),plotStyles{1},'Color',plotColors{rwd},'linewidth',linewidth);%not plotting for "task"
%         plot(mean(squeeze(freqDeconv(:,iRoi,ibin,rwd,1:end-1,5))),plotStyles{1},'Color',plotColors{rwd},'linewidth',linewidth);%not plotting for "task"

        hold on
    end
    title([roiNames{iRoi} ' bin ' num2str(ibin) ' max'])
    if ibin==1; legend({'high rwd','low rwd'},'location','northwest'); end
end

%MAX of average SF+contrast timecourse
meanContrastFreqBinTC = mean(contrastFreqDeconv,1);%average over subjects
row=row+1;
for ibin=1:nbins
    subplot(rows,cols,ibin+(row-1)*cols)
    for rwd=1:2
        plot(max(squeeze(meanContrastFreqBinTC(:,iRoi,ibin,rwd,1:end-1,:)),[],2),plotStyles{1},'Color',plotColors{rwd},'linewidth',linewidth);%not plotting for "task"
%         plot(mean(squeeze(freqDeconv(:,iRoi,ibin,rwd,1:end-1,5))),plotStyles{1},'Color',plotColors{rwd},'linewidth',linewidth);%not plotting for "task"

        hold on
    end
    title([roiNames{iRoi} ' bin ' num2str(ibin) ' max'])
    if ibin==1; legend({'high rwd','low rwd'},'location','northwest'); end
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
            plot(squeeze(freqDeconvStd(iSub,ibin,rwd,1:end-1)),plotStyles{1},'Color',plotColors{rwd},'linewidth',linewidth);%not plotting for "task"
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
            plot(squeeze(sfBetas(iSub,iRoi,ibin,rwd,1:end-1)),plotStyles{1},'Color',plotColors{rwd},'linewidth',linewidth);%not plotting for "task"
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
            plot(max(squeeze(freqDeconv(iSub,iRoi,ibin,rwd,1:end-1,:))'),plotStyles{1},'Color',plotColors{rwd},'linewidth',linewidth);%not plotting for "task"
            hold on
        end
%         title([roiNames{iRoi} ' bin ' num2str(ibin) ' std'])
        if ibin==1 & iSub==1; title('max'); legend({'high rwd','low rwd'},'location','northwest'); end
        if ibin>1
            xticks([]);           
        end
    end
end

%%
%Spatial frequency timecourses per subject
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
            plot(squeeze(freqDeconv(iSub,iRoi,ibin,rwd,1:end-1,:))',plotStyles{1},'Color',plotColors{rwd},'linewidth',linewidth);%not plotting for "task"
            hold on
            plot(squeeze(freqDeconv(iSub,iRoi,ibin,rwd,end,:))',plotStyles{2},'Color',plotColors{rwd},'linewidth',linewidth);%not plotting for "task"
        end
%         title([roiNames{iRoi} ' bin ' num2str(ibin) ' std'])
        if ibin==1 & iSub==1; legend({'high rwd','low rwd'},'location','northwest'); end
        if ibin>1
            xticks([]);           
        end
    end
end

%%
%Spatial frequency response 3d plot
i=i+1;
figure(i) 
clf
rows=2;
cols=nbins;
row=0;
%average fMRI response for SF
row=row+1;
%freqDeconv(iSub,iRoi,ibin,rwd,:,:)
% meanFreqBinTC = mean(meanSubFreqBinTC,4);%average over rwd
X=repmat(1:deconvLength, numFreqs,1);
Y=repmat(1:numFreqs, deconvLength,1)';

for ibin=1:nbins
    for rwd=1:2
        subplot(rows,cols,ibin+(rwd-1)*cols)
        Z=squeeze(meanFreqBinTC(:,:,ibin,rwd,1:numFreqs,:));
        surf(X,Y,Z);
        title([roiNames{iRoi} ' ' rwdString{rwd} ' rwd']);
        xlabel('time (TR)');
        ylabel('frequency');
        title([roiNames{iRoi} ' bin ' num2str(ibin) ' ' rwdString{rwd}]);
    end
end


%%
%Spatial frequency + contrast response 3d plot 
i=i+1;
figure(i) 
clf
rows=2;
cols=nbins;
row=0;
%average fMRI response for SF
row=row+1;
%freqDeconv(iSub,iRoi,ibin,rwd,:,:)
% meanFreqBinTC = mean(meanSubFreqBinTC,4);%average over rwd
X=repmat(1:deconvLength, numContrasts * numFreqs,1);
Y=repmat(1:(numContrasts * numFreqs), deconvLength,1)';

for ibin=1:nbins
    for rwd=1:2
        subplot(rows,cols,ibin+(rwd-1)*cols)
        Z=squeeze(meanContrastFreqDeconv(:,:,ibin,rwd,1:(numContrasts * numFreqs),:));
        surf(X,Y,Z);
%         title([roiNames{iRoi} ' ' rwdString{rwd} ' rwd']);
%         xlabel('time (TR)');
%         ylabel('frequency');
%         title([roiNames{iRoi} ' bin ' num2str(ibin) ' ' rwdString{rwd}]);
    end
end
% 
% %%
% %Spatial frequency + contrast response 3d plot - high-low contrast
% i=i+1;
% figure(i) 
% clf
% rows=2;
% cols=nbins;
% row=0;
% %average fMRI response for SF
% row=row+1;
% %freqDeconv(iSub,iRoi,ibin,rwd,:,:)
% % meanFreqBinTC = mean(meanSubFreqBinTC,4);%average over rwd
% X=repmat(1:deconvLength, numFreqs,1);
% Y=repmat(1:numFreqs, deconvLength,1)';
% 
% for ibin=1:nbins
%     for rwd=1:2
%         subplot(rows,cols,ibin+(rwd-1)*cols)
%         tempLow = meanContrastFreqDeconv(:,:,ibin,rwd,1:numFreqs,:);
%         tempHigh = meanContrastFreqDeconv(:,:,ibin,rwd,numFreqs+1:(numContrasts*numFreqs),:);
%         Z=squeeze(tempHigh-tempLow);
%         surf(X,Y,Z);
% %         title([roiNames{iRoi} ' ' rwdString{rwd} ' rwd']);
%         xlabel('time (TR)');
%         ylabel('frequency');
%         title([roiNames{iRoi} ' bin ' num2str(ibin) ' ' rwdString{rwd}]);
%     end
% end

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
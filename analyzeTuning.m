%uses files saved by savePhysio.m, saveGroupEyeData.m, saveRoiTC.m

close all
mrQuit
clear all

iRoi=1;

onlyCorrect=0;%0=all trials, 1=only correct, -1=only incorrect, 2=valid response
subtractMean = 1;
toZscore=1;%0 or 1
concatProj= 1;
% fillMissing = 1;
% bensonROIs = 1;%[1:3];
eccMin = 0.2;
eccMax = 70;
nbins = 12;
binBorders = logspace(log10(eccMin),log10(eccMax),nbins+1);
nbins = length(binBorders)-1;
for i=2:length(binBorders)
    binCenters(i-1) = (binBorders(i)+binBorders(i-1))/2;
end
binBorders;

nfreqs=5;
ncontrasts=2;
tic

dataFolder = '/Users/rothzn/rwdFmri/';
% dataFolder = 'c:\rwdFmri\';
freqs = logspace(-0.3,0.5,5);
contrasts = logspace(-0.7,0,2);


onlyCorrectStr='';
switch onlyCorrect
    case 0
        onlyCorrectStr='_all';
    case 1
        onlyCorrectStr='_correct';
    case -1
        onlyCorrectStr='_incorrect';
    case 2
        onlyCorrectStr='_valid';
end
zscoreStr='';
if toZscore
    zscoreStr = '_zscore';
end
subtractMeanStr = '';
if subtractMean
    subtractMeanStr = '_subtractMean';
end
concatProjStr = '';
if concatProj
    concatProjStr = 'proj';
end
load([dataFolder 'allData' concatProjStr onlyCorrectStr zscoreStr subtractMeanStr '.mat'],...
        'pupilTrialTime','goodSubs', 'subFolders',...
        'onlyCorrect','ncontrasts','nfreqs',...
        'subtractMean',...
        'toZscore','concatProj','bensonROIs','roiNames','eccMin','eccMax','nbins',...
        'binBorders','nbins',...
        'meanContrastPupil','meanContrastPulse','meanFreqPupil','meanFreqPulse',...
        'subBaseStimPupil','stimPupilConcat','subStimPulse','subBaseStimPupil',...
        'subStdStimPupil','subStdStimPulse','subStimRwd','stimTrials',...
        'stimTrials','subStimPupilBaseMedian','subStimPupilStdMedian', 'subStimPulseStdMedian', 'subStimRwdMedian',...
        'subBetasStimPupil','subBetasStimPulse',...
        'stimFreqTrials','stimContrastTrials',...
        'subMeanStimPulse',...
        'subBinBetasStimVox', 'binMeanStimTrial','binMeanStimTrialVox','subBinBetasStimTrial',...
        'subBinBetasFreq','subBinBetasContrast','rwdTrials',...
        'pupilMeanNull', 'pulseMeanNull','voxMeanNull','binMeanNullTrialVox','binMeanNullTrial');
% plotColors = {[0, 0.4470, 0.7410], [0.8500, 0.3250, 0.0980],[0.9290, 0.6940, 0.1250],...
%     [0.4940, 0.1840, 0.5560], [0.4660, 0.6740, 0.1880], [0.3010, 0.7450, 0.9330], [0.6350, 0.0780, 0.1840]};
plotColors = {[1 0 0], [0 0 1]};
plotStyles = {'.','-','--',':','-.','-','--',':','-.'};
linewidth = 2;
markersize=15;
errbarwidth = 1;
%%


ifig=0;
ifig=ifig+1; figure(ifig); clf
rows=3;
cols=nbins;
for i=1:3
    switch i
        case 1
            meanTrial = squeeze(mean(binMeanStimTrial(:,iRoi,:,:),1)); %mean over subjects
            ylabelStr = 'stim';
        case 2
            meanTrial = squeeze(mean(binMeanNullTrial(:,iRoi,:,:),1)); %mean over subjects
            ylabelStr = 'null';
        case 3
            meanTrial = squeeze(mean(binMeanStimTrial(:,iRoi,:,:)-binMeanNullTrial(:,iRoi,:,:),1)); %mean over subjects
            ylabelStr = 'stim-null';
    end
    ymin(i) = min(min(meanTrial(:)),0);
    ymax(i) = max(0.1+meanTrial(:));
    for ibin=1:nbins
        subplot(rows,cols,ibin+(i-1)*cols);
        plot(meanTrial(ibin,:),'k');
        axis square        
        if i==1
            title([num2str(binBorders(ibin),'%.1f') char(176) '-' num2str(binBorders(ibin+1),'%.1f') char(176)]);
        end
        xlim([0 1+size(meanTrial,2)]);
        if i<rows
            xticks([]);
        else
           xlabel('time (TR)');
           xticks([1 size(meanTrial,2)]);
        end
        if ibin==1
            ylabel(ylabelStr);
        end
    end
end
for isubplot=1:rows*cols
    subplot(rows,cols,isubplot)
    if mod(isubplot,cols)==1
%         ylabel(ylabelStr);
%         maxBeta = ceil(10*maxBeta)/10;
        yticks([floor(10*min(ymin))/10 ceil(10*max(ymax))/10]);
%         yticklabels([num2str(floor(100*min(ymin)/100) ceil(10*max(ymax)/100)]);
    else
        yticks([]);
    end
    ylim([floor(10*min(ymin))/10 ceil(10*max(ymax))/10]);
    
%     ylim([1.1*min(ymin) 1.1*max(ymax)]);
end


set(gcf,'position',[100 300 1100 300]);



ifig=ifig+1; figure(ifig); clf
rows=2;
cols=nbins;
temp = mean(subBinBetasContrast);
maxBeta = max(temp(:));
minBeta = min(0,min(temp(:)));
maxBeta = ceil(10*maxBeta)/10;
for ibin=1:nbins
    subplot(rows,cols,ibin)
    meanTuning = squeeze(mean(subBinBetasContrast(:,iRoi,ibin,:,:)));%rwd,icontrast
    stdTuning = squeeze(std(subBinBetasContrast(:,iRoi,ibin,:,:)))./sqrt(length(goodSubs));%rwd,ifreq
    for rwd=1:2
        plot(meanTuning(rwd,:)','.-','color',plotColors{rwd},'markersize',markersize);
        hold on
        %         errorbar(meanTuning(rwd,:)',stdTuning(rwd,:)','.-','color',plotColors{rwd},'linewidth',errbarwidth,'markersize',markersize);
    end
    
    xticklabels({'low','high'})
    xticks(1:ncontrasts);
    ylim([minBeta maxBeta]);
    if ibin==1
        ylabel('betas');
    else
        yticks([]);
    end
    title([num2str(binBorders(ibin),'%.1f') char(176) '-' num2str(binBorders(ibin+1),'%.1f') char(176)]);
    xlabel('contrast');
    xlim([0.5 ncontrasts+0.5]);
end
% legend('high reward','low reward');



temp = mean(subBinBetasFreq);
maxBeta = max(temp(:));
minBeta = min(0,min(temp(:)));
maxBeta = ceil(10*maxBeta)/10;
for ibin=1:nbins
    subplot(rows,cols,cols+ibin)
    meanTuning = squeeze(mean(subBinBetasFreq(:,iRoi,ibin,:,:)));%rwd,ifreq
    stdTuning = squeeze(std(subBinBetasFreq(:,iRoi,ibin,:,:)))./sqrt(length(goodSubs));%rwd,ifreq
    for rwd=1:2
        plot(meanTuning(rwd,:)','.-','color',plotColors{rwd},'markersize',markersize);
        hold on
        %         errorbar(meanTuning(rwd,:)',stdTuning(rwd,:)','.-','color',plotColors{rwd},'linewidth',errbarwidth,'markersize',markersize);
    end
    yticks([minBeta, maxBeta]);
    ylim([minBeta maxBeta]);
    if ibin==1
        ylabel('betas');
    else
        yticks([]);
    end
    xlabel('frequency');
    xlim([0 nfreqs+1]);
    xticklabels({num2str(freqs(1),'%.1f'), num2str(freqs(nfreqs),'%.1f')});
end
for isubplot=1:rows*cols
    subplot(rows,cols,isubplot)
    axis square
end
set(gcf,'position',[50 150 1100 300])
savepdf(gcf, ['~/Documents/MATLAB/rwdRapid/figures/tuning_' concatProjStr onlyCorrectStr zscoreStr subtractMeanStr '.pdf']);

%% Single Subject
% ifig=2;
% % ifig=ifig+1;
% figure(ifig); clf
% iSub = 4;
% temp = subBinBetasContrast(iSub,iRoi,:,:,:);
% maxBeta = max(temp(:));
% minBeta = min(0,min(temp(:)));
% maxBeta = ceil(10*maxBeta)/10;
% for ibin=1:nbins
%     subplot(rows,cols,ibin)
%     temp = squeeze(subBinBetasContrast(iSub,iRoi,ibin,:,:));%rwd,icontrast
%     for rwd=1:2
%         plot(temp(rwd,:)','.-','color',plotColors{rwd});
%         hold on
%     end
%     xticklabels({'low','high'})
%     xticks(1:ncontrasts);
%     ylim([minBeta maxBeta]);
%     if ibin==1
%         ylabel('betas');
%     else
%         yticks([]);
%     end
%     title([num2str(binBorders(ibin),'%.1f') char(176) '-' num2str(binBorders(ibin+1),'%.1f') char(176)]);
%     xlabel('contrast');
%     xlim([0.5 ncontrasts+0.5]);
% end
% legend('high reward','low reward');
%
%
%
% temp = subBinBetasFreq(iSub,iRoi,:,:,:);
% maxBeta = max(temp(:));
% minBeta = min(0,min(temp(:)));
% maxBeta = ceil(10*maxBeta)/10;
% for ibin=1:nbins
%     subplot(rows,cols,cols+ibin)
%     temp = squeeze(subBinBetasFreq(iSub,iRoi,ibin,:,:));%rwd,ifreq
%     for rwd=1:2
%         plot(temp(rwd,:)','.-','color',plotColors{rwd});
%         hold on
%     end
%     yticks([minBeta, maxBeta]);
%     ylim([minBeta maxBeta]);
%     if ibin==1
%        ylabel('betas');
%     else
%         yticks([]);
%     end
%     xlabel('frequency');
%     xlim([0 nfreqs+1]);
%     xticklabels({num2str(freqs(1),'%.1f'), num2str(freqs(nfreqs),'%.1f')});
% end
% for isubplot=1:rows*cols
%    subplot(rows,cols,isubplot)
%    axis square
% end
% set(gcf,'position',[50 150 1100 300])

%% 3D surface plot
% ifig=ifig+1;
% figure(ifig); clf
% rows=1;
% cols=3;
% for rwd=1:2
%     subplot(rows,cols,rwd)
%     temp = squeeze(mean(subBinBetasFreq(:,iRoi,:,rwd,:)));
%     surf(temp);
% end
% subplot(rows,cols,3)
% temp = squeeze(mean(subBinBetasFreq(:,iRoi,:,1,:)-subBinBetasFreq(:,iRoi,:,2,:)));
% surf(temp);
% set(gcf,'position',[50 150 1100 400])

%% fit tuning functions for each subject
options = optimset('MaxFunEvals',1000, 'MaxIter',1000);%used by fminsearch
initParams = [3, 1, 1, 0];

for iroi=bensonROIs
    for ibin=1:nbins
        for rwd=1:2
            for iSub=1:length(goodSubs)
                freqBetas = squeeze(subBinBetasFreq(iSub,iroi,ibin,rwd,:));
                myfunc = @(x) gaussianRms(x, freqBetas, (1:nfreqs)');
                [x, fval, exflag] = fminsearch(myfunc, initParams, options);
                x(3) = abs(x(3));
                subFreqTuning(iSub,iroi,ibin,rwd,:) = x;
                subFreqExitFlag(iSub,iroi,ibin,rwd) = exflag;
                subFreqRms(iSub,iroi,ibin,rwd) = fval;%fit RMS
            end
            %fit tuning curves to population mean responses
            freqBetas = squeeze(mean(subBinBetasFreq(:,iroi,ibin,rwd,:)));
            myfunc = @(x) gaussianRms(x, freqBetas, (1:nfreqs)');
            [x, fval, exflag] = fminsearch(myfunc, initParams, options);
            x(3) = abs(x(3));
            x(1) = abs(x(1));
            freqTuning(iroi,ibin,rwd,:) = x;
            freqExitFlag(iroi,ibin,rwd) = exflag;
            freqRms(iroi,ibin,rwd) = fval;%fit RMS
        end
    end
end
interpFreqs = 0:0.1:nfreqs+1;
for ibin=1:nbins
    subplot(rows,cols,cols+ibin)
    for rwd=1:2
        x = squeeze(freqTuning(iRoi,ibin,rwd,:));
        z = x(4) + x(3)*exp(-((interpFreqs-x(1)).^2/2/x(2)^2));
        plot(interpFreqs,z,plotStyles{2},'Color',plotColors{rwd},'linewidth',linewidth);
    end
end

%% plot fit parameters as function of eccentricity
ifig=ifig+1; figure(ifig); clf
rows=2;
cols=length(x);%number of parameters
paramNames = {'preferred freq','tuning width','amplitude','baseline'};
plotBins = 2:8;
for iparam=1:length(x)
   subplot(rows,cols,iparam)
   for rwd=1:2
      plot( squeeze(freqTuning(iRoi,plotBins,rwd,iparam)),'color',plotColors{rwd});
      hold on; axis square
      title(paramNames{iparam});
   end
   subplot(rows,cols,iparam+cols)
   for rwd=1:2
       plot(squeeze(median(subFreqTuning(:,iRoi,plotBins,rwd,iparam))),'color',plotColors{rwd}); 
       hold on; axis square
   end
end

%% bar plots of single subject fit parameters, for high and low reward, per eccentricity bin
for ibin=1:nbins
    for iparam=1:length(x)
        [h pvalParam(ibin,iparam)] = ttest(subFreqTuning(:,iRoi,ibin,1,iparam) - subFreqTuning(:,iRoi,ibin,2,iparam));
    end
end
%%
%%%%%
function [c] = gaussianRms(x, orig, xi)
% mean = x(1);
% sigma = x(2);
% amp = x(3);
% baseline = x(4);
x(3) = abs(x(3));%so that we don't have inverted tuning curves
x(1) = abs(x(1));%so that preferred frequecy is not negative
zi = x(4) + x(3)*exp(-((xi-x(1)).^2/2/x(2)^2));

c = sqrt(sum((orig-zi).^2));

end
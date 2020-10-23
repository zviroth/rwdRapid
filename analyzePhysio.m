%uses files saved by savePhysio.m, saveGroupEyeData.m, saveRoiTC.m

close all
mrQuit
clear all

iRoi=1;

onlyCorrect=2;%0=all trials, 1=only correct, -1=only incorrect, 2=valid response
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
        'pupilMeanNull', 'pulseMeanNull','voxMeanNull','binMeanNullTrialVox','binMeanNullTrial',...
        'pulseMeanNullRwd','pulseMeanStimRwd','pulseMeanRwd','pupilMeanNullRwd','pupilMeanStimRwd','pupilMeanRwd');
                               
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
rows=2;
cols=7;

j=0;

for j=1:6
    switch j
        case 1
            data = pupilMeanNullRwd([1:6 8:length(goodSubs)],:,:);
            titleStr = 'pupil null';
        case 2
            data = pupilMeanStimRwd([1:6 8:length(goodSubs)],:,:);
            titleStr = 'pupil stim';
        case 3
            data = pupilMeanRwd([1:6 8:length(goodSubs)],:,:);
            titleStr = 'pupil';
        case 4
            data = pulseMeanNullRwd;
            titleStr = 'pulse null';
        case 5
            data = pulseMeanStimRwd;
            titleStr = 'pulse stim';
        case 6
            data = pulseMeanRwd;
            titleStr = 'pulse';
    end
    subplot(rows,cols,j)
    for rwd=1:2
        plot(squeeze(nanmean(data(:,rwd,:))),'color',plotColors{rwd}); hold on;
        title(titleStr);
    end
    axis square
end

goodSubs = [1 3:10 13:14];
goodSubs = [1:14];
stdPupil = squeeze(nanstd(pupilMeanRwd(goodSubs,:,:),0,3));
meanPupil = squeeze(nanmean(pupilMeanRwd(goodSubs,:,:),3));
basePupil = squeeze(nanmean(pupilMeanRwd(goodSubs,:,1:50),3));
stdPulse = squeeze(std(pulseMeanRwd(goodSubs,:,:),0,3));
meanPulse = squeeze(mean(pulseMeanRwd(goodSubs,:,:),3));
stdPupil = stdPupil./stdPupil(:,1);
meanPupil = meanPupil./meanPupil(:,1);
basePupil = basePupil./basePupil(:,1);
stdPulse = stdPulse./stdPulse(:,1);
meanPulse = meanPulse./meanPulse(:,1);
for isub=1:length(goodSubs)
    for rwd=1:2
        meanStimTrialStdPupil(isub,rwd) = nanmean(subStdStimPupil{goodSubs(isub)}(rwdTrials{goodSubs(isub)}==rwd));
        meanStimTrialStdPulse(isub,rwd) = nanmean(subStdStimPulse{goodSubs(isub)}(rwdTrials{goodSubs(isub)}==rwd));
    end
end

for j=1:7
    switch j
        case 1
            data = stdPupil;
            titleStr = 'pupil STD';
        case 2
            data = meanPupil;
            titleStr = 'pupil Mean';
        case 3
            data = basePupil;
            titleStr = 'pupil Baseline';
        case 4
            data = stdPulse;
            titleStr = 'pulse STD';
        case 5
            data = meanPulse;
            titleStr = 'pulse Mean';
        case 6
            data = meanStimTrialStdPupil;
            titleStr = 'pupil trial STD';
        case 7
            data = meanStimTrialStdPulse;
            titleStr = 'pulse trial STD';
    end
    subplot(rows,cols,cols+j)
    plot(data');
    hold on
    plot(nanmean(data),'k.','markersize',markersize);
    hline(1);
    title(titleStr);
    axis square
    [h physioPval(j)] = ttest(data(:,1) - data(:,2));
end
set(gcf,'position',[100 200 1100 300]);
physioPval




%%
for isub=1:length(goodSubs)
   temp = squeeze(pupilMeanNullRwd(isub,:,:)); 
    mx(isub) = max(abs(temp(:)));
    s(isub) = std(temp(:));
end
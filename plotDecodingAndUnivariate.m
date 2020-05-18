%uses file created by analyzeDecoding_arousal_perm2.m
close all
clear all
tic
onlyCorrect=1;%0=all trials, 1=only correct, -1=only incorrect, 2=valid response
includeControl=0;%left and right DMN trial amplitude
includePulse = 0;
subtractMean = 0;
includePupil = 0;
includePupilBase = 0;
toZscore=1;%0 or 1
concatProj= 1;

sigHeight5 = 0.08;
sigHeight1 = 0.11;
sigHeightT = 0.05;
pvalThresh=0.05;


nfreqs=5;
ncontrasts=2;

% dataFolder = 'c:\rwdFmri\';
dataFolder = '/Users/rothzn/rwdFmri/';
ifig=0;
%%
onlyCorrectStr='';
switch onlyCorrect
%     case 0
%         onlyCorrectStr='_correct';
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
controlStr='';
if includeControl
    controlStr='_controlROI';
end
pulseStr='';
if includePulse
    pulseStr = '_pulse';
end
pupilStr='';
if includePupil
    pupilStr = '_pupilStd';
end
pupilBaseStr='';
if includePupilBase
    pupilBaseStr = '_pupilBasePupil';
end
subtractMeanStr = '';
if subtractMean
    subtractMeanStr = '_subtractMean';
end

load([dataFolder 'decodeArousalPerm' onlyCorrectStr zscoreStr controlStr pulseStr pupilStr pupilBaseStr subtractMeanStr '.mat'],...
    'goodSubs', 'subFolders',...
    'onlyCorrect', 'featureTypes','ncontrasts','nfreqs','contrastCondsTrue','freqCondsTrue',...
    'nperms', 'arousalTypes','includeControl','includePulse','subtractMean',...
    'includePupil','includePupilBase','featureTypes','toZscore','concatProj','bensonROIs','eccMin','eccMax','nbins',...
    'binBorders','nbins',...
    'freqCondsGuess','contrastCondsGuess',...
    'permGuessFreq','permGuessContrast','arousalTrials',...
    'classAccContrast','classAccFreq',...
    'freqCondsTrue','contrastCondsTrue',...    
    'permAccFreqDiff','permAccContrastDiff',...
    'permArousalTrials','realDiffFreq','realDiffContrast',...
    'permFreqMeanDiff','permContrastMeanDiff',...
    'pvalPermFreq','pvalPermContrast','pvalFreq','pvalContrast',...
    'meanContrastPupil','meanContrastPulse','meanFreqPupil','meanFreqPulse',...
    'freqCondsTrue','contrastCondsTrue',...
    'subBaseStimPupil','stimPupilConcat','subStimPulse','subBaseStimPupil',...
    'subStdStimPupil','subStdStimPulse','subStimRwd','controlStimBetas','stimTrials',...
    'stimTrials','subStimPupilBaseMedian','subStimPupilStdMedian', 'subStimPulseStdMedian', 'subStimRwdMedian',...
    'controlStimBetasMedian',...
    'subBinBetasStim','subBetasStimPupil','subBetasStimPulse',...
    'stimFreqTrials','stimContrastTrials',...
    'subMeanStimPulse');
clear permArousalTrials
npermsDecoding=nperms;
load([dataFolder 'univariateArousalPerm' onlyCorrectStr zscoreStr controlStr pulseStr pupilStr pupilBaseStr subtractMeanStr '.mat'],...
    'npermsDecoding','npermsUnivariate',...
    'arousalTypes','meanPupil','plotColors','pvalPermBetasPupil','meanPulse',...
    'pvalPermBetasPulse','subBinContrastBetas','subBinFreqBetas',...
    'pvalDiffBinContrastBetas','diffBinContrastBetas','meanContrastVar',...
    'meanFreqVar','subBinContrastBetasStd','subBinFreqBetasStd','subCnr',...
    'pvalPermBasePupil','pvalPermStdPupil','pvalPermStdPulse','pvalPermMeanPulse',...
    'pvalPermBetasPupil','pvalDiffBinContrastBetasStd');
nperms = npermsUnivariate;
%%
for arousalType=arousalTypes
    for iSub=1:length(goodSubs)
        for arousal=1:2
             rwdInd(1) = 0;
             rwdInd(2) = sum(stimTrials{iSub,1});
             rwdInd(3) = sum(stimTrials{iSub,1}) + sum(stimTrials{iSub,2});
             for rwd=1:2
                 %how many good trials in this rwd are in this arousal
                 rwdConds(arousalType,iSub,arousal,rwd) = sum(arousalTrials{arousalType,iSub,arousal}([1+rwdInd(rwd):rwdInd(rwd+1)]));
             end
        end
    end
end
plotColors = { [1 0 0], [0 0 1], [0 1 0], [0.5 1 0.2], [1 0.2 0.5] };
lineStyles = {'-','--'};
markerSize=60;
rows=1;
cols=7;
iRoi=1;
for arousalType=arousalTypes
    ifig=ifig+1; figure(ifig); clf

    subplot(rows,cols,1:2)
    
    for featureType=3
        for arousal=1:2
            plot(squeeze(nanmean(classAccFreq{arousalType,featureType,arousal})),'color',plotColors{arousal})
            hold on
        end
        for ibin=1:nbins
            for arousal=1:2
                maxAcc(arousal) = max(squeeze(nanmean(classAccFreq{arousalType,featureType,arousal}(:,iRoi,ibin))));
            end
            if pvalFreq(arousalType,featureType,iRoi,ibin)<pvalThresh
                
                if maxAcc(1)>maxAcc(2)
                    scatter(ibin,max(maxAcc)+sigHeightT,markerSize,plotColors{1});
                else
                    scatter(ibin,max(maxAcc)+sigHeightT,markerSize,plotColors{2});
                end
            end
            if pvalPermFreq(arousalType,featureType,iRoi,ibin)<pvalThresh
                scatter(ibin,max(maxAcc)+sigHeight5,markerSize,plotColors{1},'*','filled');
                if pvalPermFreq(arousalType,featureType,iRoi,ibin)<0.01
                    scatter(ibin,max(maxAcc)+sigHeight1,markerSize,plotColors{1},'*','filled');
                end
            elseif pvalPermFreq(arousalType,featureType,iRoi,ibin)>1-pvalThresh
                scatter(ibin,max(maxAcc)+sigHeight5,markerSize,plotColors{2},'*');
                if pvalPermFreq(arousalType,featureType,iRoi,ibin)>1-0.01
                    scatter(ibin,max(maxAcc)+sigHeight1,markerSize,plotColors{2},'*');
                end
            end
        end
        %     mean(squeeze(mean(classAccFreq{featureType})))
    end
    % legend('std','amp','regress');
    if arousalType==4
        legend('low','high');
    else
        legend('high','low');
    end

    
    for featureType=3
        for arousal=1:2
            plot(squeeze(nanmean(classAccContrast{arousalType,featureType,arousal})),'color',plotColors{arousal})
        end
        for ibin=1:nbins
            for arousal=1:2
                maxAcc(arousal) = max(squeeze(mean(classAccContrast{arousalType,featureType,arousal}(:,iRoi,ibin))));
            end
            if pvalContrast(arousalType,featureType,iRoi,ibin)<pvalThresh
                if maxAcc(1)>maxAcc(2)
                    scatter(ibin,max(maxAcc)+sigHeightT,markerSize,plotColors{1});
                else
                    scatter(ibin,max(maxAcc)+sigHeightT,markerSize,plotColors{2});
                end
            end
            if pvalPermContrast(arousalType,featureType,iRoi,ibin)<pvalThresh
                scatter(ibin,max(maxAcc)+sigHeight5,markerSize,plotColors{1},'*');
                if pvalPermContrast(arousalType,featureType,iRoi,ibin)<0.01
                    scatter(ibin,max(maxAcc)+sigHeight1,markerSize,plotColors{1},'*');
                end
            elseif pvalPermContrast(arousalType,featureType,iRoi,ibin)>1-pvalThresh
                scatter(ibin,max(maxAcc)+sigHeight5,markerSize,plotColors{2},'*');
                if pvalPermContrast(arousalType,featureType,iRoi,ibin)>1-0.01
                    scatter(ibin,max(maxAcc)+sigHeight1,markerSize,plotColors{2},'*');
                end
            end
        end
    end
    
    
    
    xlabel('eccentricity bin');
    ylabel('proportion correct');
    xlim([0 nbins+1]);
    % legend('frequency','contrast');
    hline(1/nfreqs);
    hline(1/ncontrasts);
    titleStr = [strrep(onlyCorrectStr,'_','') ', '];

    switch arousalType
        case 1
            arousalStr = 'pupil baseline';
        case 2
            arousalStr = 'pupil std';
        case 3
            arousalStr = 'pulse std';
        case 4
            arousalStr = 'reward';
        case 5
            arousalStr = 'DMN';
    end
    titleStr = [titleStr arousalStr ', '];
%     titleStr = onlyCorrectStr;
    if includePupil
        titleStr = [titleStr 'pupil, '];
    end
    if includePulse
        titleStr = [titleStr 'pulse, '];
    end
    if includePupilBase
        titleStr = [titleStr 'pupil baseline, '];
    end
    if subtractMean
        titleStr = [titleStr 'mean subtracted '];
    end

    title(titleStr);

    old = {' ', ','};
    new = '_';

    
    
    %number of trials for each frequency condition,by arousal
    subplot(rows,cols,3)
    for arousal=1:2
        bar((arousal-1)*nfreqs+[1:5],squeeze(mean(freqCondsTrue(arousalType,featureType,:,arousal,:))),'facecolor',plotColors{arousal},'linestyle',lineStyles{arousal});
        hold on
        er=errorbar((arousal-1)*nfreqs+[1:nfreqs],squeeze(mean(freqCondsTrue(arousalType,featureType,:,arousal,:))),squeeze(std(freqCondsTrue(arousalType,featureType,:,arousal,:))));
                er.Color = [0 0 0];
        er.LineStyle = 'none';
%         histogram(squeeze(mean(freqCondsTrue(:,arousal,:))),nfreqs+1,'facecolor',plotColors{arousal},'linestyle',lineStyles{arousal});
    end
    xlim([0 2*nfreqs+1]);
    ylabel('mean #trials');
    title('spatial frequency')
    
    %number of trials for each contrast condition,by arousal
    subplot(rows,cols,4)
    for arousal=1:2
        bar((arousal-1)*ncontrasts+[1:ncontrasts],squeeze(mean(contrastCondsTrue(arousalType,featureType,:,arousal,:))),'facecolor',plotColors{arousal},'linestyle',lineStyles{arousal});
        hold on
        er=errorbar((arousal-1)*ncontrasts+[1:ncontrasts],squeeze(mean(contrastCondsTrue(arousalType,featureType,:,arousal,:))),squeeze(std(contrastCondsTrue(arousalType,featureType,:,arousal,:))));
        er.Color = [0 0 0];
        er.LineStyle = 'none';
%         histogram(squeeze(mean(contrastCondsTrue(:,arousal,:))),ncontrasts+1,'facecolor',plotColors{arousal},'linestyle',lineStyles{arousal});
    end
    xlim([0 2*ncontrasts+1]);
    ylabel('mean #trials');
    title('contrast')
    
    
    %number of trials for each rwd condition,by arousal
    subplot(rows,cols,5)
%     rwdConds(arousalType,iSub,arousal,rwd);
    for arousal=1:2
        bar((arousal-1)*2+[1:2],squeeze(mean(rwdConds(arousalType,:,arousal,:))),'facecolor',plotColors{arousal},'linestyle',lineStyles{arousal});
        
        %         bar((arousal-1)*2+[1:2],squeeze(mean(contrastCondsTrue(arousalType,featureType,:,arousal,:))),'facecolor',plotColors{arousal},'linestyle',lineStyles{arousal});
        hold on
        er=errorbar((arousal-1)*2+[1:2],squeeze(mean(rwdConds(arousalType,:,arousal,:))),squeeze(std(rwdConds(arousalType,:,arousal,:))));
        er.Color = [0 0 0];
        er.LineStyle = 'none';
    end
    ylabel('mean #trials');
    title('reward type')
        
    %plot mean pupil timecourse for high and low arousal
    subplot(rows,cols,6)
    for arousal=1:2
       plot(squeeze(nanmean(meanPupil(arousalType,:,arousal,:))),'color',plotColors{arousal}); 
       hold on
    end
    title(['pupil (p=' num2str(pvalPermBetasPupil(arousalType),'%.4f') ')'])
 
    %plot mean pulse timecourse for high and low arousal
    subplot(rows,cols,7)
    for arousal=1:2
        plot(squeeze(nanmean(meanPulse(arousalType,:,arousal,:))),'color',plotColors{arousal});
        hold on
    end
%     title('pulse')
    title(['pulse (p=' num2str(pvalPermBetasPulse(arousalType),'%.4f') ')'])
    
    set(gcf,'position',[100 100 1100 400])
end

%%

%%
rows=2;
cols=5;
markerSize=40;
sigHeight1=0.08;
sigHeight5=0.05;
for arousalType=arousalTypes
    ifig=ifig+1; figure(ifig); clf

    
    %plot contrast betas as function of eccentricity for both arousal levels
    subplot(rows,cols,1)
    for arousal=1:2
        for icontrast=1:ncontrasts
            plot(squeeze(nanmean(subBinContrastBetas(arousalType,:,:,arousal,icontrast))),...
                'color',plotColors{arousal},'linestyle',lineStyles{icontrast});
            hold on
        end
    end
    xlim([0 nbins+1]);
    switch arousalType
        case 1
            arousalStr = 'pupil baseline';
        case 2
            arousalStr = 'pupil std';
        case 3
            arousalStr = 'pulse std';
        case 4
            arousalStr = 'reward';
        case 5
            arousalStr = 'DMN';
    end
    title(['contrast: ' arousalStr]);
    
    %plot difference between high and low contrast betas as function of eccentricity for both arousal levels
    subplot(rows,cols,2)

    arousalBetasDiff = squeeze(nanmean(subBinContrastBetas(arousalType,:,:,:,1)-subBinContrastBetas(arousalType,:,:,:,2)));
       for arousal=1:2
           plot(arousalBetasDiff(:,arousal),'color',plotColors{arousal});
        hold on
       end
    for ibin=1:nbins
        if pvalDiffBinContrastBetas(arousalType,ibin)<pvalThresh
            scatter(ibin,arousalBetasDiff(ibin,1) + sigHeight5,markerSize,plotColors{1},'*');
            if pvalDiffBinContrastBetas(arousalType,ibin)<0.01
                scatter(ibin,arousalBetasDiff(ibin,1) + sigHeight1,markerSize,plotColors{1},'*');
            end
        elseif pvalDiffBinContrastBetas(arousalType,ibin)>1-pvalThresh
            scatter(ibin,arousalBetasDiff(ibin,2) + sigHeight5,markerSize,plotColors{2},'*');
            if pvalDiffBinContrastBetas(arousalType,ibin)>1-0.01
                scatter(ibin,arousalBetasDiff(ibin,2) + sigHeight1,markerSize,plotColors{2},'*');
            end
        end
    end

    xlim([0 nbins+1]);
    title('high-low contrast');
    if arousalType==4
        legend('low','high');
    else
        legend('high','low');
    end
    
    %plot difference between betas for high and low arousal as function of
    %eccentricity for both contrasts
    subplot(rows,cols,3)
    for icontrast=1:ncontrasts
        plot(squeeze(nanmean(diffBinContrastBetas(arousalType,:,:,icontrast))),'color',plotColors{icontrast});
        hold on
    end
    xlim([0 nbins+1]);

    legend('high contrast','low contrast');
    if arousalType==4
        title('low-high arousal');
    else
        title('high-low arousal');
    end
    
    %plot contrast variability
    subplot(rows,cols,4)
   for arousal=1:2
        plot(squeeze(meanContrastVar(arousalType,:,arousal)),'color',plotColors{arousal});
        hold on
    end
    xlim([0 nbins+1]);
    title('mean variability')
    if arousalType==4
        legend('low','high');
    else
        legend('high','low');
    end
    for ibin=1:nbins
        if pvalDiffBinContrastBetasStd(arousalType,ibin)<pvalThresh
            scatter(ibin,squeeze(meanContrastVar(arousalType,ibin,1)) + sigHeight5,markerSize,plotColors{1},'*');
            if pvalDiffBinContrastBetasStd(arousalType,ibin)<0.01
                scatter(ibin,squeeze(meanContrastVar(arousalType,ibin,1)) + sigHeight1,markerSize,plotColors{1},'*');
            end
        elseif pvalDiffBinContrastBetasStd(arousalType,ibin)>1-pvalThresh
            scatter(ibin,squeeze(meanContrastVar(arousalType,ibin,2)) + sigHeight5,markerSize,plotColors{2},'*');
            if pvalDiffBinContrastBetasStd(arousalType,ibin)>1-0.01
                scatter(ibin,squeeze(meanContrastVar(arousalType,ibin,2)) + sigHeight1,markerSize,plotColors{2},'*');
            end
        end
    end
    
    %plot contrast CNR
    subplot(rows,cols,5)
    %first CNR is the mean betas divided by the mean variance
    subBinContrastMeanBetas(arousalType,:,:,:) = squeeze(mean(abs(subBinContrastBetas(arousalType,:,:,:,:)),5));%average over conditions. sub,bin,arousal
    subContrastVar(arousalType,:,:,:) = squeeze(nanmean(subBinContrastBetasStd(arousalType,:,:,:,:),5));%average over conditions
    subCnr(arousalType,:,:,:) = subBinContrastMeanBetas(arousalType,:,:,:)./subContrastVar(arousalType,:,:,:);%sub,bin,arousal
    
    %second CNR is the difference between the betas, divided by the mean variance
%     subCnr2(arousalType,:,:,:) = abs(diffBinContrastBetas(arousalType,:,:,:))./subContrastVar(arousalType,:,:,:);%sub,bin,arousal

    for arousal=1:2
        plot(squeeze(mean(subCnr(arousalType,:,:,arousal),2)),'color',plotColors{arousal},'linestyle',lineStyles{1});
        hold on
   end
    if arousalType==4
        plot(squeeze(mean(subCnr(arousalType,:,:,2) - subCnr(arousalType,:,:,1),2)), 'k','linestyle',lineStyles{1});
    else
        plot(squeeze(mean(subCnr(arousalType,:,:,1) - subCnr(arousalType,:,:,2),2)), 'k','linestyle',lineStyles{1});
    end
    title('CNR');
    
    %plot frequency betas as function of eccentricity for both arousal levels
    subplot(rows,cols,6)
    for arousal=1:2
        for ifreq=1:nfreqs
            plot(squeeze(nanmean(subBinFreqBetas(arousalType,:,:,arousal,ifreq))),...
                'color',plotColors{arousal}+(1-plotColors{arousal})/(1/nfreqs+ifreq));
            hold on
        end
    end
    xlim([0 nbins+1]);
    title('frequency');
   
    
       %plot frequency variability
    subplot(rows,cols,7)
    for arousal=1:2
        plot(squeeze(meanFreqVar(arousalType,:,arousal)),'color',plotColors{arousal});
        hold on
    end
    xlim([0 nbins+1]);
    title('mean variability')
    if arousalType==4
        legend('low','high');
    else
        legend('high','low');
    end
    
        %plot frequency CNR
    subplot(rows,cols,8)
    %first CNR is the mean betas divided by the mean variance
    subBinFreqMeanBetas(arousalType,:,:,:) = squeeze(mean(abs(subBinFreqBetas(arousalType,:,:,:,:)),5));%average over conditions. sub,bin,arousal
    subContrastVar(arousalType,:,:,:) = squeeze(nanmean(subBinFreqBetasStd(arousalType,:,:,:,:),5));%average over conditions
    subCnr(arousalType,:,:,:) = subBinContrastMeanBetas(arousalType,:,:,:)./subContrastVar(arousalType,:,:,:);%sub,bin,arousal

    for arousal=1:2
        plot(squeeze(mean(subCnr(arousalType,:,:,arousal),2)),'color',plotColors{arousal},'linestyle',lineStyles{1});
        hold on
   end
    if arousalType==4
        plot(squeeze(mean(subCnr(arousalType,:,:,2) - subCnr(arousalType,:,:,1),2)), 'k','linestyle',lineStyles{1});
    else
        plot(squeeze(mean(subCnr(arousalType,:,:,1) - subCnr(arousalType,:,:,2),2)), 'k','linestyle',lineStyles{1});
    end
    title('CNR');
    
    set(gcf,'position',[70 80 1100 400])
end
toc
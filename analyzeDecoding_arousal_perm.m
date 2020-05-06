close all

clear all
pvalThresh=0.01;
onlyCorrect=2;%0=all trials, 1=only correct, -1=only incorrect, 2=valid response
includeControl=0;%left and right DMN trial amplitude
includePulse = 0;
subtractMean = 0;
includePupil = 0;
includePupilBase = 0;
toZscore=1;%0 or 1
concatProj= 1;

nfreqs=5;
ncontrasts=2;
tic

dataFolder = 'c:\rwdFmri\';
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

load([dataFolder 'decodeArousalPerm' zscoreStr controlStr pulseStr pupilStr pupilBaseStr subtractMeanStr '.mat'],...
    'featureTypes','ncontrasts','nfreqs','contrastCondsTrue','freqCondsTrue',...
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
    'freqCondsTrue','contrastCondsTrue');

%%
ifig=0;
iRoi=1;
for arousalType=arousalTypes
    plotColors = { [1 0 0], [0 0 1], [0 1 0], [0.5 1 0.2], [1 0.2 0.5] };
    lineStyles = {'-','--'};
    ifig=ifig+1; figure(ifig); clf
    rows=1;
    cols=4;
    subplot(rows,cols,1:2)

    markerSize=80;
    for featureType=featureTypes
        for arousal=1:2
            plot(squeeze(nanmean(classAccFreq{arousalType,featureType,arousal})),lineStyles{arousal},'color',plotColors{featureType})
            hold on
        end
        for ibin=1:nbins
            for arousal=1:2
                maxAcc(arousal) = max(squeeze(nanmean(classAccFreq{arousalType,featureType,arousal}(:,iRoi,ibin))));
                %             maxAcc = max(squeeze(mean(accFreq{featureType,1}(:,iRoi,ibin))),squeeze(mean(accFreq{featureType,2}(:,iRoi,ibin))));
            end
            if pvalFreq(arousalType,featureType,iRoi,ibin)<pvalThresh

                if maxAcc(1)>maxAcc(2)
                    scatter(ibin,max(maxAcc),markerSize,plotColors{featureType},'filled');
                else
                    scatter(ibin,max(maxAcc),markerSize,plotColors{featureType});
                end
            end
            if pvalPermFreq(arousalType,featureType,iRoi,ibin)<pvalThresh
                scatter(ibin,max(maxAcc),markerSize,plotColors{featureType},'s','filled');
            elseif pvalPermFreq(arousalType,featureType,iRoi,ibin)>1-pvalThresh
                scatter(ibin,max(maxAcc),markerSize,plotColors{featureType},'s');
            end
        end
        %     mean(squeeze(mean(classAccFreq{featureType})))
    end
    % legend('std','amp','regress');
    legend('high','low');
    for featureType=featureTypes
        for arousal=1:2
            plot(squeeze(nanmean(classAccContrast{arousalType,featureType,arousal})),lineStyles{arousal},'color',plotColors{featureType})
        end
        for ibin=1:nbins
            for arousal=1:2
                maxAcc(arousal) = max(squeeze(mean(classAccContrast{arousalType,featureType,arousal}(:,iRoi,ibin))));
                %             maxAcc = max(squeeze(mean(accFreq{featureType,1}(:,iRoi,ibin))),squeeze(mean(accFreq{featureType,2}(:,iRoi,ibin))));
            end
            if pvalContrast(arousalType,featureType,iRoi,ibin)<pvalThresh
                
                if maxAcc(1)>maxAcc(2)
                    scatter(ibin,max(maxAcc),markerSize,plotColors{featureType},'filled');
                else
                    scatter(ibin,max(maxAcc),markerSize,plotColors{featureType});
                end
            end
        end
        %     mean(squeeze(mean(classAccContrast{featureType})))
    end
    xlabel('eccentricity bin');
    ylabel('proportion correct');
    % legend('frequency','contrast');
    hline(1/nfreqs);
    hline(1/ncontrasts);
    titleStr = 'decode ';
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
    
%     subplot(rows,cols,2)
%     for featureType=featureTypes
%         for arousal=1:2
%             plot(squeeze(max(classAccFreq{featureType,arousal},[],3)),lineStyles{arousal},'color',plotColors{featureType})
%             hold on
%             plot(squeeze(max(classAccContrast{featureType,arousal},[],3)),lineStyles{arousal},'color',plotColors{featureType})
%         end
%     end
%     xlabel('subject #');
%     hline(1/nfreqs);
%     hline(1/ncontrasts);
    
    old = {' ', ','};
    new = '_';
    % savefig([dataFolder replace(titleStr,old, new)]);
    for arousal=1:2
        [max(squeeze(mean(classAccContrast{arousalType,featureType,arousal},1))), max(squeeze(mean(classAccFreq{arousalType,featureType,arousal},1)))]
    end
    
    
    
    subplot(rows,cols,3)
    for arousal=1:2
        bar((arousal-1)*nfreqs+[1:5],squeeze(mean(freqCondsTrue(arousalType,featureType,:,arousal,:))),'facecolor',plotColors{arousal},'linestyle',lineStyles{arousal});
        hold on
        er=errorbar((arousal-1)*nfreqs+[1:nfreqs],squeeze(mean(freqCondsTrue(arousalType,featureType,:,arousal,:))),squeeze(std(freqCondsTrue(arousalType,featureType,:,arousal,:))));
                er.Color = [0 0 0];
        er.LineStyle = 'none';
%         histogram(squeeze(mean(freqCondsTrue(:,arousal,:))),nfreqs+1,'facecolor',plotColors{arousal},'linestyle',lineStyles{arousal});
        
    end
    subplot(rows,cols,4)
    for arousal=1:2
        bar((arousal-1)*ncontrasts+[1:ncontrasts],squeeze(mean(contrastCondsTrue(arousalType,featureType,:,arousal,:))),'facecolor',plotColors{arousal},'linestyle',lineStyles{arousal});
        hold on
        er=errorbar((arousal-1)*ncontrasts+[1:ncontrasts],squeeze(mean(contrastCondsTrue(arousalType,featureType,:,arousal,:))),squeeze(std(contrastCondsTrue(arousalType,featureType,:,arousal,:))));
        er.Color = [0 0 0];
        er.LineStyle = 'none';
%         histogram(squeeze(mean(contrastCondsTrue(:,arousal,:))),ncontrasts+1,'facecolor',plotColors{arousal},'linestyle',lineStyles{arousal});
        
    end
    set(gcf,'position',[100 100 1100 400])
end



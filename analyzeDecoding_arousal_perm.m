close all

clear all
sigHeight5 = 0.08;
sigHeight1 = 0.11;
sigHeightT = 0.05;
pvalThresh=0.05;
onlyCorrect=1;%0=all trials, 1=only correct, -1=only incorrect, 2=valid response
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

% dataFolder = 'c:\rwdFmri\';
dataFolder = '/Users/rothzn/rwdFmri/';
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
    'controlStimBetasMedian');
clear permArousalTrials
%%
ifig=0;
iRoi=1;
for arousalType=arousalTypes
    %permutation test for difference in pupil baseline, pupil std, pulse
    %std, 
    for iSub=1:length(goodSubs)
        for arousal=1:2
            meanBasePupil(arousalType,iSub,arousal) = nanmean(subBaseStimPupil{iSub}(arousalTrials{arousalType,iSub,arousal}));
             meanStdPupil(arousalType,iSub,arousal) = nanmean(subStdStimPupil{iSub}(arousalTrials{arousalType,iSub,arousal}));
             meanStdPulse(arousalType,iSub,arousal) = nanmean(subStdStimPulse{iSub}(arousalTrials{arousalType,iSub,arousal}));
             meanPupil(arousalType,iSub,arousal,:) = nanmean(stimPupilConcat{iSub}(arousalTrials{arousalType,iSub,arousal},:));
             meanPulse(arousalType,iSub,arousal,:) = nanmean(subStimPulse{iSub}(:,arousalTrials{arousalType,iSub,arousal}),2);
             
             rwdInd(1) = 0;
             rwdInd(2) = sum(stimTrials{iSub,1});
             rwdInd(3) = sum(stimTrials{iSub,1}) + sum(stimTrials{iSub,2});
             for rwd=1:2
                 %how many good trials in this rwd are in this arousal
                 rwdConds(arousalType,iSub,arousal,rwd) = sum(arousalTrials{arousalType,iSub,arousal}([1+rwdInd(rwd):rwdInd(rwd+1)]));
             end
%                  rwdConds(arousalType,iSub,arousal,rwd) = sum(arousalTrials{arousalType,iSub,arousal});
                 %number of total trials in rwd:
                 %size(stimTrials{iSub,rwd})
                 %number of good stim trials in rwd:
                 %sum(stimTrials{iSub,rwd})
                 %total good stim trials, across rwd/arousal: 
                 %size(arousalTrials{arousalType,iSub,1}) =
                 %size(arousalTrials{arousalType,iSub,2}) =
                 %sum(stimTrials{iSub,1}) + sum(stimTrials{iSub,2})
                 % from concatenated rwd of good trials, which trials are this arousal:
                 % arousalTrials{arousalType,iSub,arousal}
                 % number of arousal trials in rwd:
                 % sum(arousalTrials{arousalType,iSub,arousal}(1:length(stimTrials{iSub,1})
                 %
                 %sum(stimTrials{iSub,rwd});%number of trials for this rwd
        end
        ntrials = length(arousalTrials{arousalType,iSub,1}+arousalTrials{arousalType,iSub,2});
        ind(1) = 0;
        ind(2) = sum(arousalTrials{arousalType,iSub,1});
        ind(3) = ntrials;
        for iperm=1:nperms
            randorder = randperm(ntrials);

            for arousal=1:2
                permArousalTrials{arousalType,iSub,arousal}(iperm,:) = randorder(1+ind(arousal):ind(arousal+1));
                permMeanBasePupil(arousalType,iSub,arousal,iperm) = nanmean(subBaseStimPupil{iSub}(permArousalTrials{arousalType,iSub,arousal}(iperm,:)));
                permMeanStdPupil(arousalType,iSub,arousal,iperm) = nanmean(subStdStimPupil{iSub}(permArousalTrials{arousalType,iSub,arousal}(iperm,:)));
                permMeanStdPulse(arousalType,iSub,arousal,iperm) = nanmean(subStdStimPulse{iSub}(permArousalTrials{arousalType,iSub,arousal}(iperm,:)));
            end
        end
        temp = squeeze(permMeanBasePupil(arousalType,iSub,1,:))' - squeeze(permMeanBasePupil(arousalType,iSub,2,:));
        permDiffBasePupil(arousalType,iSub,:) = temp(:);
        temp = squeeze(permMeanStdPupil(arousalType,iSub,1,:))' - squeeze(permMeanStdPupil(arousalType,iSub,2,:));
        permDiffStdPupil(arousalType,iSub,:) = temp(:);
        temp = squeeze(permMeanStdPulse(arousalType,iSub,1,:))' - squeeze(permMeanStdPulse(arousalType,iSub,2,:));
        permDiffStdPulse(arousalType,iSub,:) = temp(:);
    end
    diffBasePupil(arousalType,:) = squeeze(meanBasePupil(arousalType,:,1)) - squeeze(meanBasePupil(arousalType,:,2));
    diffStdPupil(arousalType,:) = squeeze(meanStdPupil(arousalType,:,1)) - squeeze(meanStdPupil(arousalType,:,2));
    diffStdPulse(arousalType,:) = squeeze(meanStdPulse(arousalType,:,1)) - squeeze(meanStdPulse(arousalType,:,2));
    
    pvalPermBasePupil(arousalType) = sum(squeeze(nanmean(permDiffBasePupil(arousalType,:,:),2)) >= nanmean(diffBasePupil(arousalType,:)))/(nperms^2);
    pvalPermStdPupil(arousalType) = sum(squeeze(nanmean(permDiffStdPupil(arousalType,:,:),2)) >= nanmean(diffStdPupil(arousalType,:)))/(nperms^2);
    pvalPermStdPulse(arousalType) = sum(squeeze(nanmean(permDiffStdPulse(arousalType,:,:),2)) >= nanmean(diffStdPulse(arousalType,:)))/(nperms^2);
    
    plotColors = { [1 0 0], [0 0 1], [0 1 0], [0.5 1 0.2], [1 0.2 0.5] };
    lineStyles = {'-','--'};
    ifig=ifig+1; figure(ifig); clf
    rows=1;
    cols=7;
    subplot(rows,cols,1:2)

    markerSize=60;
%     for featureType=featureTypes
%         for arousal=1:2
%             plot(squeeze(nanmean(classAccFreq{arousalType,featureType,arousal})),lineStyles{arousal},'color',plotColors{featureType})
%             hold on
%         end
%         for ibin=1:nbins
%             for arousal=1:2
%                 maxAcc(arousal) = max(squeeze(nanmean(classAccFreq{arousalType,featureType,arousal}(:,iRoi,ibin))));
%                 %             maxAcc = max(squeeze(mean(accFreq{featureType,1}(:,iRoi,ibin))),squeeze(mean(accFreq{featureType,2}(:,iRoi,ibin))));
%             end
%             if pvalFreq(arousalType,featureType,iRoi,ibin)<pvalThresh
%                 
%                 if maxAcc(1)>maxAcc(2)
%                     scatter(ibin,max(maxAcc),markerSize,plotColors{featureType},'filled');
%                 else
%                     scatter(ibin,max(maxAcc),markerSize,plotColors{featureType});
%                 end
%             end
%             if pvalPermFreq(arousalType,featureType,iRoi,ibin)<pvalThresh
%                 scatter(ibin,max(maxAcc),markerSize,plotColors{featureType},'v','filled');
%             elseif pvalPermFreq(arousalType,featureType,iRoi,ibin)>1-pvalThresh
%                 scatter(ibin,max(maxAcc),markerSize,plotColors{featureType},'v');
%             end
%         end
%         %     mean(squeeze(mean(classAccFreq{featureType})))
%     end
    
    for featureType=3
        for arousal=1:2
            plot(squeeze(nanmean(classAccFreq{arousalType,featureType,arousal})),'color',plotColors{arousal})
            hold on
        end
        for ibin=1:nbins
            for arousal=1:2
                maxAcc(arousal) = max(squeeze(nanmean(classAccFreq{arousalType,featureType,arousal}(:,iRoi,ibin))));
                %             maxAcc = max(squeeze(mean(accFreq{featureType,1}(:,iRoi,ibin))),squeeze(mean(accFreq{featureType,2}(:,iRoi,ibin))));
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
    legend('high','low');
%     for featureType=featureTypes
%         for arousal=1:2
%             plot(squeeze(nanmean(classAccContrast{arousalType,featureType,arousal})),lineStyles{arousal},'color',plotColors{featureType})
%         end
%         for ibin=1:nbins
%             for arousal=1:2
%                 maxAcc(arousal) = max(squeeze(mean(classAccContrast{arousalType,featureType,arousal}(:,iRoi,ibin))));
%                 %             maxAcc = max(squeeze(mean(accFreq{featureType,1}(:,iRoi,ibin))),squeeze(mean(accFreq{featureType,2}(:,iRoi,ibin))));
%             end
%             if pvalContrast(arousalType,featureType,iRoi,ibin)<pvalThresh
%                 if maxAcc(1)>maxAcc(2)
%                     scatter(ibin,max(maxAcc),markerSize,plotColors{featureType},'filled');
%                 else
%                     scatter(ibin,max(maxAcc),markerSize,plotColors{featureType});
%                 end
%             end
%             if pvalPermContrast(arousalType,featureType,iRoi,ibin)<pvalThresh
%                 scatter(ibin,max(maxAcc),markerSize,plotColors{featureType},'v','filled');
%             elseif pvalPermContrast(arousalType,featureType,iRoi,ibin)>1-pvalThresh
%                 scatter(ibin,max(maxAcc),markerSize,plotColors{featureType},'v');
%             end
%         end
%         %     mean(squeeze(mean(classAccContrast{featureType})))
%     end
    
    for featureType=3
        for arousal=1:2
            plot(squeeze(nanmean(classAccContrast{arousalType,featureType,arousal})),'color',plotColors{arousal})
        end
        for ibin=1:nbins
            for arousal=1:2
                maxAcc(arousal) = max(squeeze(mean(classAccContrast{arousalType,featureType,arousal}(:,iRoi,ibin))));
                %             maxAcc = max(squeeze(mean(accFreq{featureType,1}(:,iRoi,ibin))),squeeze(mean(accFreq{featureType,2}(:,iRoi,ibin))));
            end
            if pvalContrast(arousalType,featureType,iRoi,ibin)<pvalThresh
                if maxAcc(1)>maxAcc(2)
                    scatter(ibin,max(maxAcc)+sigHeightT,markerSize,plotColors{1});
                else
                    scatter(ibin,max(maxAcc)+sigHeightT,markerSize,plotColors{2});
                end
            end
            if pvalPermContrast(arousalType,featureType,iRoi,ibin)<pvalThresh
                scatter(ibin,max(maxAcc)+sigHeight5,markerSize,plotColors{1},'*','filled');
                if pvalPermContrast(arousalType,featureType,iRoi,ibin)<0.01
                    scatter(ibin,max(maxAcc)+sigHeight1,markerSize,plotColors{1},'*','filled');
                end
            elseif pvalPermContrast(arousalType,featureType,iRoi,ibin)>1-pvalThresh
                scatter(ibin,max(maxAcc)+sigHeight5,markerSize,plotColors{2},'*');
                if pvalPermContrast(arousalType,featureType,iRoi,ibin)>1-0.01
                    scatter(ibin,max(maxAcc)+sigHeight1,markerSize,plotColors{2},'*');
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
%     for arousal=1:2
%         [max(squeeze(nanmean(classAccContrast{arousalType,featureType,arousal},1))), max(squeeze(nanmean(classAccFreq{arousalType,featureType,arousal},1)))]
%     end
    
    
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
    
    
    %number of trials for each rwd condition,by arousal
    subplot(rows,cols,5)
%     rwdConds(arousalType,iSub,arousal,rwd);
    for arousal=1:2
        bar((arousal-1)*2+[1:2],squeeze(mean(rwdConds(arousalType,:,arousal,:))),'facecolor',plotColors{arousal},'linestyle',lineStyles{arousal});
        
        %         bar((arousal-1)*2+[1:2],squeeze(mean(contrastCondsTrue(arousalType,featureType,:,arousal,:))),'facecolor',plotColors{arousal},'linestyle',lineStyles{arousal});
        hold on
        er=errorbar((arousal-1)*ncontrasts+[1:ncontrasts],squeeze(mean(rwdConds(arousalType,:,arousal,:))),squeeze(std(rwdConds(arousalType,:,arousal,:))));
        %         er=errorbar((arousal-1)*ncontrasts+[1:ncontrasts],squeeze(mean(contrastCondsTrue(arousalType,featureType,:,arousal,:))),squeeze(std(contrastCondsTrue(arousalType,featureType,:,arousal,:))));
        er.Color = [0 0 0];
        er.LineStyle = 'none';
        %         histogram(squeeze(mean(contrastCondsTrue(:,arousal,:))),ncontrasts+1,'facecolor',plotColors{arousal},'linestyle',lineStyles{arousal});
    end
    
        
    %plot mean pupil timecourse for high and low arousal
    subplot(rows,cols,6)
    for arousal=1:2
       plot(squeeze(nanmean(meanPupil(arousalType,:,arousal,:))),'color',plotColors{arousal}); 
       hold on
    end
    title('pupil')
%     
%     %plot mean pupil STD for high and low arousal
%     subplot(rows,cols,6)
%     
    %plot mean pulse timecourse for high and low arousal
    subplot(rows,cols,7)
    for arousal=1:2
        plot(squeeze(nanmean(meanPulse(arousalType,:,arousal,:))),'color',plotColors{arousal});
        hold on
    end
    title('pulse')
    
    set(gcf,'position',[100 100 1100 400])
end

 pvalPermBasePupil 
 pvalPermStdPupil
 pvalPermStdPulse

%uses file created by saveDecoding_arousal_perm.m
close all

clear all
npermsUnivariate=10000;
onlyCorrect=1;%0=all trials, 1=only correct, -1=only incorrect, 2=valid response
for onlyCorrect=[0 2]
    subtractMean = 0;
    
    % for onlyCorrect=[0 1 2]
    includeControl=0;%left and right DMN trial amplitude
    includePulse = 0;
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
        'controlStimBetasMedian',...
        'subBinBetasStim','subBetasStimPupil','subBetasStimPulse',...
        'stimFreqTrials','stimContrastTrials',...
        'subMeanStimPulse');
    clear permArousalTrials
    npermsDecoding=nperms;
    nperms = npermsUnivariate;
    
    tic
    %%
    ifig=0;
    iRoi=1;
    
    for iSub=1:length(goodSubs)
        for ibin=1:nbins
            subBinBetas{iSub,iRoi,ibin} = nanmean(subBinBetasStim{iSub,iRoi,ibin},1);
        end
    end
    
    for arousalType=arousalTypes
        %permutation test for difference in pupil baseline, pupil std, pulse
        %std,
        for iSub=1:length(goodSubs)
            for arousal=1:2
                meanBasePupil(arousalType,iSub,arousal) = nanmean(subBaseStimPupil{iSub}(arousalTrials{arousalType,iSub,arousal}));
                meanStdPupil(arousalType,iSub,arousal) = nanmean(subStdStimPupil{iSub}(arousalTrials{arousalType,iSub,arousal}));
                meanStdPulse(arousalType,iSub,arousal) = nanmean(subStdStimPulse{iSub}(arousalTrials{arousalType,iSub,arousal}));
                meanMeanPulse(arousalType,iSub,arousal) = nanmean(subMeanStimPulse{iSub}(arousalTrials{arousalType,iSub,arousal}));
                
                meanPupil(arousalType,iSub,arousal,:) = nanmean(stimPupilConcat{iSub}(arousalTrials{arousalType,iSub,arousal},:));
                meanPulse(arousalType,iSub,arousal,:) = nanmean(subStimPulse{iSub}(:,arousalTrials{arousalType,iSub,arousal}),2);
                
                meanBetasPupil(arousalType,iSub,arousal) = nanmean(subBetasStimPupil{iSub}(arousalTrials{arousalType,iSub,arousal}));
                meanBetasPulse(arousalType,iSub,arousal) = nanmean(subBetasStimPulse{iSub}(arousalTrials{arousalType,iSub,arousal}));
                
                for ibin=1:nbins
                    for icontrast=1:ncontrasts
                        subBinContrastBetas(arousalType,iSub,ibin,arousal,icontrast) = ...
                            nanmean(subBinBetas{iSub,iRoi,ibin}(stimContrastTrials{iSub}(:)==icontrast & arousalTrials{arousalType,iSub,arousal}(:)==1));
                        subBinContrastBetasStd(arousalType,iSub,ibin,arousal,icontrast) = ...
                            nanstd(subBinBetas{iSub,iRoi,ibin}(stimContrastTrials{iSub}(:)==icontrast & arousalTrials{arousalType,iSub,arousal}(:)==1));
                        
                    end
                    for ifreq=1:nfreqs
                        subBinFreqBetas(arousalType,iSub,ibin,arousal,ifreq) = ...
                            nanmean(subBinBetas{iSub,iRoi,ibin}(stimFreqTrials{iSub}(:)==ifreq & arousalTrials{arousalType,iSub,arousal}(:)==1));
                        subBinFreqBetasStd(arousalType,iSub,ibin,arousal,ifreq) = ...
                            nanstd(subBinBetas{iSub,iRoi,ibin}(stimFreqTrials{iSub}(:)==ifreq & arousalTrials{arousalType,iSub,arousal}(:)==1));
                        
                    end
                    %                 subBinMeanBetas(arousalType,iSub,ibin,arousal) = nanmean(subBinBetas{iSub,iRoi,ibin}(arousalTrials{arousalType,iSub,arousal}));
                end
                
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
                    permMeanMeanPulse(arousalType,iSub,arousal,iperm) = nanmean(subMeanStimPulse{iSub}(permArousalTrials{arousalType,iSub,arousal}(iperm,:)));
                    permMeanBetasPupil(arousalType,iSub,arousal,iperm) = nanmean(subBetasStimPupil{iSub}(permArousalTrials{arousalType,iSub,arousal}(iperm,:)));
                    permMeanBetasPulse(arousalType,iSub,arousal,iperm) = nanmean(subBetasStimPulse{iSub}(permArousalTrials{arousalType,iSub,arousal}(iperm,:)));
                    %                 for ibin=1:nbins
                    %                     permBinMeanBetas(arousalType,iSub,ibin,arousal,iperm) = nanmean(subBinBetas{iSub,iRoi,ibin}(permArousalTrials{arousalType,iSub,arousal}(iperm,:)));
                    %                 end
                    %                 permArousalTrialsBinary = zeros(size(arousalTrials{arousalType,iSub,1}));
                    %                 permArousalTrialsBinary(permArousalTrials{arousalType,iSub,arousal}(iperm,:)) = ones;
                    temp = arousalTrials{arousalType,iSub,arousal}(:)';
                    for ibin=1:nbins
                        
                        for icontrast=1:ncontrasts
                            
                            permBinContrastBetas(arousalType,iSub,ibin,arousal,icontrast,iperm)= ...
                                nanmean(subBinBetas{iSub,iRoi,ibin}(stimContrastTrials{iSub}==icontrast...
                                & temp(randorder)==1));
                            permBinContrastBetasStd(arousalType,iSub,ibin,arousal,icontrast,iperm)= ...
                                nanstd(subBinBetas{iSub,iRoi,ibin}(stimContrastTrials{iSub}==icontrast...
                                & temp(randorder)==1));
                        end
                        for ifreq=1:nfreqs
                            permBinFreqBetas(arousalType,iSub,ibin,arousal,ifreq,iperm)= ...
                                nanmean(subBinBetas{iSub,iRoi,ibin}(stimFreqTrials{iSub}==ifreq...
                                & temp(randorder)==1));
                            permBinFreqBetasStd(arousalType,iSub,ibin,arousal,ifreq,iperm)= ...
                                nanstd(subBinBetas{iSub,iRoi,ibin}(stimFreqTrials{iSub}==ifreq...
                                & temp(randorder)==1));
                        end
                    end
                end
            end
            
            %         permBinContrastBetasDiff(arousalType,iSub,:,:,1,:) = permBinContrastBetas(arousalType,iSub,:,:,1,:)-permBinContrastBetas(arousalType,iSub,:,:,2,:);
            permBinContrastBetasDiff(arousalType,iSub,:,:,1,:) = abs(permBinContrastBetas(arousalType,iSub,:,:,1,:)-permBinContrastBetas(arousalType,iSub,:,:,2,:));
            permBinDiffContrastBetasDiff(arousalType,iSub,:,:) = squeeze(permBinContrastBetasDiff(arousalType,iSub,:,1,1,:) - permBinContrastBetasDiff(arousalType,iSub,:,2,1,:));
            
            permBinMeanContrastBetas(arousalType,iSub,:,:,:) = mean(permBinContrastBetas(arousalType,iSub,:,:,:,:),5);
            permBinDiffMeanContrastBetas(arousalType,iSub,:,:) = squeeze(permBinMeanContrastBetas(arousalType,iSub,:,1,:) - permBinMeanContrastBetas(arousalType,iSub,:,2,:));
            permBinMeanFreqBetas(arousalType,iSub,:,:,:) = mean(permBinFreqBetas(arousalType,iSub,:,:,:,:),5);
            permBinDiffMeanFreqBetas(arousalType,iSub,:,:) = squeeze(permBinMeanFreqBetas(arousalType,iSub,:,1,:) - permBinMeanContrastBetas(arousalType,iSub,:,2,:));
            
            %
            for ibin=1:nbins
                for arousal=1:2
                    for iperm=1:nperms
                        temp = squeeze(permBinFreqBetas(arousalType,iSub,ibin,arousal,:,iperm));
                        allDiff = temp - temp';
                        trilOnes = tril(ones(size(allDiff)),-1);
                        %                     permBinFreqBetasDiff(arousalType,iSub,ibin,arousal,iperm) = mean(allDiff(trilOnes>0));
                        permBinFreqBetasDiff(arousalType,iSub,ibin,arousal,iperm) = mean(abs(allDiff(trilOnes>0)));
                    end
                end
            end
            %difference between mean beta-difference, between arousals
            permBinDiffFreqBetasDiff(arousalType,iSub,:,:) = squeeze(permBinFreqBetasDiff(arousalType,iSub,:,1,:) - permBinFreqBetasDiff(arousalType,iSub,:,2,:));
            
            permBinDiffContrastBetasStd(arousalType,iSub,:,:) = squeeze(mean(permBinContrastBetasStd(arousalType,iSub,:,1,:,:),5) - mean(permBinContrastBetasStd(arousalType,iSub,:,2,:,:),5));%difference between arousals for mean std
            permBinDiffFreqBetasStd(arousalType,iSub,:,:) = squeeze(mean(permBinFreqBetasStd(arousalType,iSub,:,1,:,:),5) - mean(permBinFreqBetasStd(arousalType,iSub,:,2,:,:),5));%difference between arousals for mean std
            
            %          permBinContrastBetasDiff(arousalType,iSub,ibin,arousal,1,iperm);
            %         permBinContrastCnr(arousalType,iSub,ibin,arousal,iperm) = squeeze(permBinContrastBetasDiff(arousalType,iSub,ibin,arousal,1,:))./permBinContrastBetasStd(arousalType,iSub,ibin,arousal,icontrast,iperm)
            
            %         permDiffContrastCnr(arousalType,iSub,:,:) = squeeze(permBinContrastCnr(arousalType,iSub,:,1,:)-permBinContrastCnr(arousalType,iSub,:,2,:));
            %         permDiffFreqCnr(arousalType,iSub,:,:) = squeeze(permBinFreqCnr(arousalType,iSub,:,1,:)-permBinFreqCnr(arousalType,iSub,:,2,:));
            %         permBinContrastBetas(arousalType,iSub,ibin,arousal,icontrast,iperm)
            
            permDiffBasePupil(arousalType,iSub,:) = squeeze(permMeanBasePupil(arousalType,iSub,1,:) - permMeanBasePupil(arousalType,iSub,2,:));
            permDiffStdPupil(arousalType,iSub,:) = squeeze(permMeanStdPupil(arousalType,iSub,1,:) - permMeanStdPupil(arousalType,iSub,2,:));
            permDiffStdPulse(arousalType,iSub,:) = squeeze(permMeanStdPulse(arousalType,iSub,1,:) - permMeanStdPulse(arousalType,iSub,2,:));
            permDiffMeanPulse(arousalType,iSub,:) = squeeze(permMeanMeanPulse(arousalType,iSub,1,:) - permMeanMeanPulse(arousalType,iSub,2,:));
            permDiffBetasPupil(arousalType,iSub,:) = squeeze(permMeanBetasPupil(arousalType,iSub,1,:) - permMeanBetasPupil(arousalType,iSub,2,:));
            permDiffBetasPulse(arousalType,iSub,:) = squeeze(permMeanBetasPulse(arousalType,iSub,1,:) - permMeanBetasPulse(arousalType,iSub,2,:));
            
        end
        
        
        
        %mean difference between freq betas, for each arousal
        for ibin=1:nbins
            for arousal=1:2
                diffBinFreqBetas(arousalType,iSub,ibin,arousal) = nanstd(squeeze(subBinFreqBetas(arousalType,iSub,ibin,arousal,:)));
                %             diffBinFreqBetas(arousalType,iSub,ibin,arousal) = nanmean(abs(squeeze(subBinFreqBetas(arousalType,iSub,ibin,arousal,:))));
            end
        end
        
        
        %mean contrast std for each arousal
        diffBinContrastBetasStd(arousalType,:,:) = squeeze(nanmean(subBinContrastBetasStd(arousalType,:,:,1,:),5)-nanmean(subBinContrastBetasStd(arousalType,:,:,2,:),5));
        
        %mean freq std for each arousal
        diffBinFreqBetasStd(arousalType,:,:) = squeeze(nanmean(subBinFreqBetasStd(arousalType,:,:,1,:),5)-nanmean(subBinFreqBetasStd(arousalType,:,:,2,:),5));
        
        
        %     subBinContrastBetas(arousalType,iSub,ibin,arousal,icontrast)
        temp = squeeze(nanmean(subBinFreqBetasStd(arousalType,:,:,:,:)));%average over subjects
        meanFreqVar(arousalType,:,:) = squeeze(nanmean(temp,3));%average over freqs.
        temp = squeeze(nanmean(subBinContrastBetasStd(arousalType,:,:,:,:)));%average over subjects
        meanContrastVar(arousalType,:,:) = squeeze(nanmean(temp,3));%average over contrasts.
        
        diffBasePupil(arousalType,:) = squeeze(meanBasePupil(arousalType,:,1)) - squeeze(meanBasePupil(arousalType,:,2));
        diffStdPupil(arousalType,:) = squeeze(meanStdPupil(arousalType,:,1)) - squeeze(meanStdPupil(arousalType,:,2));
        diffStdPulse(arousalType,:) = squeeze(meanStdPulse(arousalType,:,1)) - squeeze(meanStdPulse(arousalType,:,2));
        diffMeanPulse(arousalType,:) = squeeze(meanMeanPulse(arousalType,:,1)) - squeeze(meanMeanPulse(arousalType,:,2));
        diffBetasPupil(arousalType,:) = squeeze(meanBetasPupil(arousalType,:,1)) - squeeze(meanBetasPupil(arousalType,:,2));
        diffBetasPulse(arousalType,:) = squeeze(meanBetasPulse(arousalType,:,1)) - squeeze(meanBetasPulse(arousalType,:,2));
        
        
        
        pvalPermBasePupil(arousalType) = sum(squeeze(nanmean(permDiffBasePupil(arousalType,:,:),2)) >= nanmean(diffBasePupil(arousalType,:)))/nperms;
        pvalPermStdPupil(arousalType) = sum(squeeze(nanmean(permDiffStdPupil(arousalType,:,:),2)) >= nanmean(diffStdPupil(arousalType,:)))/nperms;
        pvalPermStdPulse(arousalType) = sum(squeeze(nanmean(permDiffStdPulse(arousalType,:,:),2)) >= nanmean(diffStdPulse(arousalType,:)))/nperms;
        pvalPermMeanPulse(arousalType) = sum(squeeze(nanmean(permDiffMeanPulse(arousalType,:,:),2)) >= nanmean(diffMeanPulse(arousalType,:)))/nperms;
        pvalPermBetasPupil(arousalType) = sum(squeeze(nanmean(permDiffBetasPupil(arousalType,:,:),2)) >= nanmean(diffBetasPupil(arousalType,:)))/nperms;
        pvalPermBetasPulse(arousalType) = sum(squeeze(nanmean(permDiffBetasPulse(arousalType,:,:),2)) >= nanmean(diffBetasPulse(arousalType,:)))/nperms;
        
        
        %high-low contrast for each arousal
        diffBinContrastBetas(arousalType,:,:,:) = abs(squeeze(subBinContrastBetas(arousalType,:,:,:,1)-subBinContrastBetas(arousalType,:,:,:,2)));
        
        
        %SNR is the mean betas divided by the mean variance
        subBinContrastMeanBetas(arousalType,:,:,:) = squeeze(mean(abs(subBinContrastBetas(arousalType,:,:,:,:)),5));%average over conditions. sub,bin,arousal
        subContrastVar(arousalType,:,:,:) = squeeze(nanmean(subBinContrastBetasStd(arousalType,:,:,:,:),5));%average over conditions
        subContrastSnr(arousalType,:,:,:) = subBinContrastMeanBetas(arousalType,:,:,:)./subContrastVar(arousalType,:,:,:);%sub,bin,arousal
        subContrastSnrDiff(arousalType,:,:) = squeeze(subContrastSnr(arousalType,:,:,1)-subContrastSnr(arousalType,:,:,2));
        
        %CNR is the difference in betas divided by the mean variance
        subContrastCnr(arousalType,:,:,:) = diffBinContrastBetas(arousalType,:,:,:)./subContrastVar(arousalType,:,:,:);%sub,bin,arousal
        subContrastCnrDiff(arousalType,:,:) = squeeze(subContrastCnr(arousalType,:,:,1)-subContrastCnr(arousalType,:,:,2));
        
        
        %SNR is the mean betas divided by the mean variance
        subBinFreqMeanBetas(arousalType,:,:,:) = squeeze(mean(abs(subBinFreqBetas(arousalType,:,:,:,:)),5));%average over conditions. sub,bin,arousal
        subFreqVar(arousalType,:,:,:) = squeeze(nanmean(subBinFreqBetasStd(arousalType,:,:,:,:),5));%average over conditions
        subFreqSnr(arousalType,:,:,:) = subBinFreqMeanBetas(arousalType,:,:,:)./subFreqVar(arousalType,:,:,:);%sub,bin,arousal
        subFreqSnrDiff(arousalType,:,:) = squeeze(subFreqSnr(arousalType,:,:,1)-subFreqSnr(arousalType,:,:,2));
        
        %CNR is the difference in betas divided by the mean variance
        subFreqCnr(arousalType,:,:,:) = diffBinFreqBetas(arousalType,:,:,:)./subFreqVar(arousalType,:,:,:);%sub,bin,arousal
        subFreqCnrDiff(arousalType,:,:) = squeeze(subFreqCnr(arousalType,:,:,1)-subFreqCnr(arousalType,:,:,2));
        
        %difference between mean betas per arousal, for each bin and subject
        diffBinContrastBetasDiff(arousalType,:,:) = squeeze(diffBinContrastBetas(arousalType,:,:,1) - diffBinContrastBetas(arousalType,:,:,2));%isub,ibin
        diffBinMeanContrastBetas(arousalType,:,:) = squeeze(subBinContrastMeanBetas(arousalType,:,:,1) - subBinContrastMeanBetas(arousalType,:,:,2));
        
        diffBinFreqBetasDiff(arousalType,:,:) = squeeze(diffBinFreqBetas(arousalType,:,:,1) - diffBinFreqBetas(arousalType,:,:,2));
        diffBinMeanFreqBetas(arousalType,:,:) = squeeze(subBinFreqMeanBetas(arousalType,:,:,1) - subBinFreqMeanBetas(arousalType,:,:,2));
        
        
        %high-low contrast, then difference between arousals
        pvalDiffBinContrastBetas(arousalType,:) = sum(squeeze(nanmean(permBinDiffContrastBetasDiff(arousalType,:,:,:),2))' >= squeeze(nanmean(diffBinContrastBetasDiff(arousalType,:,:)))')/nperms;
        %average std across contrast, then difference between arousals
        pvalDiffBinContrastBetasStd(arousalType,:) = sum(squeeze(nanmean(permBinDiffContrastBetasStd(arousalType,:,:,:),2))' >= squeeze(nanmean(diffBinContrastBetasStd(arousalType,:,:)))')/nperms;
        %mean across contrasts, then diff between arousals
        pvalDiffBinMeanContrastBetas(arousalType,:) = sum(squeeze(nanmean(permBinDiffMeanContrastBetas(arousalType,:,:,:),2))' >= squeeze(nanmean(diffBinMeanContrastBetas(arousalType,:,:)))')/nperms;
        
        %std of freq betas, then difference between arousals
        pvalDiffBinFreqBetas(arousalType,:) = sum(squeeze(nanmean(permBinDiffFreqBetasDiff(arousalType,:,:,:),2))' >= squeeze(nanmean(diffBinFreqBetasDiff(arousalType,:,:)))')/nperms;
        %average std across freqs, then difference between arousals
        pvalDiffBinFreqBetasStd(arousalType,:) = sum(squeeze(nanmean(permBinDiffFreqBetasStd(arousalType,:,:,:),2))' >= squeeze(nanmean(diffBinFreqBetasStd(arousalType,:,:)))')/nperms;
        %mean across freqs, then diff between arousals
        pvalDiffBinMeanFreqBetas(arousalType,:) = sum(squeeze(nanmean(permBinDiffMeanFreqBetas(arousalType,:,:,:),2))' >= squeeze(nanmean(diffBinMeanFreqBetas(arousalType,:,:)))')/nperms;
        
        
        %can we use permBinMeanContrastBetas instead?
        %     permBinContrastMeanBetas(arousalType,:,:,:,:) = squeeze(mean(permBinContrastBetas(arousalType,:,:,:,:,:),5));
        %     permBinFreqMeanBetas(arousalType,:,:,:,:) = squeeze(mean(permBinFreqBetas(arousalType,:,:,:,:,:),5));
        %     permBinContrastMeanBetas(arousalType,iSub,ibin,arousal,iperm)
        
        %perm CNR
        permBinContrastCnr(arousalType,:,:,:,:) = squeeze(permBinContrastBetasDiff(arousalType,:,:,:,1,:))./...
            squeeze(nanmean(permBinContrastBetasStd(arousalType,:,:,:,:,:),5));
        permBinFreqCnr(arousalType,:,:,:,:) = squeeze(permBinFreqBetasDiff(arousalType,:,:,:,:))./...
            squeeze(mean(permBinFreqBetasStd(arousalType,:,:,:,:,:),5));
        
        permBinContrastCnrDiff(arousalType,:,:,:) = squeeze(permBinContrastCnr(arousalType,:,:,1,:)-permBinContrastCnr(arousalType,:,:,2,:));
        permBinFreqCnrDiff(arousalType,:,:,:) = squeeze(permBinFreqCnr(arousalType,:,:,1,:)-permBinFreqCnr(arousalType,:,:,2,:));
        
        
        %perm SNR
        permBinContrastSnr(arousalType,:,:,:,:) = squeeze(permBinMeanContrastBetas(arousalType,:,:,:,:))./...
            squeeze(nanmean(permBinContrastBetasStd(arousalType,:,:,:,:,:),5));
        permBinFreqSnr(arousalType,:,:,:,:) = squeeze(permBinMeanFreqBetas(arousalType,:,:,:,:))./...
            squeeze(mean(permBinFreqBetasStd(arousalType,:,:,:,:,:),5));
        
        permBinContrastSnrDiff(arousalType,:,:,:) = squeeze(permBinContrastSnr(arousalType,:,:,1,:)-permBinContrastSnr(arousalType,:,:,2,:));
        permBinFreqSnrDiff(arousalType,:,:,:) = squeeze(permBinFreqSnr(arousalType,:,:,1,:)-permBinFreqSnr(arousalType,:,:,2,:));
        
        
        %CNR pvals
        pvalContrastCnrDiff(arousalType,:) = sum(squeeze(nanmean(permBinContrastCnrDiff(arousalType,:,:,:,:),2))' >= squeeze(nanmean(subContrastCnrDiff(arousalType,:,:),2))')/nperms;
        pvalFreqCnrDiff(arousalType,:) = sum(squeeze(nanmean(permBinFreqCnrDiff(arousalType,:,:,:,:),2))' >= squeeze(nanmean(subFreqCnrDiff(arousalType,:,:),2))')/nperms;
        %     permBinContrastCnr(arousalType,iSub,ibin,arousal,iperm)
        
        %SNR pvals
        pvalContrastSnrDiff(arousalType,:) = sum(squeeze(nanmean(permBinContrastSnrDiff(arousalType,:,:,:,:),2))' >= squeeze(nanmean(subContrastSnrDiff(arousalType,:,:),2))')/nperms;
        pvalFreqSnrDiff(arousalType,:) = sum(squeeze(nanmean(permBinFreqSnrDiff(arousalType,:,:,:,:),2))' >= squeeze(nanmean(subFreqSnrDiff(arousalType,:,:),2))')/nperms;
        
        
    end
    %%
    plotColors = { [1 0 0], [0 0 1], [0 1 0], [0.5 1 0.2], [1 0.2 0.5] };
    lineStyles = {'-','--'};
    markerSize=60;
    rows=1;
    cols=7;
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
        end
        legend('high','low');
        
        
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
            %         er=errorbar((arousal-1)*ncontrasts+[1:ncontrasts],squeeze(mean(contrastCondsTrue(arousalType,featureType,:,arousal,:))),squeeze(std(contrastCondsTrue(arousalType,featureType,:,arousal,:))));
            er.Color = [0 0 0];
            er.LineStyle = 'none';
            %         histogram(squeeze(mean(contrastCondsTrue(:,arousal,:))),ncontrasts+1,'facecolor',plotColors{arousal},'linestyle',lineStyles{arousal});
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
    
    meanContrastCnr = squeeze(mean(subContrastCnr,2));
    meanFreqCnr = squeeze(mean(subFreqCnr,2));
    meanContrastSnr = squeeze(mean(subContrastSnr,2));
    meanFreqSnr = squeeze(mean(subFreqSnr,2));
    
    meanFreqBetasDiff = squeeze(nanmean(diffBinFreqBetas,2));%arousalTYpe,ibin,arousal
    meanContrastBetasDiff = squeeze(nanmean(diffBinContrastBetas,2));
    
    meanContrastMeanBetas = squeeze(nanmean(subBinContrastMeanBetas,2));
    meanFreqMeanBetas = squeeze(nanmean(subBinFreqMeanBetas,2));
    %%
    rows=2;
    cols=6;
    markerSize=40;
    sigHeight1=0.08;
    sigHeight5=0.05;
    for arousalType=arousalTypes
        ifig=ifig+1; figure(ifig); clf
        isubplot=0;
        %plot contrast betas as function of eccentricity for both arousal levels
        isubplot = isubplot+1;
        subplot(rows,cols,isubplot)
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
        
        %plot mean contrast betas as function of eccentricity for both arousal levels
        isubplot = isubplot+1;
        subplot(rows,cols,isubplot)
        for arousal=1:2
            plot(squeeze(meanContrastMeanBetas(arousalType,:,arousal)),'color',plotColors{arousal});
            hold on
        end
        xlim([0 nbins+1]);
        title('mean contrast');
        plotPvals(nbins, squeeze(meanContrastMeanBetas(arousalType,:,:)), squeeze(pvalDiffBinMeanContrastBetas(arousalType,:)),  sigHeight5, sigHeight1,markerSize, plotColors);
        
        
        %plot difference between high and low contrast betas as function of eccentricity for both arousal levels
        isubplot = isubplot+1;
        subplot(rows,cols,isubplot)
        
        for arousal=1:2
            plot(meanContrastBetasDiff(arousalType,:,arousal),'color',plotColors{arousal});
            hold on
        end
        xlim([0 nbins+1]);
        title('high-low contrast');
        plotPvals(nbins, squeeze(meanContrastBetasDiff(arousalType,:,:)), squeeze(pvalDiffBinContrastBetas(arousalType,:)),  sigHeight5, sigHeight1,markerSize, plotColors);
        
        
        %plot contrast variability
        isubplot = isubplot+1;
        subplot(rows,cols,isubplot)
        for arousal=1:2
            plot(squeeze(meanContrastVar(arousalType,:,arousal)),'color',plotColors{arousal});
            hold on
        end
        xlim([0 nbins+1]);
        title('mean variability')
        legend('high','low');
        plotPvals(nbins, squeeze(meanContrastVar(arousalType,:,:)), squeeze(pvalDiffBinContrastBetasStd(arousalType,:)),  sigHeight5, sigHeight1,markerSize, plotColors);
        
        %contrast SNR
        isubplot = isubplot+1;
        subplot(rows,cols,isubplot)
        for arousal=1:2
            plot(squeeze(meanContrastSnr(arousalType,:,arousal)),'color',plotColors{arousal},'linestyle',lineStyles{1});
            hold on
            %         plot(squeeze(mean(subCnr2(arousalType,:,:,arousal),2)),'color',plotColors{arousal},'linestyle',lineStyles{2});
        end
        title('contrast SNR');
        plotPvals(nbins, squeeze(meanContrastSnr(arousalType,:,:)), squeeze(pvalContrastSnrDiff(arousalType,:)),  sigHeight5, sigHeight1,markerSize, plotColors);
        
        
        %contrast CNR
        isubplot = isubplot+1;
        subplot(rows,cols,isubplot)
        for arousal=1:2
            plot(squeeze(meanContrastCnr(arousalType,:,arousal)),'color',plotColors{arousal},'linestyle',lineStyles{1});
            hold on
            %         plot(squeeze(mean(subCnr2(arousalType,:,:,arousal),2)),'color',plotColors{arousal},'linestyle',lineStyles{2});
        end
        title('contrast CNR');
        plotPvals(nbins, squeeze(meanContrastCnr(arousalType,:,:)), squeeze(pvalContrastCnrDiff(arousalType,:)),  sigHeight5, sigHeight1,markerSize, plotColors);
        
        
        
        %plot frequency betas as function of eccentricity for both arousal levels
        isubplot = isubplot+1;
        subplot(rows,cols,isubplot)
        for arousal=1:2
            for ifreq=1:nfreqs
                plot(squeeze(nanmean(subBinFreqBetas(arousalType,:,:,arousal,ifreq))),...
                    'color',plotColors{arousal}+(1-plotColors{arousal})/(1/nfreqs+ifreq));
                hold on
            end
        end
        xlim([0 nbins+1]);
        title('frequency');
        
        %plot mean frequency betas as function of eccentricity for both arousal levels
        isubplot = isubplot+1;
        subplot(rows,cols,isubplot)
        for arousal=1:2
            plot(squeeze(meanFreqMeanBetas(arousalType,:,arousal)),'color',plotColors{arousal});
            hold on
        end
        plotPvals(nbins, squeeze(meanFreqMeanBetas(arousalType,:,:)), squeeze(pvalDiffBinMeanFreqBetas(arousalType,:)),  sigHeight5, sigHeight1,markerSize, plotColors);
        xlim([0 nbins+1]);
        title('mean freq');
        
        
        %plot difference between frequency betas as function of eccentricity for both arousal levels
        isubplot = isubplot+1;
        subplot(rows,cols,isubplot)
        for arousal=1:2
            plot(meanFreqBetasDiff(arousalType,:,arousal),'color',plotColors{arousal});
            hold on
        end
        title('diff between freqs');
        plotPvals(nbins, squeeze(meanFreqBetasDiff(arousalType,:,:)), squeeze(pvalDiffBinFreqBetas(arousalType,:)),  sigHeight5, sigHeight1,markerSize, plotColors);
        
        
        
        %plot frequency variability
        isubplot = isubplot+1;
        subplot(rows,cols,isubplot)
        for arousal=1:2
            plot(squeeze(meanFreqVar(arousalType,:,arousal)),'color',plotColors{arousal});
            hold on
        end
        xlim([0 nbins+1]);
        title('mean variability')
        plotPvals(nbins, squeeze(meanFreqVar(arousalType,:,:)), squeeze(pvalDiffBinFreqBetasStd(arousalType,:)),  sigHeight5, sigHeight1,markerSize, plotColors);
        
        
        
        %plot freq SNR
        isubplot = isubplot+1;
        subplot(rows,cols,isubplot)
        
        for arousal=1:2
            plot(squeeze(meanFreqSnr(arousalType,:,arousal)),'color',plotColors{arousal},'linestyle',lineStyles{1});
            hold on
        end
        title('freq SNR');
        plotPvals(nbins, squeeze(meanFreqSnr(arousalType,:,:)), squeeze(pvalFreqSnrDiff(arousalType,:)),  sigHeight5, sigHeight1,markerSize, plotColors);
        
        
        %plot freq CNR
        isubplot = isubplot+1;
        subplot(rows,cols,isubplot)
        
        for arousal=1:2
            plot(squeeze(meanFreqCnr(arousalType,:,arousal)),'color',plotColors{arousal},'linestyle',lineStyles{1});
            hold on
            %         plot(squeeze(mean(subCnr2(arousalType,:,:,arousal),2)),'color',plotColors{arousal},'linestyle',lineStyles{2});
        end
        title('freq CNR');
        plotPvals(nbins, squeeze(meanFreqCnr(arousalType,:,:)), squeeze(pvalFreqCnrDiff(arousalType,:)),  sigHeight5, sigHeight1,markerSize, plotColors);
        
        
        set(gcf,'position',[70 80 1100 400])
    end
    
    %%
    pvalPermBasePupil
    pvalPermStdPupil
    pvalPermStdPulse
    pvalPermMeanPulse
    pvalPermBetasPupil
    pvalPermBetasPulse
    
    pvalDiffBinContrastBetas
    pvalDiffBinContrastBetasStd
    
    
    save([dataFolder 'univariateArousalPerm' onlyCorrectStr zscoreStr controlStr pulseStr pupilStr pupilBaseStr subtractMeanStr '.mat'],...
        'npermsDecoding','npermsUnivariate',...
        'arousalTypes','meanPupil','plotColors','pvalPermBetasPupil','meanPulse',...
        'pvalPermBetasPulse','subBinContrastBetas','subBinFreqBetas',...
        'pvalDiffBinContrastBetas','diffBinContrastBetas','meanContrastVar',...
        'pvalDiffBinFreqBetas','meanFreqVar',...
        'pvalDiffBinFreqBetasStd','pvalDiffBinContrastBetasStd',...
        'subBinContrastBetasStd','subBinFreqBetasStd','subContrastCnr','subFreqCnr',...
        'pvalPermBasePupil','pvalPermStdPupil','pvalPermStdPulse','pvalPermMeanPulse',...
        'pvalPermBetasPupil',...
        'rwdConds',...
        'pvalContrastCnrDiff','pvalFreqCnrDiff','meanContrastCnr','meanFreqCnr',...
        'pvalContrastSnrDiff','pvalFreqSnrDiff','meanContrastSnr','meanFreqSnr',...
        'meanContrastBetasDiff','pvalDiffBinContrastBetas','meanContrastMeanBetas','pvalDiffBinMeanContrastBetas',...
        'meanFreqBetasDiff','pvalDiffBinFreqBetas','meanFreqMeanBetas','pvalDiffBinMeanFreqBetas');
    
    
    toc
    % end
end

% [ pvalPermBasePupil;  pvalPermStdPupil;  pvalPermStdPulse; pvalPermBetasPupil; pvalPermBetasPulse]

function [] = plotPvals(nbins, data, pvals, sigHeight5, sigHeight1,markerSize, plotColors)
markerSize=30;
sigHeight5 = (max(data(:)) - min(data(:)))*0.15;
sigHeight1 = (max(data(:)) - min(data(:)))*0.22;
for ibin=1:nbins
    if pvals(ibin)<0.05
        scatter(ibin,data(ibin,1) + sigHeight5,markerSize,plotColors{1},'*');
        if pvals(ibin)<0.01
            scatter(ibin,data(ibin,1) + sigHeight1,markerSize,plotColors{1},'*');
        end
    elseif pvals(ibin)>0.95
        scatter(ibin,data(ibin,2) + sigHeight5,markerSize,plotColors{2},'*');
        if pvals(ibin)>0.99
            scatter(ibin,data(ibin,2) + sigHeight1,markerSize,plotColors{2},'*');
        end
    end
end
end
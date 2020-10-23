%uses files saved by savePhysio.m, saveGroupEyeData.m, saveRoiTC.m

close all
mrQuit
clear all


% fillMissing = 1;
bensonROIs = [1:3];
eccMin = 1;
eccMax = 70 +eccMin;
nbins = 12;
binBorders = logspace(log10(eccMin),log10(eccMax),nbins+1);
nbins = length(binBorders)-1;
for i=2:length(binBorders)
    binCenters(i-1) = (binBorders(i)+binBorders(i-1))/2;
end
binBorders = binBorders - eccMin;


nfreqs=5;
ncontrasts=2;
tic

% dataFolder = '/Volumes/MH02086153MACDT-Drobo/rwdFmri/';
dataFolder = '/Users/rothzn/rwdFmri/';
% dataFolder = 'c:\rwdFmri\';
freqs = logspace(-0.3,0.5,5);
contrasts = logspace(-0.7,0,2);



onlyCorrect=0;%0=all trials, 1=only correct, -1=only incorrect, 2=valid response
for subtractMean = 0:1%subtract mean null trial from all trials
    for toZscore=0:1 %0 or 1
        for concatProj= 0:1
            for onlyCorrect=[0 1 2]
                
                %Load eyetracking data
                load([dataFolder 'rwdRapidEyeData.mat'], 'subFolders', 'samplerate',  ...
                    'trialsPerRun', 'trialLength', ...
                    'numSubs', ...
                    'sacRwd','binnedSac', 'smoothSac', ...
                    'rwdPupil','meanPupil','taskTimes',...
                    'sacRun','pupilRun', 'numRuns','startTimes','endTimes','eyetrackerTime',...
                    'nullTrials','nullPupil','stimPupil','meanNullPupil','meanStimPupil',...
                    'sacRunSmooth', 'sacRunSmoothTrialZeroFilled','sacRunSmoothTrial','L',...
                    'trialCorrectness', 'trialResponse','trialCorrectResponse');
                eyeSubFolders = subFolders;
                eyeNullTrials = nullTrials;
                origNullPupil = nullPupil;
                origStimPupil = stimPupil;
                pupilTrialTime = 7500;
                for iSub=1:size(nullTrials,1)
                    for rwd=1:2
                        stimTrials{iSub,rwd} = 1-nullTrials{iSub,rwd};
                    end
                end
                eyeStimTrials = stimTrials;
                clear stimTrials
                baselineT = 25;
                
                
                %Load physiology
                ConcatProjStr = '';
                if concatProj
                    ConcatProjStr = 'concatProj';
                end
                load([dataFolder 'rwdRapid_physio.mat'], 'concatInfo',  ...
                    'subFolders','trialLength',...
                    'ecgselect','ecgSampleRate','ecgTrial','ecgRunLength','ecgInterpMethod',...
                    'ecg','ecgPulseRate','interpPulseRate','ecgPpgAmp',...
                    'respselect','resp',...
                    'rwdPulseTrials','rwdPpgTrials','rwdRvTrials','meanPulse','meanPpg','meanRV',...
                    'designMatPulse','designMatResp','designMatRespPulse','designMatPpg','designMatRespPulsePpg');
                
                
                %Load fMRI data
                zScoreString = '';
                if toZscore
                    zScoreString = '_zscored';
                end
                concatProjStr = '';
                if concatProj
                    concatProjStr = 'proj';
                end
                load([dataFolder 'roiTC_' zScoreString concatProjStr '.mat'], 'subFolders', 'roiNames', ...
                    'numRuns','numTRs','concatInfo',...
                    'frames', 'junkedFrames', 'TR', 'trialsPerRun', 'trialLength', 'nVolumes',...
                    'eccen','ang','areas',...
                    'numVox','roiTC','nullTrialsTRs',...
                    'nullTrialsRun','nullTrials','contrastTrialsRun','freqTrialsRun','contrastTrials','freqTrials',...
                    'nullTseries','nullTrialTseries','stimTseries','stimTrialTseries');
                
                % goodSubs = [1:12 14];
                goodSubs = [1:14];
                % goodSubs = 1:2;
                %%
                nullPupilRwd = origNullPupil;
                stimPupilRwd = origStimPupil;
                for iSub=1:length(goodSubs)
                    
                    for rwd=1:2
                        %find trials with correct type of response
                        temp = trialCorrectness{goodSubs(iSub),rwd}(:,2:end);%first trial of each run is junked, as in fMRI data
                        rwdTrialCorrectness{iSub,rwd} = temp(:);%convert all trials to one long vector
                        allTrialResponse = [];
                        for r=1:size(trialResponse,3)
                            allTrialResponse = [allTrialResponse; trialResponse{iSub,rwd,r}(2:end)];
                        end
                        rwdTrialResponse{iSub,rwd} = allTrialResponse;
                        
                        switch onlyCorrect
                            case 0 %all trials
                                correctTrials{iSub,rwd} = ones(size(rwdTrialCorrectness{iSub,rwd}));
                            case 1 %only correct trials
                                correctTrials{iSub,rwd} = rwdTrialCorrectness{iSub,rwd}==1;
                            case -1 %only incorrect trials
                                correctTrials{iSub,rwd} = rwdTrialCorrectness{iSub,rwd}==0 & rwdTrialResponse{iSub,rwd}>0;
                            case 2 %all trials with valid response
                                correctTrials{iSub,rwd} = rwdTrialResponse{iSub,rwd}>0;
                        end
                        %convert trial flags to TR flags
                        correctTrials{iSub,rwd} = correctTrials{iSub,rwd}';
                        temp = repmat(correctTrials{iSub,rwd},trialLength,1);
                        correctTRs{iSub,rwd} = temp(:);
                        
                        %this is for all data
                        stimTrials{iSub,rwd} = 1-nullTrials{goodSubs(iSub),rwd};
                        freqTrials{goodSubs(iSub),rwd} = freqTrials{goodSubs(iSub),rwd}.*stimTrials{iSub,rwd}.*correctTrials{iSub,rwd};%set null trials to 0
                        contrastTrials{goodSubs(iSub),rwd} = contrastTrials{goodSubs(iSub),rwd}.*stimTrials{iSub,rwd}.*correctTrials{iSub,rwd};%set null trials to 0
                        
                        %for pupil data
                        nonJunkedStimTrials{iSub,rwd} = eyeStimTrials{goodSubs(iSub),rwd};
                        nonJunkedStimTrials{iSub,rwd}(1:17:end) = [];
                        nonJunkedNullTrials{iSub,rwd} = eyeNullTrials{goodSubs(iSub),rwd};
                        nonJunkedNullTrials{iSub,rwd}(1:17:end) = [];
                        
                    end
                    
                    %concatenate across rwd, separately for each subject
                    stimPupilConcat{iSub} = [stimPupilRwd{goodSubs(iSub),1}(correctTrials{iSub,1}(nonJunkedStimTrials{goodSubs(iSub),1}==1)==1,1:pupilTrialTime); ...
                        stimPupilRwd{goodSubs(iSub),2}(correctTrials{iSub,2}(nonJunkedStimTrials{goodSubs(iSub),2}==1)==1,1:pupilTrialTime)];
                    nullPupilConcat{iSub} = [nullPupilRwd{goodSubs(iSub),1}(correctTrials{iSub,1}(nonJunkedNullTrials{goodSubs(iSub),1}==1)==1,1:pupilTrialTime); ...
                        nullPupilRwd{goodSubs(iSub),2}(correctTrials{iSub,2}(nonJunkedNullTrials{goodSubs(iSub),2}==1)==1,1:pupilTrialTime)];
                    
                    %subtract mean null trial
                    pupilMeanNull(iSub,:) = squeeze(nanmean(nullPupilConcat{iSub},1));%mean across trials
                    subPupilMeanNull = repmat(pupilMeanNull(iSub,:),[size(stimPupilConcat{iSub},1),1]);%repeat X stim trials
                    if subtractMean
                        stimPupilConcat{iSub} = stimPupilConcat{iSub} - subPupilMeanNull;
                    end
                    %average across all null, stim, and both trials
                    for rwd=1:2
                        allPupilRwd = [nullPupilRwd{goodSubs(iSub),rwd}(correctTrials{iSub,rwd}(nonJunkedNullTrials{goodSubs(iSub),rwd}==1)==1,1:pupilTrialTime); ...
                            stimPupilRwd{goodSubs(iSub),rwd}(correctTrials{iSub,rwd}(nonJunkedStimTrials{goodSubs(iSub),rwd}==1)==1,1:pupilTrialTime)];
                       pupilMeanRwd(iSub, rwd,:) = squeeze(nanmean(allPupilRwd,1));
                       pupilMeanNullRwd(iSub, rwd,:) = squeeze(nanmean(nullPupilRwd{goodSubs(iSub),rwd}(correctTrials{iSub,rwd}(nonJunkedNullTrials{goodSubs(iSub),rwd}==1)==1,1:pupilTrialTime),1));%mean across trials
                       pupilMeanStimRwd(iSub, rwd,:) = squeeze(nanmean(stimPupilRwd{goodSubs(iSub),rwd}(correctTrials{iSub,rwd}(nonJunkedStimTrials{goodSubs(iSub),rwd}==1)==1,1:pupilTrialTime),1));%mean across trials
                    end
                    %pupil amplitude measures and stats - STIM
                    subStdStimPupil{iSub} = nanstd(stimPupilConcat{iSub},0,2);%std amplitude per trial
                    %     subBaseStimPupil{iSub} = nanmean(stimPupilConcat{iSub}(:,1:baselineT),2);%per trial
                    temp = fillmissing(stimPupilConcat{iSub},'linear',2);
                    subBaseStimPupil{iSub} = nanmean(temp(:,1:baselineT),2);%per trial
                    f=fft(fillmissing(stimPupilConcat{iSub},'linear',2),[],2);
                    subPhStimPupil{iSub} = angle(f(:,2));
                    subAmpStimPupil{iSub} = abs(f(:,2));
                    
                    %pupil amplitude measures and stats - NULL
                    subStdNullPupil{iSub} = nanstd(nullPupilConcat{iSub},0,2);%per trial
                    subBaseNullPupil{iSub} = nanmean(nullPupilConcat{iSub}(:,1:baselineT),2);%per trial
                    f=fft(fillmissing(nullPupilConcat{iSub},'linear',2),[],2);
                    subPhNullPupil{iSub} = angle(f(:,2));
                    subAmpNullPupil{iSub} = abs(f(:,2));
                    
                    subStimPupilStdMedian(iSub) = nanmedian(subStdStimPupil{iSub});
                    subStimPupilBaseMedian(iSub) = nanmedian(subBaseStimPupil{iSub});
                    
                    %compute mean pupil for stim trials
                    meanStimPupilTrial(iSub,:) = squeeze(nanmean(fillmissing(stimPupilConcat{iSub},'linear',2),1));
                    %regress mean trial
                    A = zscore(meanStimPupilTrial(iSub,:));
                    B = fillmissing(stimPupilConcat{iSub},'linear',2);
                    %     B = stimPupilConcat{iSub};
                    subBetasStimPupil{iSub} = A'\B';%A'x=B
                    
                    %pulse stats
                    for rwd=1:2
                        stimTrials{iSub,rwd} = stimTrials{iSub,rwd}.*correctTrials{iSub,rwd};
                        nullTrials{goodSubs(iSub),rwd} = nullTrials{goodSubs(iSub),rwd}.*correctTrials{iSub,rwd};
                    end
                    subStimPulse{iSub} = [rwdPulseTrials{goodSubs(iSub),1}(:,stimTrials{iSub,1}==1) rwdPulseTrials{goodSubs(iSub),2}(:,stimTrials{iSub,2}==1)];%concatenate pulse of null trials
                    subNullPulse{iSub} = [rwdPulseTrials{goodSubs(iSub),1}(:,nullTrials{goodSubs(iSub),1}==1) rwdPulseTrials{goodSubs(iSub),2}(:,nullTrials{goodSubs(iSub),2}==1)];%concatenate pulse of null trials
                    
                    
                    %subtract mean null trial
                    pulseMeanNull(iSub,:) = squeeze(nanmean(subNullPulse{iSub},2));%mean across trials
                    for rwd=1:2
                        pulseMeanRwd(iSub,rwd,:) =  squeeze(nanmean(rwdPulseTrials{goodSubs(iSub),rwd},2));%mean across trials
                       pulseMeanNullRwd(iSub,rwd,:) =  squeeze(nanmean(rwdPulseTrials{goodSubs(iSub),rwd}(:,nullTrials{goodSubs(iSub),rwd}==1),2));%mean across trials
                       pulseMeanStimRwd(iSub,rwd,:) =  squeeze(nanmean(rwdPulseTrials{goodSubs(iSub),rwd}(:,nullTrials{goodSubs(iSub),rwd}==0),2));%mean across trials
                    end
                    subPulseMeanNull = repmat(pulseMeanNull(iSub,:),[size(subStimPulse{iSub},2),1]);%repeat X stim trials
                    if subtractMean
                        subStimPulse{iSub} = subStimPulse{iSub} - subPulseMeanNull';
                    end
                    
                    subStdStimPulse{iSub} = std(fillmissing(subStimPulse{iSub},'linear',2));
                    subStimPulseStdMedian(iSub) = nanmedian(subStdStimPulse{iSub});
                    
%                     subStimRwd{iSub} = 1:length(subStdStimPulse{iSub});
                    %     subStimRwdMedian(iSub) = length(freqTrials{goodSubs(iSub),1})+0.5;%median(subRwd{iSub});
                    subStimRwdMedian(iSub) = sum(stimTrials{iSub,1})+0.5;
                    %     subStimRwdMedian(iSub) = sum(nullTrialsTRs{goodSubs(iSub),1}==0)/trialLength+0.5;%after the last trial of high reward
                    
                    %     subStimPulseStd{iSub} = std(subStimPulse{iSub});
                    subMeanStimPulse{iSub} = nanmean(subStimPulse{iSub});
                    
                    %compute mean pulse for stim trials
                    meanStimPulseTrial(iSub,:) = squeeze(nanmean(fillmissing(subStimPulse{iSub},'linear',2),2));
                    %regress mean trial
                    A = zscore(meanStimPulseTrial(iSub,:));
                    B = fillmissing(subStimPulse{iSub},'linear',2);
                    %     B = stimPupilConcat{iSub};
                    subBetasStimPulse{iSub} = A'\B;%A'x=B
                    
                    subNullPulseStd{iSub} = std(subNullPulse{iSub});
                    subNullPulseMean{iSub} = mean(subNullPulse{iSub});
                    %             f=fft(fillmissing(subNullPulse{iSub},'linear',2),[],2);
                    
                    
                    subNullResp{iSub} = [rwdRvTrials{goodSubs(iSub),1}(:,nullTrials{goodSubs(iSub),1}==1) rwdRvTrials{goodSubs(iSub),2}(:,nullTrials{goodSubs(iSub),2}==1)];%concatenate pulse of null trials
                    subNullRespStd{iSub} = std(subNullResp{iSub});
                    subNullRespMean{iSub} = mean(subNullResp{iSub});
                    
                    
                    %get mean pupil for high and low contrast, and for frequencies
                    
                    stimContrastTrials{iSub} = [contrastTrials{goodSubs(iSub),1}(stimTrials{iSub,1}==1) contrastTrials{goodSubs(iSub),2}(stimTrials{iSub,2}==1)];
                    rwdTrials{iSub} = [ones(size(contrastTrials{goodSubs(iSub),1}(stimTrials{iSub,1}==1))) 2*ones(size(contrastTrials{goodSubs(iSub),2}(stimTrials{iSub,2}==1)))];
                    
                    for icontrast=1:ncontrasts
                        contrastPupil{iSub,icontrast} = stimPupilConcat{iSub}(stimContrastTrials{iSub}==icontrast,:);
                        meanContrastPupil(iSub,icontrast,:) = nanmean(contrastPupil{iSub,icontrast});
                        contrastPulse{iSub,icontrast} = subStimPulse{iSub}(:,stimContrastTrials{iSub}==icontrast);
                        meanContrastPulse(iSub,icontrast,:) = nanmean(contrastPulse{iSub,icontrast},2);
                    end
                    
                    stimFreqTrials{iSub} = [freqTrials{goodSubs(iSub),1}(stimTrials{iSub,1}==1) freqTrials{goodSubs(iSub),2}(stimTrials{iSub,2}==1)];
                    for ifreq=1:nfreqs
                        freqPupil{iSub,ifreq} = stimPupilConcat{iSub}(stimFreqTrials{iSub}==ifreq,:);
                        meanFreqPupil(iSub,ifreq,:) = nanmean(freqPupil{iSub,ifreq});
                        freqPulse{iSub,ifreq} = subStimPulse{iSub}(:,stimFreqTrials{iSub}==ifreq);
                        meanFreqPulse(iSub,ifreq,:) = nanmean(freqPulse{iSub,ifreq},2);
                    end
                    
                end
                
                %%
                
                
                controlROIs = 4;%left and right DMN
                for iSub=1:length(goodSubs)
                    %         for iRoi=1:length(controlROIs)
                    %             roiTseriesL = [roiTC{goodSubs(iSub),controlROIs(iRoi),1}.tSeries roiTC{goodSubs(iSub),controlROIs(iRoi),2}.tSeries];
                    %             roiTseriesR = [roiTC{goodSubs(iSub),controlROIs(iRoi)+1,1}.tSeries roiTC{goodSubs(iSub),controlROIs(iRoi)+1,2}.tSeries];
                    %             roiTseries = [roiTseriesL; roiTseriesR];
                    %
                    %             for rwd=1:2
                    %                 controlTseriesRwdL = roiTC{goodSubs(iSub),controlROIs(iRoi),rwd}.tSeries;
                    %                 controlTseriesRwdR = roiTC{goodSubs(iSub),controlROIs(iRoi),rwd}.tSeries;
                    %                 controlTseriesRwd = [controlTseriesRwdL; controlTseriesRwdR];
                    %                 controlNullTseriesRwd{iSub,iRoi,rwd} = controlTseriesRwd(:,nullTrialsTRs{goodSubs(iSub),rwd}==1 & correctTRs{iSub,rwd});
                    %                 controlStimTseriesRwd{iSub,iRoi,rwd} = controlTseriesRwd(:,nullTrialsTRs{goodSubs(iSub),rwd}==0 & correctTRs{iSub,rwd});
                    %             end
                    %
                    %             %concatenate across rwd, average across voxels
                    %             controlNullTseries{iSub,iRoi} = [ nanmean(controlNullTseriesRwd{iSub,iRoi,1},1) nanmean(controlNullTseriesRwd{iSub,iRoi,2},1)];
                    %             controlStimTseries{iSub,iRoi} = [ nanmean(controlStimTseriesRwd{iSub,iRoi,1},1) nanmean(controlStimTseriesRwd{iSub,iRoi,2},1)];
                    %             %reshape into trials
                    %             controlNullTrialTseries{iSub,iRoi} = reshape(controlNullTseries{iSub,iRoi},trialLength,[]);%T, trials
                    %             controlStimTrialTseries{iSub,iRoi} = reshape(controlStimTseries{iSub,iRoi},trialLength,[]);%T, trials
                    %
                    %             %subtract mean null trial
                    %             controlMeanNull{iSub,iRoi} = squeeze(mean(controlNullTrialTseries{iSub,iRoi},2));%mean across trials
                    %             subControlMeanNull = repmat(controlMeanNull{iSub,iRoi},[1,size(controlStimTrialTseries{iSub,iRoi},2)]);%repeat X stim trials
                    %             if subtractMean
                    %                 controlStimTrialTseries{iSub,iRoi} = controlStimTrialTseries{iSub,iRoi} - subControlMeanNull;
                    %             end
                    %
                    %             controlStimStd{iSub,iRoi} = squeeze(std(controlStimTrialTseries{iSub,iRoi},0,1));
                    %             f = fft(controlStimTrialTseries{iSub,iRoi});
                    %             controlStimAmp{iSub,iRoi} = squeeze(abs(f(2,:)));
                    %             controlMeanStimTrial{iSub,iRoi} = squeeze(nanmean(controlStimTrialTseries{iSub,iRoi},2));%T
                    %             ntrials = size(controlStimTrialTseries{iSub,iRoi},2);
                    %             %         for ivox=1:sum(goodVoxels)
                    %             A = controlMeanStimTrial{iSub,iRoi};
                    %             A = zscore(A);
                    %             B = controlStimTrialTseries{iSub,iRoi};
                    %             controlStimBetas{iSub,iRoi} = A\B;%A'x=B
                    %             %         end
                    %             %         subStimControlStdMedian{iSub,iRoi} = nanmedian(controlStimStd{iSub,iRoi});
                    %
                    %             controlStimBetasMedian(iSub) = median(controlStimBetas{iSub,1});
                    %         end
                    for iRoi=bensonROIs%1:length(roiNames) %length(bensonROIs)
                        for ibin=1:nbins
                            roiTseries = [roiTC{goodSubs(iSub),iRoi,1}.tSeries roiTC{goodSubs(iSub),iRoi,2}.tSeries];
                            goodVoxels = ~isnan(mean(roiTseries,2));
                            
                            
                            binVoxels = eccen{goodSubs(iSub),iRoi}>binBorders(ibin) & eccen{goodSubs(iSub),iRoi}<=binBorders(ibin+1);
                            %                                 binVoxels = binVoxels & areas{goodSubs(iSub),iRoi}==1;%ONLY V1
                            binVoxels = goodVoxels & binVoxels & areas{goodSubs(iSub),iRoi}>0;%ONLY V1,V2,V3
                            numBinVoxels(iSub,iRoi,ibin) = sum(binVoxels);
                            
                            for rwd=1:2
                                binTseriesRwd = roiTC{goodSubs(iSub),iRoi,rwd}.tSeries(binVoxels,:);
                                binNullTseriesRwd{iSub,iRoi,ibin,rwd} = binTseriesRwd(:,nullTrialsTRs{goodSubs(iSub),rwd}==1 & correctTRs{iSub,rwd});
                                binStimTseriesRwd{iSub,iRoi,ibin,rwd} = binTseriesRwd(:,nullTrialsTRs{goodSubs(iSub),rwd}==0 & correctTRs{iSub,rwd});
                            end
                            
                            %concatenate across rwd
                            binNullTseries{iSub,iRoi,ibin} = [ binNullTseriesRwd{iSub,iRoi,ibin,1} binNullTseriesRwd{iSub,iRoi,ibin,2}];
                            binStimTseries{iSub,iRoi,ibin} = [ binStimTseriesRwd{iSub,iRoi,ibin,1} binStimTseriesRwd{iSub,iRoi,ibin,2}];
                            %reshape into trials
                            binNullTrialTseries{iSub,iRoi,ibin} = reshape(binNullTseries{iSub,iRoi,ibin},numBinVoxels(iSub,iRoi,ibin),trialLength,[]);%vox, T, trials
                            binStimTrialTseries{iSub,iRoi,ibin} = reshape(binStimTseries{iSub,iRoi,ibin},numBinVoxels(iSub,iRoi,ibin),trialLength,[]);%vox, T, trials
                            
                            %subtract mean null trial
                            voxMeanNull{iSub,iRoi,ibin} = squeeze(mean(binNullTrialTseries{iSub,iRoi,ibin},3));%mean across trials
                            subBinMeanNull = repmat(voxMeanNull{iSub,iRoi,ibin},[1,1,size(binStimTrialTseries{iSub,iRoi,ibin},3)]);%repeat X stim trials
                            if subtractMean
                                binStimTrialTseries{iSub,iRoi,ibin} = binStimTrialTseries{iSub,iRoi,ibin} - subBinMeanNull;
                            end
                            
                            %get fMRI response amplitude
                            %compute mean stim trial per voxel
                            binMeanStimTrialVox{iSub,iRoi,ibin} = squeeze(nanmean(binStimTrialTseries{iSub,iRoi,ibin},3));%vox, T
                            binMeanNullTrialVox{iSub,iRoi,ibin} = squeeze(nanmean(binNullTrialTseries{iSub,iRoi,ibin},3));%vox, T
                            %regress mean trial
                            ntrials = size(binStimTrialTseries{iSub,iRoi,ibin},3);
                            subBinBetasStimVox{iSub,iRoi,ibin}=[];
                            for ivox=1:numBinVoxels(iSub,iRoi,ibin)
                                
                                A = zscore(binMeanStimTrialVox{iSub,iRoi,ibin}(ivox,:));
                                B = squeeze(binStimTrialTseries{iSub,iRoi,ibin}(ivox,:,:));
                                subBinBetasStimVox{iSub,iRoi,ibin}(ivox,:) = A'\B;%A'x=B
                                %                 Y = squeeze(binStimTrialTseries{iSub,iRoi,ibin}(ivox,:,:));
                                %                 b = regress(Y,X);%Y=X*b
                            end
                            %compute mean trial across entire bin
                            binMeanStimTrial(iSub,iRoi,ibin,:) = mean(binMeanStimTrialVox{iSub,iRoi,ibin},1);%mean stim trial
                            binMeanNullTrial(iSub,iRoi,ibin,:) = mean(binMeanNullTrialVox{iSub,iRoi,ibin},1);%mean null trial
                            A = squeeze(zscore(binMeanStimTrial(iSub,iRoi,ibin,:)));
                            B = squeeze(mean(binStimTrialTseries{iSub,iRoi,ibin},1));%mean across voxels, per trial
                            subBinBetasStimTrial{iSub,iRoi,ibin} = A\B;%one beta per trial per bin
                            
                            %             stimFreqTrials{iSub} = [freqTrials{goodSubs(iSub),1}(stimTrials{goodSubs(iSub),1}==1) freqTrials{goodSubs(iSub),2}(stimTrials{goodSubs(iSub),2}==1)];
                            for rwd=1:2
                                for ifreq=1:max(freqTrials{goodSubs(iSub),rwd})
                                    meanBinFreqTrial = squeeze(mean(binStimTrialTseries{iSub,iRoi,ibin}(:,:,stimFreqTrials{iSub}==ifreq & rwdTrials{iSub}==rwd),1));%mean over voxels. TxnumTrials
                                    binMeanFreqTrial(iSub,iRoi,ibin,rwd,ifreq,:) = mean(meanBinFreqTrial,2);
                                    B = squeeze(binMeanFreqTrial(iSub,iRoi,ibin,rwd,ifreq,:));
                                    subBinBetasFreq(iSub,iRoi,ibin,rwd,ifreq) = A\B;
                                end
                                
                                for icontrast=1:max(contrastTrials{goodSubs(iSub),rwd})
                                    meanBinContrastTrial = squeeze(mean(binStimTrialTseries{iSub,iRoi,ibin}(:,:,stimContrastTrials{iSub}==icontrast & rwdTrials{iSub}==rwd),1));%mean over voxels. TxnumTrials
                                    binMeanContrastTrial(iSub,iRoi,ibin,rwd,icontrast,:) = mean(meanBinContrastTrial,2);
                                    B = squeeze(binMeanContrastTrial(iSub,iRoi,ibin,rwd,icontrast,:));
                                    subBinBetasContrast(iSub,iRoi,ibin,rwd,icontrast) = A\B;
                                end
                            end
                        end
                    end
                end
                
                histEdgesFreq = [1:nfreqs+1]-0.5;
                histEdgesContrast = [1:ncontrasts+1]-0.5;
                
                
                
                %     %arousal
                %     clear classAccFreq classAccContrast permAccFreq permAccContrast
                %     for arousalType=arousalTypes
                %
                %         for iSub=1:length(goodSubs)
                %
                %             switch arousalType
                %                 case 1
                %                     trialArousal = subBaseStimPupil{iSub}';
                %                     medianArousal = subStimPupilBaseMedian(iSub);
                %                 case 2
                %                     trialArousal = subStdStimPupil{iSub}';
                %                     medianArousal = subStimPupilStdMedian(iSub);
                %                 case 3
                %                     trialArousal = subStdStimPulse{iSub}';
                %                     medianArousal = subStimPulseStdMedian(iSub);
                %                 case 4
                %                     trialArousal = -subStimRwd{iSub}';
                %                     medianArousal = -subStimRwdMedian(iSub);
                %                 case 5
                %                     trialArousal = controlStimBetas{iSub,1}';
                %                     medianArousal = controlStimBetasMedian(iSub);
                %             end
                %
                %             for arousal=1:2
                %                 if arousal==1
                %                     arousalTrials{arousalType,iSub,arousal} = trialArousal>medianArousal;
                %                 else
                %                     arousalTrials{arousalType,iSub,arousal} = trialArousal<medianArousal;
                %                 end
                %                 [iSub, arousal, sum(arousalTrials{arousalType,iSub,arousal})];
                %
                %             end
                %             for arousal=1:2
                %                 permArousalTrials{iSub,arousal}=[];
                %                 for iperm=1:nperms
                %                     permArousalTrials{iSub,arousal}(iperm,:) = randperm(sum(arousalTrials{arousalType,iSub,arousal}));
                %                     %                 permArousalTrials{iSub,arousal}(iperm,:) = find(arousalTrials{iSub,arousal}(randperm(length(arousalTrials{iSub,arousal}))));
                %                 end
                %             end
                %
                %         end
                %     end
                toc
                
                %%
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
                
                
                %%
                save([dataFolder 'allData' concatProjStr onlyCorrectStr zscoreStr subtractMeanStr '.mat'],...
                    'pupilTrialTime','goodSubs', 'subFolders',...
                    'onlyCorrect','ncontrasts','nfreqs',...
                    'subtractMean',...
                    'toZscore','concatProj','bensonROIs','roiNames','eccMin','eccMax','nbins',...
                    'binBorders','nbins',...
                    'meanContrastPupil','meanContrastPulse','meanFreqPupil','meanFreqPulse',...
                    'stimPupilConcat','subBaseStimPupil','subStimPulse',...
                    'subStdStimPupil','subStdStimPulse','stimTrials',...
                    'stimTrials','subStimPupilBaseMedian','subStimPupilStdMedian', 'subStimPulseStdMedian', 'subStimRwdMedian',...
                    'subBetasStimPupil','subBetasStimPulse',...
                    'stimFreqTrials','stimContrastTrials',...
                    'subMeanStimPulse',...
                    'subBinBetasStimVox', 'binMeanStimTrial','binMeanStimTrialVox','subBinBetasStimTrial',...
                    'subBinBetasFreq','subBinBetasContrast','rwdTrials',...
                    'pupilMeanNull', 'pulseMeanNull','voxMeanNull','binMeanNullTrialVox','binMeanNullTrial',...
                    'pulseMeanNullRwd','pulseMeanStimRwd','pulseMeanRwd','pupilMeanNullRwd','pupilMeanStimRwd','pupilMeanRwd');
                
            end
        end
    end
end
close all
mrQuit
clear all
trialsPerRun=16;
% runLength = 160;
includePulse = 1;
subtractMeanNull = 1;%slightly improves max accuracy for 16,32 trials
includePupil = 1;
includePupilBase = 1;
featureTypes=1:3;
% featureType = 1;%1-std, 2-amp, 3-regression.
toZscore=1;%0 or 1
concatProj= 1;

bensonROIs = 1;%[1:3];
eccMin = 0.2;
eccMax = 70;
nbins = 12;
binBorders = logspace(log10(eccMin),log10(eccMax),nbins+1);
nbins = length(binBorders)-1;
for i=2:length(binBorders)
    binCenters(i-1) = (binBorders(i)+binBorders(i-1))/2;
end
binBorders;

saveFigs=1;
nfreqs=5;
ncontrasts=2;
tic

% dataFolder = '/Volumes/MH02086153MACDT-Drobo/rwdFmri/';
dataFolder = '/Users/rothzn/rwdFmri/';
freqs = logspace(-0.3,0.5,5);
contrasts = logspace(-0.7,0,2);
%Load eyetracking data
load([dataFolder 'rwdRapidEyeData.mat'], 'subFolders', 'samplerate',  ...
    'trialLength', ...
    'numSubs', 'onlyCorrect', ...
    'sacRwd','binnedSac', 'smoothSac', ...
    'rwdPupil','meanPupil','taskTimes',...
    'sacRun','pupilRun', 'numRuns','startTimes','endTimes','eyetrackerTime',...
    'nullTrials','nullPupil','stimPupil','meanNullPupil','meanStimPupil',...
    'sacRunSmooth', 'sacRunSmoothTrialZeroFilled','sacRunSmoothTrial','L');




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

baselineT = 25;


%Load physiology
ConcatProjStr = '';
if concatProj
    ConcatProjStr = 'concatProj';
end
load([dataFolder 'rwdRapid_physio.mat'], 'concatInfo',  ...
    'subFolders','trialLength',...
    'ecgselect','ecgSampleRate','ecgTrial','ecgRunLength','ecgInterpMethod',...
    'ecg','ecgPulseRate','interpPulseRate',...
    'respselect','resp',...
    'rwdPulseTrials','rwdRvTrials','meanPulse','meanRV',...
    'designMatPulse','designMatResp','designMatRespPulse');


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
    'frames', 'junkedFrames', 'TR', 'trialLength', 'nVolumes',...
    'eccen','ang','areas',...
    'numVox','roiTC','nullTrialsTRs',...
    'nullTrialsRun','nullTrials','contrastTrialsRun','freqTrialsRun','contrastTrials','freqTrials',...
    'nullTseries','nullTrialTseries','stimTseries','stimTrialTseries');



for iSub=1:size(nullTrials,1)
    for rwd=1:2
        trialRunRwd{iSub,rwd} = ceil([1:size(nullTrials{iSub,rwd},2)]/trialsPerRun);
        stimTrials{iSub,rwd} = 1-nullTrials{iSub,rwd};
        freqTrials{iSub,rwd} = freqTrials{iSub,rwd}.*stimTrials{iSub,rwd};%set null trials to 0
        contrastTrials{iSub,rwd} = contrastTrials{iSub,rwd}.*stimTrials{iSub,rwd};%set null trials to 0
    end
    trialRun{iSub} = [trialRunRwd{iSub,1} trialRunRwd{iSub,2}];
end
goodSubs = [1:12 14];
goodSubs = [1 4:5 7:12 14];
% goodSubs = 1:2;

%% get pupil data
nullPupilRwd = origNullPupil;
stimPupilRwd = origStimPupil;
for iSub=1:length(goodSubs)
    
    for rwd=1:2
        %remove the first trial of each run
        junkedIndices = zeros(size(eyeStimTrials{goodSubs(iSub),rwd}));
        junkedIndices(1:17:end) = ones;
        stimTrialIndices = find(eyeStimTrials{goodSubs(iSub),rwd}==1);
        junkedStimTrials = find(mod(stimTrialIndices,17)==1);
        stimPupilRwd{goodSubs(iSub),rwd}(junkedStimTrials,:) = [];
        
        junkedIndices = zeros(size(eyeNullTrials{goodSubs(iSub),rwd}));
        junkedIndices(1:17:end) = ones;
        nullTrialIndices = find(eyeNullTrials{goodSubs(iSub),rwd}==1);
        junkedNullTrials = find(mod(nullTrialIndices,17)==1);
        nullPupilRwd{goodSubs(iSub),rwd}(junkedNullTrials,:) = [];
        
        
    end
    
    %concatenate across rwd, separately for each subject
    stimPupilConcat{iSub} = [stimPupilRwd{goodSubs(iSub),1}(:,1:pupilTrialTime); stimPupilRwd{goodSubs(iSub),2}(:,1:pupilTrialTime)];
    nullPupilConcat{iSub} = [nullPupilRwd{goodSubs(iSub),1}(:,1:pupilTrialTime); nullPupilRwd{goodSubs(iSub),2}(:,1:pupilTrialTime)];
    
    %get run number for each stim trial
    trialRunStim{iSub} = [trialRunRwd{goodSubs(iSub),1}(stimTrials{goodSubs(iSub),1}==1) max(trialRunRwd{goodSubs(iSub),1})+trialRunRwd{goodSubs(iSub),2}(stimTrials{goodSubs(iSub),2}==1)];
    
    %get condition per trial
    stimContrastTrials{iSub} = [contrastTrials{goodSubs(iSub),1}(stimTrials{goodSubs(iSub),1}==1) contrastTrials{goodSubs(iSub),2}(stimTrials{goodSubs(iSub),2}==1)];
    stimFreqTrials{iSub} = [freqTrials{goodSubs(iSub),1}(stimTrials{goodSubs(iSub),1}==1) freqTrials{goodSubs(iSub),2}(stimTrials{goodSubs(iSub),2}==1)];
    
    
    %pulse data
    subStimPulse{iSub} = [rwdPulseTrials{goodSubs(iSub),1}(:,stimTrials{goodSubs(iSub),1}==1) rwdPulseTrials{goodSubs(iSub),2}(:,stimTrials{goodSubs(iSub),2}==1)];%concatenate pulse of null trials
    subNullPulse{iSub} = [rwdPulseTrials{goodSubs(iSub),1}(:,nullTrials{goodSubs(iSub),1}==1) rwdPulseTrials{goodSubs(iSub),2}(:,nullTrials{goodSubs(iSub),2}==1)];%concatenate pulse of null trials
    
    %subtract mean null trial
    pulseMeanNull(iSub,:) = squeeze(nanmean(subNullPulse{iSub},2));%mean across trials
    subPulseMeanNull = repmat(pulseMeanNull(iSub,:),[size(subStimPulse{iSub},2),1]);%repeat X stim trials
    if subtractMeanNull
        subStimPulse{iSub} = subStimPulse{iSub} - subPulseMeanNull';
    end
end

%% get fMRI data
for iSub=1:length(goodSubs)
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
                binNullTseriesRwd{iSub,iRoi,ibin,rwd} = binTseriesRwd(:,nullTrialsTRs{goodSubs(iSub),rwd}==1);
                binStimTseriesRwd{iSub,iRoi,ibin,rwd} = binTseriesRwd(:,nullTrialsTRs{goodSubs(iSub),rwd}==0);
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
            if subtractMeanNull
                binStimTrialTseries{iSub,iRoi,ibin} = binStimTrialTseries{iSub,iRoi,ibin} - subBinMeanNull;
            end
        end
    end
end


%% RUN AVERAGE
for iSub=1:length(goodSubs)
    %average BOLD, pupil and pulse for each freq/contrast, within each run
    oldStimPupilConcat = stimPupilConcat{iSub};
    stimPupilConcat{iSub} = [];
    oldStimPulse= subStimPulse{iSub};
    subStimPulse{iSub} = [];
    nruns=max(trialRunStim{iSub});
%     nfreqs=max(stimFreqTrials{iSub});
    for irun=1:nruns
        for icontrast=1:max(stimContrastTrials{iSub})
%             for  ifreq=1:nfreqs
%                 ind = ifreq + (irun-1)*nfreqs;
                ind = icontrast + (irun-1)*ncontrasts;
                trialsInd = trialRunStim{iSub}==irun & stimContrastTrials{iSub}==icontrast;
                stimPupilConcat{iSub}(ind,:) = nanmean(oldStimPupilConcat(trialsInd,:),1);
                subStimPulse{iSub}(:,ind) = nanmean(oldStimPulse(:,trialsInd),2);
%                 stimFreqRuns{iSub}(ind) = ifreq;
                stimContrastRuns{iSub}(ind) = icontrast;
            end
%         end
    end
    
    
    for iRoi=bensonROIs%1:length(roiNames) %length(bensonROIs)
        for ibin=1:nbins
            oldStimTrialTseries = binStimTrialTseries{iSub,iRoi,ibin};
            binStimTrialTseries{iSub,iRoi,ibin}= [];
            for irun=1:nruns
                for icontrast=1:ncontrasts
%                     for  ifreq=1:nfreqs
%                         ind = ifreq + (irun-1)*nfreqs;
                        ind = icontrast + (irun-1)*ncontrasts;
%                         trialsInd = trialRunStim{iSub}==irun & stimFreqTrials{iSub}==ifreq;
                        trialsInd = trialRunStim{iSub}==irun & stimContrastTrials{iSub}==icontrast;
                        binStimTrialTseries{iSub,iRoi,ibin}(:,:,ind) = nanmean(oldStimTrialTseries(:,:,trialsInd),3);
%                     end
                end
            end
        end
    end
end




%% get physio statistics
for iSub=1:length(goodSubs)
    
    
    %subtract mean null trial
    pupilMeanNull(iSub,:) = squeeze(nanmean(nullPupilConcat{iSub},1));%mean across trials
    subPupilMeanNull = repmat(pupilMeanNull(iSub,:),[size(stimPupilConcat{iSub},1),1]);%repeat X stim trials
    if subtractMeanNull
        stimPupilConcat{iSub} = stimPupilConcat{iSub} - subPupilMeanNull;
    end
    
    %pupil amplitude measures and stats - STIM
    subStdStimPupil{iSub} = nanstd(stimPupilConcat{iSub},0,2);%per trial
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
    A = meanStimPupilTrial(iSub,:);
    B = fillmissing(stimPupilConcat{iSub},'linear',2);
    %     B = stimPupilConcat{iSub};
    subBetasStimPupil{iSub} = A'\B';%A'x=B
    
    
    %PULSE
    
    subStimPulseStd{iSub} = std(fillmissing(subStimPulse{iSub},'linear',2));
    %     subStimPulseStd{iSub} = std(subStimPulse{iSub});
    subStimPulseMean{iSub} = mean(subStimPulse{iSub});
    f=fft(fillmissing(subStimPulse{iSub},'linear',2),[],2);
    %     f=fft(subStimPulse{iSub});
    subStimPulsePh{iSub} = angle(f(2,:));
    subStimPulseAmp{iSub} = abs(f(2,:));
    %compute mean pulse for stim trials
    meanStimPulseTrial(iSub,:) = squeeze(nanmean(fillmissing(subStimPulse{iSub},'linear',2),2));
    %regress mean trial
    A = meanStimPulseTrial(iSub,:);
    B = fillmissing(subStimPulse{iSub},'linear',2);
    %     B = stimPupilConcat{iSub};
    subBetasStimPulse{iSub} = A'\B;%A'x=B
    
    subNullPulseStd{iSub} = std(subNullPulse{iSub});
    subNullPulseMean{iSub} = mean(subNullPulse{iSub});
    %             f=fft(fillmissing(subNullPulse{iSub},'linear',2),[],2);
    f=fft(subNullPulse{iSub});
    subNullPulsePh{iSub} = angle(f(2,:));
    subNullPulseAmp{iSub} = abs(f(2,:));
    
    
    subNullResp{iSub} = [rwdRvTrials{goodSubs(iSub),1}(:,nullTrials{goodSubs(iSub),1}==1) rwdRvTrials{goodSubs(iSub),2}(:,nullTrials{goodSubs(iSub),2}==1)];%concatenate pulse of null trials
    subNullRespStd{iSub} = std(subNullResp{iSub});
    subNullRespMean{iSub} = mean(subNullResp{iSub});
    f=fft(subNullResp{iSub});
    subNullRespPh{iSub} = angle(f(2,:));
    subNullRespAmp{iSub} = abs(f(2,:));
    
    
    %get mean pupil for high and low contrast, and for frequencies
    for icontrast=1:max(stimContrastRuns{iSub})
        contrastPupil{iSub,icontrast} = stimPupilConcat{iSub}(stimContrastRuns{iSub}==icontrast,:);
        meanContrastPupil(iSub,icontrast,:) = nanmean(contrastPupil{iSub,icontrast});
    end
%     for ifreq=1:max(stimFreqRuns{iSub})
%         freqPupil{iSub,ifreq} = stimPupilConcat{iSub}(stimFreqRuns{iSub}==ifreq,:);
%         meanFreqPupil(iSub,ifreq,:) = nanmean(freqPupil{iSub,ifreq});
%     end
    
    
end




%% get fMRI statistics
for iSub=1:length(goodSubs)
    for iRoi=bensonROIs%1:length(roiNames) %length(bensonROIs)
        for ibin=1:nbins
            %get fMRI response amplitude
            %STD
            subBinStdNull{iSub,iRoi,ibin} = squeeze(std(binNullTrialTseries{iSub,iRoi,ibin},0,2));%vox, trial
            subBinStdStim{iSub,iRoi,ibin} = squeeze(std(binStimTrialTseries{iSub,iRoi,ibin},0,2));%vox, trial
            %FFT amplitude
            f = fft(binStimTrialTseries{iSub,iRoi,ibin},[],2);
            subBinAmpStim{iSub,iRoi,ibin} = squeeze(abs(f(:,2,:)));
            
            %compute mean stim trial per voxel
            binMeanStimTrial{iSub,iRoi,ibin} = squeeze(nanmean(binStimTrialTseries{iSub,iRoi,ibin},3));%vox, T
            %regress mean trial
            ntrials = size(binStimTrialTseries{iSub,iRoi,ibin},3);
            for ivox=1:numBinVoxels(iSub,iRoi,ibin)
                
                A = binMeanStimTrial{iSub,iRoi,ibin}(ivox,:);
                B = squeeze(binStimTrialTseries{iSub,iRoi,ibin}(ivox,:,:));
                subBinBetasStim{iSub,iRoi,ibin}(ivox,:) = A'\B;%A'x=B
                %                 Y = squeeze(binStimTrialTseries{iSub,iRoi,ibin}(ivox,:,:));
                %                 b = regress(Y,X);%Y=X*b
            end
            
            %             subBinBetasStim{iSub,iRoi,ibin} = x;
            
            %             stimFreqTrials{iSub} = [freqTrials{goodSubs(iSub),1}(stimTrials{goodSubs(iSub),1}==1) freqTrials{goodSubs(iSub),2}(stimTrials{goodSubs(iSub),2}==1)];
%             for ifreq=1:max(freqTrials{goodSubs(iSub),1})
%                 subBinFreqStd{iSub,iRoi,ibin,ifreq} =  subBinStdStim{iSub,iRoi,ibin}(:,stimFreqRuns{iSub}==ifreq,:);%vox,trial
%             end
            
            %             stimContrastTrials{iSub} =  [contrastTrials{goodSubs(iSub),1}(stimTrials{goodSubs(iSub),1}==1) contrastTrials{goodSubs(iSub),2}(stimTrials{goodSubs(iSub),2}==1)];
            for icontrast=1:max(contrastTrials{goodSubs(iSub),1})
                subBinContrastStd{iSub,iRoi,ibin,icontrast} =  subBinStdStim{iSub,iRoi,ibin}(:,stimContrastRuns{iSub}==icontrast,:);%vox,trial
            end
            %arousal
            for arousal=1:2
                if arousal==1
                    arousalTrials = subBaseStimPupil{iSub}'>subStimPupilBaseMedian(iSub);
                else
                    arousalTrials = subBaseStimPupil{iSub}'<subStimPupilBaseMedian(iSub);
                end
%                 for ifreq=1:max(freqTrials{goodSubs(iSub),1})
%                     subBinFreqStdArousal{iSub,iRoi,ibin,ifreq,arousal} =  subBinStdStim{iSub,iRoi,ibin}(:,stimFreqRuns{iSub}==ifreq & arousalTrials,:);%vox,trial
%                 end
                for icontrast=1:max(contrastTrials{goodSubs(iSub),1})
                    subBinContrastStdArousal{iSub,iRoi,ibin,icontrast,arousal} =  subBinStdStim{iSub,iRoi,ibin}(:,stimContrastRuns{iSub}==icontrast & arousalTrials,:);%vox,trial
                end
            end
        end
    end
end

ifig=0;
%% DECODING
for featureType=featureTypes
    for iSub=1:length(goodSubs)
        nruns=max(trialRunStim{iSub});
        for iRoi=bensonROIs%1:length(roiNames) %length(bensonROIs)
            for ibin=1:nbins
                
                %decode frequency and contrast
                switch featureType
                    case 1
                        trialDataFmri = subBinStdStim{iSub,iRoi,ibin};
                        trialDataPupil = subStdStimPupil{iSub}';
                        trialDataPulse = subStimPulseStd{iSub};
                    case 2
                        trialDataFmri = subBinAmpStim{iSub,iRoi,ibin};
                        trialDataPupil = subAmpStimPupil{iSub}';
                        trialDataPulse = subStimPulseAmp{iSub};
                    case 3
                        trialDataFmri = subBinBetasStim{iSub,iRoi,ibin};
                        trialDataPupil = subBetasStimPupil{iSub};
                        trialDataPulse = subBetasStimPulse{iSub};
                end
                
                
                trialData =  trialDataFmri;
                if includePupil
                    trialData = [trialData; trialDataPupil];
                end
                if includePupilBase
                    trialData = [trialData; subBaseStimPupil{iSub}'];
                end
                if includePulse
                    trialData = [trialData; trialDataPulse];
                end
                
                
                ntrials = size(trialData,2);
                %                 ind = ifreq + (icontrast-1)*nfreqs + (irun-1)*ncontrasts*nfreqs;
                
                for irun=1:nruns
%                     indFirst = (irun-1)*nfreqs + 1;
                    indFirst = (irun-1)*ncontrasts + 1;
                    indLast = irun*ncontrasts;
%                     indLast = irun*nfreqs;
                    trainData =  trialData(:,[1:indFirst-1 indLast+1:ntrials]);
%                     trainLabelsFreq = stimFreqRuns{iSub}([1:indFirst-1 indLast+1:ntrials]);
                    testData = trialData(:,indFirst:indLast);
%                     classGuessFreq{iSub,iRoi,ibin}(indFirst:indLast) = classify(testData', trainData', trainLabelsFreq, 'diagLinear');
                    
                    trainLabelsContrast = stimContrastRuns{iSub}([1:indFirst-1 indLast+1:ntrials]);
%                     %     testData = trialData(:,indFirst:indLast);
                    classGuessContrast{iSub,iRoi,ibin}(indFirst:indLast) = classify(testData', trainData', trainLabelsContrast, 'diagLinear');
                    
                    
                end
                %                 classConfMatFreq{iSub,iRoi,ibin} = confusionmat(stimFreqTrials{iSub}, classGuessFreq{iSub,iRoi,ibin});
%                 classAccFreq{featureType}(iSub,iRoi,ibin) = sum(stimFreqRuns{iSub}==classGuessFreq{iSub,iRoi,ibin})/ntrials;%make sure not to include NaN trials?
                
                %                 classConfMatContrast{iSub,iRoi,ibin} = confusionmat(stimContrastTrials{iSub}, classGuessContrast{iSub,iRoi,ibin});
                classAccContrast{featureType}(iSub,iRoi,ibin) = sum(stimContrastRuns{iSub}==classGuessContrast{iSub,iRoi,ibin})/ntrials;%make sure not to include NaN trials?
            end
        end
    end
end



%%
plotColors = { [1 0 0], [0 0 1], [0 1 0], [0.5 1 0.2], [1 0.2 0.5] };
ifig=ifig+1; figure(ifig); clf
rows=1;
cols=2;
subplot(rows,cols,1)
for featureType=featureTypes
%     plot(squeeze(mean(classAccFreq{featureType})),'color',plotColors{featureType})
    plot(squeeze(mean(classAccContrast{featureType})),'color',plotColors{featureType})
    hold on
    %     mean(squeeze(mean(classAccFreq{featureType})))
end
legend('std','amp','regress');
% for featureType=featureTypes
%     plot(squeeze(mean(classAccContrast{featureType})),'color',plotColors{featureType})
%     hold on
%     %     mean(squeeze(mean(classAccContrast{featureType})))
% end
xlabel('eccentricity bin');
ylabel('proportion correct');
% legend('frequency','contrast');
% hline(1/nfreqs);
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
if subtractMeanNull
    titleStr = [titleStr 'mean subtracted '];
end


title(titleStr);

subplot(rows,cols,2)
for featureType=featureTypes
%     plot(squeeze(mean(classAccFreq{featureType},3)),'color',plotColors{featureType})
        plot(squeeze(mean(classAccContrast{featureType},3)),'color',plotColors{featureType})
    hold on

end
xlabel('subject #');
% legend('frequency','contrast');
% hline(1/nfreqs);
hline(1/ncontrasts);

old = {' ', ','};
new = '_';
% savefig([dataFolder replace(titleStr,old, new)]);

% max(squeeze(mean(classAccFreq{featureType},1)))
max(squeeze(mean(classAccContrast{featureType},1)))
% savefig(replace(titleStr,old, new));
toc

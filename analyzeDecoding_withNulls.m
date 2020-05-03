close all
mrQuit
clear all
arousalType=3;

includePulse = 0;
subtractMean = 0;
includePupil = 0;
includePupilBase = 0;
featureTypes=3;
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
nfreqs=6;
ncontrasts=3;
tic

% dataFolder = '/Volumes/MH02086153MACDT-Drobo/rwdFmri/';
dataFolder = '/Users/rothzn/rwdFmri/';
freqs = logspace(-0.3,0.5,5);
contrasts = logspace(-0.7,0,2);
%Load eyetracking data
load([dataFolder 'rwdRapidEyeData.mat'], 'subFolders', 'samplerate',  ...
    'trialsPerRun', 'trialLength', ...
    'numSubs', 'onlyCorrect', ...
    'sacRwd','binnedSac', 'smoothSac', ...
    'rwdPupil','meanPupil','taskTimes',...
    'sacRun','pupilRun', 'numRuns','startTimes','endTimes','eyetrackerTime',...
    'nullTrials','nullPupil','stimPupil','meanNullPupil','meanStimPupil',...
    'sacRunSmooth', 'sacRunSmoothTrialZeroFilled','sacRunSmoothTrial','L');
origPupil = rwdPupil;
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
    'frames', 'junkedFrames', 'TR', 'trialsPerRun', 'trialLength', 'nVolumes',...
    'eccen','ang','areas',...
    'numVox','roiTC','nullTrialsTRs',...
    'nullTrialsRun','nullTrials','contrastTrialsRun','freqTrialsRun','contrastTrials','freqTrials',...
    'nullTseries','nullTrialTseries','stimTseries','stimTrialTseries');
% for iSub=1:size(nullTrials,1)
%     for rwd=1:2
%         stimTrials{iSub,rwd} = 1-nullTrials{iSub,rwd};
%         freqTrials{iSub,rwd} = freqTrials{iSub,rwd}.*stimTrials{iSub,rwd};%set null trials to 0
%         contrastTrials{iSub,rwd} = contrastTrials{iSub,rwd}.*stimTrials{iSub,rwd};%set null trials to 0
%     end
% end
goodSubs = [1:12 14];
% goodSubs = 1:2;
%%
pupilRwd = origPupil;
% nullPupilRwd = origNullPupil;
% stimPupilRwd = origStimPupil;
for iSub=1:length(goodSubs)
    
    for rwd=1:2
        %remove the first trial of each run
        junkedIndices = zeros(size(eyeStimTrials{goodSubs(iSub),rwd}));
        junkedIndices(1:17:end) = ones;
        pupilRwd{goodSubs(iSub),rwd}(junkedIndices>0,:) = [];
        
%         junkedIndices = zeros(size(eyeStimTrials{goodSubs(iSub),rwd}));
%         junkedIndices(1:17:end) = ones;
%         stimTrialIndices = find(eyeStimTrials{goodSubs(iSub),rwd}==1);
%         junkedStimTrials = find(mod(stimTrialIndices,17)==1);
%         stimPupilRwd{goodSubs(iSub),rwd}(junkedStimTrials,:) = [];
        
%         junkedIndices = zeros(size(eyeNullTrials{goodSubs(iSub),rwd}));
%         junkedIndices(1:17:end) = ones;
%         nullTrialIndices = find(eyeNullTrials{goodSubs(iSub),rwd}==1);
%         junkedNullTrials = find(mod(nullTrialIndices,17)==1);
%         nullPupilRwd{goodSubs(iSub),rwd}(junkedNullTrials,:) = [];
        
        %get measure of pupil size per trial
%         stdStimPupil{iSub,rwd} = nanstd(stimPupilRwd{goodSubs(iSub),rwd},0,2);%per trial
%         baseStimPupil{iSub,rwd} = nanmean(stimPupilRwd{goodSubs(iSub),rwd}(:,1:baselineT),2);%per trial
%         f=fft(fillmissing(stimPupilRwd{goodSubs(iSub),rwd},'linear',2),[],2);
%         phStimPupil{iSub,rwd} = angle(f(:,2));
%         ampStimPupil{iSub,rwd} = abs(f(:,2));
%         
%         stdNullPupil{iSub,rwd} = nanstd(nullPupilRwd{goodSubs(iSub),rwd},0,2);%per trial
%         baseNullPupil{iSub,rwd} = nanmean(nullPupilRwd{goodSubs(iSub),rwd}(:,1:baselineT),2);%per trial
%         f=fft(fillmissing(nullPupilRwd{goodSubs(iSub),rwd},'linear',2),[],2);
%         phNullPupil{iSub,rwd} = angle(f(:,2));
%         ampNullPupil{iSub,rwd} = abs(f(:,2));
        
    end
    
    %concatenate across rwd, separately for each subject
    
    pupilConcat{iSub} = [pupilRwd{goodSubs(iSub),1}(:,1:pupilTrialTime); pupilRwd{goodSubs(iSub),2}(:,1:pupilTrialTime)];
%     stimPupilConcat{iSub} = [stimPupilRwd{goodSubs(iSub),1}(:,1:pupilTrialTime); stimPupilRwd{goodSubs(iSub),2}(:,1:pupilTrialTime)];
%     nullPupilConcat{iSub} = [nullPupilRwd{goodSubs(iSub),1}(:,1:pupilTrialTime); nullPupilRwd{goodSubs(iSub),2}(:,1:pupilTrialTime)];

    %subtract mean null trial
    pupilMean(iSub,:) = squeeze(nanmean(pupilConcat{iSub},1));%mean across trials
    subPupilMean = repmat(pupilMean(iSub,:),[size(pupilConcat{iSub},1),1]);%repeat X stim trials
    if subtractMean
        pupilConcat{iSub} = pupilConcat{iSub} - subPupilMean;
    end

%     pupilMeanNull(iSub,:) = squeeze(nanmean(nullPupilConcat{iSub},1));%mean across trials
%     subPupilMeanNull = repmat(pupilMeanNull(iSub,:),[size(stimPupilConcat{iSub},1),1]);%repeat X stim trials
%     if subtractMeanNull
%         stimPupilConcat{iSub} = stimPupilConcat{iSub} - subPupilMeanNull;
%     end
    
%pupil amplitude measures and stats
    subStdPupil{iSub} = nanstd(pupilConcat{iSub},0,2);%per trial
%     subBaseStimPupil{iSub} = nanmean(stimPupilConcat{iSub}(:,1:baselineT),2);%per trial
    temp = fillmissing(pupilConcat{iSub},'linear',2);
    subBasePupil{iSub} = nanmean(temp(:,1:baselineT),2);%per trial
    f=fft(fillmissing(pupilConcat{iSub},'linear',2),[],2);
    subPhPupil{iSub} = angle(f(:,2));
    subAmpPupil{iSub} = abs(f(:,2));
   
    subPupilStdMedian(iSub) = nanmedian(subStdPupil{iSub});
    subPupilBaseMedian(iSub) = nanmedian(subBasePupil{iSub});
    
    %compute mean pupil
    meanPupilTrial(iSub,:) = squeeze(nanmean(fillmissing(pupilConcat{iSub},'linear',2),1));
    %regress mean trial
    A = meanPupilTrial(iSub,:);
    B = fillmissing(pupilConcat{iSub},'linear',2);
    subBetasPupil{iSub} = A'\B';%A'x=B

    %pulse stats
    subPulse{iSub} = [rwdPulseTrials{goodSubs(iSub),1} rwdPulseTrials{goodSubs(iSub),2}];%concatenate pulse of all trials

    %subtract mean null trial
    pulseMean(iSub,:) = squeeze(nanmean(subPulse{iSub},2));%mean across trials
    subPulseMean = repmat(pulseMean(iSub,:),[size(subPulse{iSub},2),1]);%repeat X stim trials
%     subPulseMean{iSub} = mean(subPulse{iSub});
    if subtractMean
        subPulse{iSub} = subPulse{iSub} - subPulseMean';
    end
    
    subPulseStd{iSub} = std(fillmissing(subPulse{iSub},'linear',2));
%     subStimPulseStd{iSub} = std(subStimPulse{iSub});

    f=fft(fillmissing(subPulse{iSub},'linear',2),[],2);
%     f=fft(subStimPulse{iSub});
    subPulsePh{iSub} = angle(f(2,:));
    subPulseAmp{iSub} = abs(f(2,:));
        %compute mean pulse trials
    meanPulseTrial(iSub,:) = squeeze(nanmean(fillmissing(subPulse{iSub},'linear',2),2));
    %regress mean trial
    A = meanPulseTrial(iSub,:);
    B = fillmissing(subPulse{iSub},'linear',2);
%     B = stimPupilConcat{iSub};
    subBetasPulse{iSub} = A'\B;%A'x=B
    
%         subNullPulseStd{iSub} = std(subNullPulse{iSub});
%     subNullPulseMean{iSub} = mean(subNullPulse{iSub});
    %             f=fft(fillmissing(subNullPulse{iSub},'linear',2),[],2);
%     f=fft(subNullPulse{iSub});
%     subNullPulsePh{iSub} = angle(f(2,:));
%     subNullPulseAmp{iSub} = abs(f(2,:));
    
    
    subResp{iSub} = [rwdRvTrials{goodSubs(iSub),1} rwdRvTrials{goodSubs(iSub),2}];%concatenate respiration 
    subRespStd{iSub} = std(subResp{iSub});
    subRespMean{iSub} = mean(subResp{iSub});
    f=fft(subResp{iSub});
    subRespPh{iSub} = angle(f(2,:));
    subRespAmp{iSub} = abs(f(2,:));
    
    
    
    %get mean pupil for high and low contrast, and for frequencies

    contrastTrialsConcat{iSub} = [contrastTrials{goodSubs(iSub),1} contrastTrials{goodSubs(iSub),2}];
    contrastTrialsConcat{iSub} = contrastTrialsConcat{iSub}+1;%so null trials are 1, not 0. 
    for icontrast=1:max(contrastTrialsConcat{iSub})
        contrastPupil{iSub,icontrast} = pupilConcat{iSub}(contrastTrialsConcat{iSub}==icontrast,:);
        meanContrastPupil(iSub,icontrast,:) = nanmean(contrastPupil{iSub,icontrast});
    end
    
    freqTrialsConcat{iSub} = [freqTrials{goodSubs(iSub),1} freqTrials{goodSubs(iSub),2}];
    freqTrialsConcat{iSub} = freqTrialsConcat{iSub}+1;%so null trials are 1, not 0.
    for ifreq=1:max(freqTrialsConcat{iSub})
        freqPupil{iSub,ifreq} = pupilConcat{iSub}(freqTrialsConcat{iSub}==ifreq,:);
        meanFreqPupil(iSub,ifreq,:) = nanmean(freqPupil{iSub,ifreq});
    end
    
end

%%



for iSub=1:length(goodSubs)
    for iRoi=bensonROIs%1:length(roiNames) %length(bensonROIs)
        for ibin=1:nbins
            roiTseries = [roiTC{goodSubs(iSub),iRoi,1}.tSeries roiTC{goodSubs(iSub),iRoi,2}.tSeries];
            goodVoxels = ~isnan(mean(roiTseries,2));
            
            
            binVoxels = eccen{goodSubs(iSub),iRoi}>binBorders(ibin) & eccen{goodSubs(iSub),iRoi}<=binBorders(ibin+1);
            %                                 binVoxels = binVoxels & areas{goodSubs(iSub),iRoi}==1;%ONLY V1
            binVoxels = goodVoxels & binVoxels & areas{goodSubs(iSub),iRoi}>0;%ONLY V1,V2,V3
            numBinVoxels(iSub,iRoi,ibin) = sum(binVoxels);
            
            %             %remove voxels that include NaNs
            %             binTseries = [roiTC{goodSubs(iSub),iRoi,1}.tSeries(binVoxels,:) roiTC{goodSubs(iSub),iRoi,2}.tSeries(binVoxels,:)];
            %             binTseries = binTseries(~isnan(mean(binTseries,2)),:);
            
            %             %update number of voxels
            %             binVoxInd = binVoxInd(~isnan(mean(binTseries,2)));
            %             numBinVoxels(iSub,iRoi,ibin) = length(binVoxInd);
            
            for rwd=1:2
%                 binTseriesRwd = roiTC{goodSubs(iSub),iRoi,rwd}.tSeries(binVoxels,:);
%                 binNullTseriesRwd{iSub,iRoi,ibin,rwd} = binTseriesRwd(:,nullTrialsTRs{goodSubs(iSub),rwd}==1);
                binTseriesRwd{iSub,iRoi,ibin,rwd} = roiTC{goodSubs(iSub),iRoi,rwd}.tSeries(binVoxels,:);
            end
            
            %concatenate across rwd
%             binNullTseries{iSub,iRoi,ibin} = [ binNullTseriesRwd{iSub,iRoi,ibin,1} binNullTseriesRwd{iSub,iRoi,ibin,2}];
            binTseries{iSub,iRoi,ibin} = [ binTseriesRwd{iSub,iRoi,ibin,1} binTseriesRwd{iSub,iRoi,ibin,2}];
            %reshape into trials
%             binNullTrialTseries{iSub,iRoi,ibin} = reshape(binNullTseries{iSub,iRoi,ibin},numBinVoxels(iSub,iRoi,ibin),trialLength,[]);%vox, T, trials
            binTrialTseries{iSub,iRoi,ibin} = reshape(binTseries{iSub,iRoi,ibin},numBinVoxels(iSub,iRoi,ibin),trialLength,[]);%vox, T, trials
            
            %subtract mean null trial
            voxMeanTrial{iSub,iRoi,ibin} = squeeze(mean(binTrialTseries{iSub,iRoi,ibin},3));%mean across trials
            subBinMeanTrial = repmat(voxMeanTrial{iSub,iRoi,ibin},[1,1,size(binTrialTseries{iSub,iRoi,ibin},3)]);%repeat X stim trials
            if subtractMean
                binTrialTseries{iSub,iRoi,ibin} = binTrialTseries{iSub,iRoi,ibin} - subBinMeanTrial;
            end
            
            %get fMRI response amplitude
            %STD
%             subBinStdNull{iSub,iRoi,ibin} = squeeze(std(binNullTrialTseries{iSub,iRoi,ibin},0,2));%vox, trial
            subBinStd{iSub,iRoi,ibin} = squeeze(std(binTrialTseries{iSub,iRoi,ibin},0,2));%vox, trial
            %FFT amplitude
            f = fft(binTrialTseries{iSub,iRoi,ibin},[],2);
            subBinAmp{iSub,iRoi,ibin} = squeeze(abs(f(:,2,:)));
            
            %compute mean stim trial per voxel
            binMeanTrial{iSub,iRoi,ibin} = squeeze(nanmean(binTrialTseries{iSub,iRoi,ibin},3));%vox, T
            %regress mean trial
            ntrials = size(binTrialTseries{iSub,iRoi,ibin},3);
            for ivox=1:numBinVoxels(iSub,iRoi,ibin)
                
                A = binMeanTrial{iSub,iRoi,ibin}(ivox,:);
                A = zscore(A,0,2);
                B = squeeze(binTrialTseries{iSub,iRoi,ibin}(ivox,:,:));
                subBinBetas{iSub,iRoi,ibin}(ivox,:) = A'\B;%A'x=B
                %                 Y = squeeze(binStimTrialTseries{iSub,iRoi,ibin}(ivox,:,:));
                %                 b = regress(Y,X);%Y=X*b
            end
            
            %             subBinBetasStim{iSub,iRoi,ibin} = x;
            
            %             stimFreqTrials{iSub} = [freqTrials{goodSubs(iSub),1}(stimTrials{goodSubs(iSub),1}==1) freqTrials{goodSubs(iSub),2}(stimTrials{goodSubs(iSub),2}==1)];
            for ifreq=1:max(freqTrials{goodSubs(iSub),1})
                subBinFreqStd{iSub,iRoi,ibin,ifreq} =  subBinStd{iSub,iRoi,ibin}(:,freqTrialsConcat{iSub}==ifreq,:);%vox,trial
            end
            
            %             stimContrastTrials{iSub} =  [contrastTrials{goodSubs(iSub),1}(stimTrials{goodSubs(iSub),1}==1) contrastTrials{goodSubs(iSub),2}(stimTrials{goodSubs(iSub),2}==1)];
            for icontrast=1:max(contrastTrials{goodSubs(iSub),1})
                subBinContrastStd{iSub,iRoi,ibin,icontrast} =  subBinStd{iSub,iRoi,ibin}(:,contrastTrialsConcat{iSub}==icontrast,:);%vox,trial
            end

            switch arousalType
                case 1
                    trialArousal = subBasePupil{iSub}';
                    medianArousal = subPupilBaseMedian(iSub);
                case 2
                    trialArousal = subStdPupil{iSub}';
                    medianArousal = subPupilStdMedian(iSub);
                case 3
                    trialArousal = subStdPulse{iSub}';
                    medianArousal = subPulseStdMedian(iSub);   
            end
                
            for arousal=1:2
                if arousal==1
                    arousalTrials{iSub,arousal} = trialArousal>medianArousal;
                else
                    arousalTrials{iSub,arousal} = trialArousal<medianArousal;
                end
            end
            
        end
    end
end

ifig=0;
%% DECODING
for featureType=featureTypes
    for iSub=1:length(goodSubs)
        for iRoi=bensonROIs%1:length(roiNames) %length(bensonROIs)
            for ibin=1:nbins
                
                %decode frequency and contrast
                switch featureType
                    case 1
                        trialDataFmri = subBinStd{iSub,iRoi,ibin};
                        trialDataPupil = subStdPupil{iSub}';
                        trialDataPulse = subPulseStd{iSub};
                    case 2
                        trialDataFmri = subBinAmp{iSub,iRoi,ibin};
                        trialDataPupil = subAmpPupil{iSub}';
                        trialDataPulse = subPulseAmp{iSub};
                    case 3
                        trialDataFmri = subBinBetas{iSub,iRoi,ibin};
                        trialDataPupil = subBetasPupil{iSub};
                        trialDataPulse = subBetasPulse{iSub};
                end
                
                
                trialData =  trialDataFmri;
                if includePupil
                    trialData = [trialData; trialDataPupil];
                end
                if includePupilBase
                    trialData = [trialDataFmri; subBasePupil{iSub}'];
                end
                if includePulse
                    trialData = [trialData; trialDataPulse];
                end
                
                
                
                
                
                
                
                
                
                
                
                for arousal=1:2
                    arousalData = trialData(:,arousalTrials{iSub,arousal});
                    arousalFreqLabels = stimFreqTrials{iSub}(arousalTrials{iSub,arousal});
                    arousalContrastLabels = stimContrastTrials{iSub}(arousalTrials{iSub,arousal});
                    ntrials = size(arousalData,2);
                    for itrial=1:ntrials
                        trainData =  arousalData(:,[1:itrial-1 itrial+1:ntrials]);
                        trainLabelsFreq = arousalFreqLabels([1:itrial-1 itrial+1:ntrials]);
                        testData = arousalData(:,itrial);
                        classGuessFreq{iSub,iRoi,ibin,arousal}(itrial) = classify(testData', trainData', trainLabelsFreq, 'diagLinear');
                        trainLabelsContrast = arousalContrastLabels([1:itrial-1 itrial+1:ntrials]);
                        classGuessContrast{iSub,iRoi,ibin,arousal}(itrial) = classify(testData', trainData', trainLabelsContrast, 'diagLinear');
                    end
                    classAccFreq{featureType,arousal}(iSub,iRoi,ibin) = sum(arousalFreqLabels==classGuessFreq{iSub,iRoi,ibin,arousal})/ntrials;%make sure not to include NaN trials?
                    classAccContrast{featureType,arousal}(iSub,iRoi,ibin) = sum(arousalContrastLabels==classGuessContrast{iSub,iRoi,ibin,arousal})/ntrials;%make sure not to include NaN trials?
                end
                
                
                
                
                
                
                
                
                ntrials = size(trialData,2);
                for itrial=1:ntrials
                    trainData =  trialData(:,[1:itrial-1 itrial+1:ntrials]);
                    trainLabelsFreq = freqTrialsConcat{iSub}([1:itrial-1 itrial+1:ntrials]);
                    testData = trialData(:,itrial);
                    %                testLabelFreq = stimFreqTrials{iSub}(itrial);
                    classGuessFreq{iSub,iRoi,ibin}(itrial) = classify(testData', trainData', trainLabelsFreq, 'diagLinear');
%                     classGuessFreq{iSub,iRoi,ibin}(itrial) = classify(testDataFreq', trainDataFreq', trainLabelsFreq, 'diagquadratic');%
                    
%                     trainData =  trialData(:,[1:itrial-1 itrial+1:ntrials]);
                    trainLabelsContrast = contrastTrialsConcat{iSub}([1:itrial-1 itrial+1:ntrials]);
%                     testData = trialData(:,itrial);
                    classGuessContrast{iSub,iRoi,ibin}(itrial) = classify(testData', trainData', trainLabelsContrast, 'diagLinear');
%                     classGuessContrast{iSub,iRoi,ibin}(itrial) = classify(testDataContrast', trainDataContrast', trainLabelsContrast, 'diagquadratic');
                end
%                 classConfMatFreq{iSub,iRoi,ibin} = confusionmat(stimFreqTrials{iSub}, classGuessFreq{iSub,iRoi,ibin});
                stimAccFreq{featureType}(iSub,iRoi,ibin) = sum(freqTrialsConcat{iSub}(freqTrialsConcat{iSub}>1)==classGuessFreq{iSub,iRoi,ibin}(freqTrialsConcat{iSub}>1))/sum(freqTrialsConcat{iSub}>1);%make sure not to include NaN trials?
                accFreq{featureType}(iSub,iRoi,ibin) = sum(freqTrialsConcat{iSub}==classGuessFreq{iSub,iRoi,ibin})/ntrials;%make sure not to include NaN trials?

%                 classConfMatContrast{iSub,iRoi,ibin} = confusionmat(stimContrastTrials{iSub}, classGuessContrast{iSub,iRoi,ibin});
               stimAccContrast{featureType}(iSub,iRoi,ibin) = sum(contrastTrialsConcat{iSub}(contrastTrialsConcat{iSub}>1)==classGuessContrast{iSub,iRoi,ibin}(contrastTrialsConcat{iSub}>1))/sum(contrastTrialsConcat{iSub}>1);%make sure not to include NaN trials?                
                accContrast{featureType}(iSub,iRoi,ibin) = sum(contrastTrialsConcat{iSub}==classGuessContrast{iSub,iRoi,ibin})/ntrials;%make sure not to include NaN trials?                


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
    plot(squeeze(mean(accFreq{featureType})),'color',plotColors{featureType})
    hold on
%     mean(squeeze(mean(classAccFreq{featureType})))
end
legend('std','amp','regress');
for featureType=featureTypes
    plot(squeeze(mean(accContrast{featureType})),'color',plotColors{featureType})
%     mean(squeeze(mean(classAccContrast{featureType})))
end
xlabel('eccentricity bin');
ylabel('proportion correct');
% legend('frequency','contrast');
% hline(1/nfreqs);
hline(0.25);
% hline(1/ncontrasts);
hline((1-0.25)/2);
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

subplot(rows,cols,2)
for featureType=featureTypes
    plot(squeeze(mean(accFreq{featureType},3)),'color',plotColors{featureType})
    hold on
    plot(squeeze(mean(accContrast{featureType},3)),'color',plotColors{featureType})
end
xlabel('subject #');
% legend('frequency','contrast');
% hline(1/nfreqs);
hline(0.25);
% hline(1/ncontrasts);
hline((1-0.25)/2);

old = {' ', ','};
new = '_';
savefig([dataFolder replace(titleStr,old, new)]);

[max(squeeze(mean(accContrast{featureType},1))), max(squeeze(mean(accFreq{featureType},1)))]

% savefig(replace(titleStr,old, new));



% ifig=ifig+1; figure(ifig); clf
% rows=1;
% cols=2;
% subplot(rows,cols,1)
% for featureType=featureTypes
%     plot(squeeze(mean(stimAccFreq{featureType})),'color',plotColors{featureType})
%     hold on
% %     mean(squeeze(mean(classAccFreq{featureType})))
% end
% legend('std','amp','regress');
% for featureType=featureTypes
%     plot(squeeze(mean(stimAccContrast{featureType})),'color',plotColors{featureType})
% %     mean(squeeze(mean(classAccContrast{featureType})))
% end
% xlabel('eccentricity bin');
% ylabel('proportion correct');
% % legend('frequency','contrast');
% hline(1/nfreqs);
% hline(1/ncontrasts);
% titleStr = 'decode ';
% if includePupil
%     titleStr = [titleStr 'pupil, '];
% end
% if includePulse
%     titleStr = [titleStr 'pulse, '];
% end
% if includePupilBase
%     titleStr = [titleStr 'pupil baseline, '];
% end
% if subtractMean
%     titleStr = [titleStr 'mean subtracted '];
% end
% 
% 
% title(titleStr);
% 
% subplot(rows,cols,2)
% for featureType=featureTypes
%     plot(squeeze(mean(stimAccFreq{featureType},3)),'color',plotColors{featureType})
%     hold on
%     plot(squeeze(mean(stimAccContrast{featureType},3)),'color',plotColors{featureType})
% end
% xlabel('subject #');
% % legend('frequency','contrast');
% hline(1/nfreqs);
% hline(1/ncontrasts);
% 
% old = {' ', ','};
% new = '_';
% savefig([dataFolder replace(titleStr,old, new)]);
% 
% [max(squeeze(mean(stimAccContrast{featureType},1))), max(squeeze(mean(stimAccFreq{featureType},1)))]




toc

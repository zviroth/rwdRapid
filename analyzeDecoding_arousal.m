close all
mrQuit
clear all
ifig=0;
arousalTypes = 1:5;%1-pupil baseline, 2-pupil std, 3-pulse std, 4-rwd, 5-DMN betas
includeControl=0;%left and right DMN trial amplitude
includePulse = 0;
subtractMean = 1;
includePupil = 0;
includePupilBase = 0;
featureTypes = 1:3;
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
    'trialsPerRun', 'trialLength', ...
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
    'frames', 'junkedFrames', 'TR', 'trialsPerRun', 'trialLength', 'nVolumes',...
    'eccen','ang','areas',...
    'numVox','roiTC','nullTrialsTRs',...
    'nullTrialsRun','nullTrials','contrastTrialsRun','freqTrialsRun','contrastTrials','freqTrials',...
    'nullTseries','nullTrialTseries','stimTseries','stimTrialTseries');
for iSub=1:size(nullTrials,1)
    for rwd=1:2
        stimTrials{iSub,rwd} = 1-nullTrials{iSub,rwd};
        freqTrials{iSub,rwd} = freqTrials{iSub,rwd}.*stimTrials{iSub,rwd};%set null trials to 0
        contrastTrials{iSub,rwd} = contrastTrials{iSub,rwd}.*stimTrials{iSub,rwd};%set null trials to 0
    end
end
goodSubs = [1:12 14];
% goodSubs = 1:2;
%%
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
    
    stimPupilConcat{iSub} = [stimPupilRwd{goodSubs(iSub),1}(:,1:pupilTrialTime); stimPupilRwd{goodSubs(iSub),2}(:,1:pupilTrialTime)];
    nullPupilConcat{iSub} = [nullPupilRwd{goodSubs(iSub),1}(:,1:pupilTrialTime); nullPupilRwd{goodSubs(iSub),2}(:,1:pupilTrialTime)];
    
    %subtract mean null trial
    pupilMeanNull(iSub,:) = squeeze(nanmean(nullPupilConcat{iSub},1));%mean across trials
    subPupilMeanNull = repmat(pupilMeanNull(iSub,:),[size(stimPupilConcat{iSub},1),1]);%repeat X stim trials
    if subtractMean
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
    
    
    
    
    %pulse stats
    subStimPulse{iSub} = [rwdPulseTrials{goodSubs(iSub),1}(:,stimTrials{goodSubs(iSub),1}==1) rwdPulseTrials{goodSubs(iSub),2}(:,stimTrials{goodSubs(iSub),2}==1)];%concatenate pulse of null trials
    subNullPulse{iSub} = [rwdPulseTrials{goodSubs(iSub),1}(:,nullTrials{goodSubs(iSub),1}==1) rwdPulseTrials{goodSubs(iSub),2}(:,nullTrials{goodSubs(iSub),2}==1)];%concatenate pulse of null trials
    
    %subtract mean null trial
    pulseMeanNull(iSub,:) = squeeze(nanmean(subNullPulse{iSub},2));%mean across trials
    subPulseMeanNull = repmat(pulseMeanNull(iSub,:),[size(subStimPulse{iSub},2),1]);%repeat X stim trials
    if subtractMean
        subStimPulse{iSub} = subStimPulse{iSub} - subPulseMeanNull';
    end
    
    subStdStimPulse{iSub} = std(fillmissing(subStimPulse{iSub},'linear',2));
    subStimPulseStdMedian(iSub) = nanmedian(subStdStimPulse{iSub});
    
    subStimRwd{iSub} = 1:length(subStdStimPulse{iSub});
%     subStimRwdMedian(iSub) = length(freqTrials{goodSubs(iSub),1})+0.5;%median(subRwd{iSub});
    subStimRwdMedian(iSub) = sum(nullTrialsTRs{goodSubs(iSub),1}==0)/trialLength+0.5;%after the last trial of high reward
    
    
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
    
    stimContrastTrials{iSub} = [contrastTrials{goodSubs(iSub),1}(stimTrials{goodSubs(iSub),1}==1) contrastTrials{goodSubs(iSub),2}(stimTrials{goodSubs(iSub),2}==1)];
    for icontrast=1:ncontrasts
        contrastPupil{iSub,icontrast} = stimPupilConcat{iSub}(stimContrastTrials{iSub}==icontrast,:);
        meanContrastPupil(iSub,icontrast,:) = nanmean(contrastPupil{iSub,icontrast});
        contrastPulse{iSub,icontrast} = subStimPulse{iSub}(:,stimContrastTrials{iSub}==icontrast);
        meanContrastPulse(iSub,icontrast,:) = nanmean(contrastPulse{iSub,icontrast},2);
    end
    
    stimFreqTrials{iSub} = [freqTrials{goodSubs(iSub),1}(stimTrials{goodSubs(iSub),1}==1) freqTrials{goodSubs(iSub),2}(stimTrials{goodSubs(iSub),2}==1)];
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
    for iRoi=1:length(controlROIs)
        roiTseries = [roiTC{goodSubs(iSub),controlROIs(iRoi),1}.tSeries roiTC{goodSubs(iSub),controlROIs(iRoi),2}.tSeries];
        %         goodVoxels = ~isnan(mean(roiTseries,2));
        
        %         roiTseries = roiTseries(goodVoxels,:);
        for rwd=1:2
            controlTseriesRwd = roiTC{goodSubs(iSub),controlROIs(iRoi),rwd}.tSeries;
            controlNullTseriesRwd{iSub,iRoi,rwd} = controlTseriesRwd(:,nullTrialsTRs{goodSubs(iSub),rwd}==1);
            controlStimTseriesRwd{iSub,iRoi,rwd} = controlTseriesRwd(:,nullTrialsTRs{goodSubs(iSub),rwd}==0);
        end
        %concatenate across rwd
        controlNullTseries{iSub,iRoi} = [ nanmean(controlNullTseriesRwd{iSub,iRoi,1},1) nanmean(controlNullTseriesRwd{iSub,iRoi,2},1)];
        controlStimTseries{iSub,iRoi} = [ nanmean(controlStimTseriesRwd{iSub,iRoi,1},1) nanmean(controlStimTseriesRwd{iSub,iRoi,2},1)];
        %reshape into trials
        controlNullTrialTseries{iSub,iRoi} = reshape(controlNullTseries{iSub,iRoi},trialLength,[]);%T, trials
        controlStimTrialTseries{iSub,iRoi} = reshape(controlStimTseries{iSub,iRoi},trialLength,[]);%T, trials
        
        %subtract mean null trial
        controlMeanNull{iSub,iRoi} = squeeze(mean(controlNullTrialTseries{iSub,iRoi},2));%mean across trials
        subControlMeanNull = repmat(controlMeanNull{iSub,iRoi},[1,size(controlStimTrialTseries{iSub,iRoi},2)]);%repeat X stim trials
        if subtractMean
            controlStimTrialTseries{iSub,iRoi} = controlStimTrialTseries{iSub,iRoi} - subControlMeanNull;
        end
        
        controlStimStd{iSub,iRoi} = squeeze(std(controlStimTrialTseries{iSub,iRoi},0,1));
        f = fft(controlStimTrialTseries{iSub,iRoi});
        controlStimAmp{iSub,iRoi} = squeeze(abs(f(2,:)));
        controlMeanStimTrial{iSub,iRoi} = squeeze(nanmean(controlStimTrialTseries{iSub,iRoi},2));%T
        ntrials = size(controlStimTrialTseries{iSub,iRoi},2);
        %         for ivox=1:sum(goodVoxels)
        A = controlMeanStimTrial{iSub,iRoi};
        A = zscore(A);
        B = controlStimTrialTseries{iSub,iRoi};
        controlStimBetas{iSub,iRoi} = A\B;%A'x=B
        %         end
        %         subStimControlStdMedian{iSub,iRoi} = nanmedian(controlStimStd{iSub,iRoi});
        
        controlStimBetasMedian(iSub) = median(controlStimBetas{iSub,1});
    end
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
            if subtractMean
                binStimTrialTseries{iSub,iRoi,ibin} = binStimTrialTseries{iSub,iRoi,ibin} - subBinMeanNull;
            end
            
            
            
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
            for ifreq=1:max(freqTrials{goodSubs(iSub),1})
                subBinFreqStd{iSub,iRoi,ibin,ifreq} =  subBinStdStim{iSub,iRoi,ibin}(:,stimFreqTrials{iSub}==ifreq,:);%vox,trial
            end
            
            %             stimContrastTrials{iSub} =  [contrastTrials{goodSubs(iSub),1}(stimTrials{goodSubs(iSub),1}==1) contrastTrials{goodSubs(iSub),2}(stimTrials{goodSubs(iSub),2}==1)];
            for icontrast=1:max(contrastTrials{goodSubs(iSub),1})
                subBinContrastStd{iSub,iRoi,ibin,icontrast} =  subBinStdStim{iSub,iRoi,ibin}(:,stimContrastTrials{iSub}==icontrast,:);%vox,trial
            end
        end
    end
end

histEdgesFreq = [1:nfreqs+1]-0.5;
histEdgesContrast = [1:ncontrasts+1]-0.5;



%arousal
for arousalType=arousalTypes
    for iSub=1:length(goodSubs)
        
        switch arousalType
            case 1
                trialArousal = subBaseStimPupil{iSub}';
                medianArousal = subStimPupilBaseMedian(iSub);
            case 2
                trialArousal = subStdStimPupil{iSub}';
                medianArousal = subStimPupilStdMedian(iSub);
            case 3
                trialArousal = subStdStimPulse{iSub}';
                medianArousal = subStimPulseStdMedian(iSub);
            case 4
                trialArousal = subStimRwd{iSub}';
                medianArousal = subStimRwdMedian(iSub);
            case 5
                trialArousal = controlStimBetas{iSub,1}';
                medianArousal = controlStimBetasMedian(iSub);
        end
        
        for arousal=1:2
            if arousal==1
                arousalTrials{iSub,arousal} = trialArousal>medianArousal;
            else
                arousalTrials{iSub,arousal} = trialArousal<medianArousal;
            end
            
            
        end
        
    end
    
    
    % ifig=0;
    %% DECODING
    for featureType=featureTypes
        for iSub=1:length(goodSubs)
            for iRoi=bensonROIs%1:length(roiNames) %length(bensonROIs)
                for ibin=1:nbins
                    
                    %decode frequency and contrast
                    switch featureType
                        case 1
                            trialDataFmri = subBinStdStim{iSub,iRoi,ibin};
                            trialDataPupil = subStdStimPupil{iSub}';
                            trialDataPulse = subStdStimPulse{iSub};
                            for controlRoi=1:length(controlROIs)
                                trialDataControl{controlRoi} = controlStimStd{iSub,controlRoi};
                            end
                        case 2
                            trialDataFmri = subBinAmpStim{iSub,iRoi,ibin};
                            trialDataPupil = subAmpStimPupil{iSub}';
                            trialDataPulse = subStimPulseAmp{iSub};
                            for controlRoi=1:length(controlROIs)
                                trialDataControl{controlRoi} = controlStimAmp{iSub,controlRoi};
                            end
                        case 3
                            trialDataFmri = subBinBetasStim{iSub,iRoi,ibin};
                            trialDataPupil = subBetasStimPupil{iSub};
                            trialDataPulse = subBetasStimPulse{iSub};
                            for controlRoi=1:length(controlROIs)
                                trialDataControl{controlRoi} = controlStimBetas{iSub,controlRoi};
                            end
                    end
                    
                    
                    trialData =  trialDataFmri;
                    if includePupil
                        trialData = [trialData; trialDataPupil];
                    end
                    if includePupilBase
                        trialData = [trialDataFmri; subBaseStimPupil{iSub}'];
                    end
                    if includePulse
                        trialData = [trialData; trialDataPulse];
                    end
                    if includeControl
                        for controlRoi=1:length(controlROIs)
                            trialData = [trialData; trialDataControl{controlRoi}];
                        end
                    end
                    for arousal=1:2
                        arousalData = trialData(:,arousalTrials{iSub,arousal});
                        arousalFreqLabels = stimFreqTrials{iSub}(arousalTrials{iSub,arousal});
                        arousalContrastLabels = stimContrastTrials{iSub}(arousalTrials{iSub,arousal});
                        ntrials = size(arousalData,2);
                        classGuessFreq{iSub,iRoi,ibin,arousal}=[];
                        classGuessContrast{iSub,iRoi,ibin,arousal}=[];
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
                        
                        
                        freqCondsTrue(iSub,arousal,:) = histcounts(stimFreqTrials{iSub}(arousalTrials{iSub,arousal}),histEdgesFreq);
                        contrastCondsTrue(iSub,arousal,:) = histcounts(stimContrastTrials{iSub}(arousalTrials{iSub,arousal}),histEdgesContrast);
                        
                        freqCondsGuess(featureType,iSub,iRoi,ibin,arousal,:) = histcounts(classGuessFreq{iSub,iRoi,ibin,arousal},histEdgesFreq);
                        contrastCondsGuess(featureType,iSub,iRoi,ibin,arousal,:) = histcounts(classGuessContrast{iSub,iRoi,ibin,arousal},histEdgesContrast);
                    end
                    
                end
            end
        end
    end
    

    
    %%
    plotColors = { [1 0 0], [0 0 1], [0 1 0], [0.5 1 0.2], [1 0.2 0.5] };
    lineStyles = {'-','--'};
    ifig=ifig+1; figure(ifig); clf
    rows=1;
    cols=4;
    subplot(rows,cols,1)
    pvalThresh=0.05;
    markerSize=50;
    for featureType=featureTypes
        for arousal=1:2
            plot(squeeze(mean(classAccFreq{featureType,arousal})),lineStyles{arousal},'color',plotColors{featureType})
            hold on
        end
        
        %     mean(squeeze(mean(classAccFreq{featureType})))
    end
    % legend('std','amp','regress');
    legend('high','low');
    for featureType=featureTypes
        for arousal=1:2
            plot(squeeze(mean(classAccContrast{featureType,arousal})),lineStyles{arousal},'color',plotColors{featureType})
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
    
    subplot(rows,cols,2)
    for featureType=featureTypes
        for arousal=1:2
            plot(squeeze(mean(classAccFreq{featureType,arousal},3)),lineStyles{arousal},'color',plotColors{featureType})
            hold on
            plot(squeeze(mean(classAccContrast{featureType,arousal},3)),lineStyles{arousal},'color',plotColors{featureType})
        end
    end
    xlabel('subject #');
    % legend('frequency','contrast');
    hline(1/nfreqs);
    hline(1/ncontrasts);
    
    old = {' ', ','};
    new = '_';
    % savefig([dataFolder replace(titleStr,old, new)]);
    for arousal=1:2
        [max(squeeze(mean(classAccContrast{featureType,arousal},1))), max(squeeze(mean(classAccFreq{featureType,arousal},1)))]
    end
    
    
    
    subplot(rows,cols,3)
    for arousal=1:2
        bar((arousal-1)*nfreqs+[1:5],squeeze(mean(freqCondsTrue(:,arousal,:))),'facecolor',plotColors{arousal},'linestyle',lineStyles{arousal});
        hold on
        er=errorbar((arousal-1)*nfreqs+[1:nfreqs],squeeze(mean(freqCondsTrue(:,arousal,:))),squeeze(std(freqCondsTrue(:,arousal,:))));
                er.Color = [0 0 0];
        er.LineStyle = 'none';
%         histogram(squeeze(mean(freqCondsTrue(:,arousal,:))),nfreqs+1,'facecolor',plotColors{arousal},'linestyle',lineStyles{arousal});
        
    end
    subplot(rows,cols,4)
    for arousal=1:2
        bar((arousal-1)*ncontrasts+[1:ncontrasts],squeeze(mean(contrastCondsTrue(:,arousal,:))),'facecolor',plotColors{arousal},'linestyle',lineStyles{arousal});
        hold on
        er=errorbar((arousal-1)*ncontrasts+[1:ncontrasts],squeeze(mean(contrastCondsTrue(:,arousal,:))),squeeze(std(contrastCondsTrue(:,arousal,:))));
        er.Color = [0 0 0];
        er.LineStyle = 'none';
%         histogram(squeeze(mean(contrastCondsTrue(:,arousal,:))),ncontrasts+1,'facecolor',plotColors{arousal},'linestyle',lineStyles{arousal});
        
    end
    set(gcf,'position',[100 100 1100 400])
end
% savefig(replace(titleStr,old, new));
toc

figure
rows=1;
cols=4;
subplot(rows,cols,1)
plot(squeeze(mean(meanContrastPupil))')
subplot(rows,cols,2)
plot(squeeze(mean(meanContrastPulse))')
subplot(rows,cols,3)
plot(squeeze(mean(meanFreqPupil))')
subplot(rows,cols,4)
plot(squeeze(mean(meanFreqPulse))')
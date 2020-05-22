close all
mrQuit
clear all

onlyCorrect=0;%0=all trials, 1=only correct, -1=only incorrect, 2=valid response
% for onlyCorrect=[0 1 2]
ifig=0;
arousalTypes = [1:2 4:5];%1-pupil baseline, 2-pupil std, 3-pulse std, 4-rwd, 5-DMN betas
% arousalTypes = 4;
subtractMean = 0;

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
plotColors = { [1 0 0], [0 0 1], [0 1 0], [0.5 1 0.2], [1 0.2 0.5] };
lineStyles = {'-','--'};

% dataFolder = '/Volumes/MH02086153MACDT-Drobo/rwdFmri/';
dataFolder = '/Users/rothzn/rwdFmri/';
% dataFolder = 'c:\rwdFmri\';
freqs = logspace(-0.3,0.5,5);
contrasts = logspace(-0.7,0,2);
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
% for iSub=1:size(nullTrials,1)
%     for rwd=1:2
%         stimTrials{iSub,rwd} = 1-nullTrials{iSub,rwd};
%         freqTrials{iSub,rwd} = freqTrials{iSub,rwd}.*stimTrials{iSub,rwd};%set null trials to 0
%         contrastTrials{iSub,rwd} = contrastTrials{iSub,rwd}.*stimTrials{iSub,rwd};%set null trials to 0
%     end
% end

% goodSubs = [1:12 14];
goodSubs = [1:14];
% goodSubs = 1:2;
%%
nullPupilRwd = origNullPupil;
stimPupilRwd = origStimPupil;
for iSub=1:length(goodSubs)
    
    for rwd=1:2
        %find trials with correct type of response
        temp = trialCorrectness{goodSubs(iSub),rwd}(:,2:end);%first trial is junked
        rwdTrialCorrectness{iSub,rwd} = temp(:);
        allTrialResponse = [];
        for r=1:size(trialResponse,3)
            allTrialResponse = [allTrialResponse; trialResponse{iSub,rwd,r}(2:end)];
        end
        rwdTrialResponse{iSub,rwd} = allTrialResponse;
        %         temp = trialResponse{goodSubs(iSub),rwd}(:,2:end);%first trial is junked
        %         rwdTrialResponse{iSub,rwd} = temp(:);
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
        correctTrials{iSub,rwd} = correctTrials{iSub,rwd}';
        temp = repmat(correctTrials{iSub,rwd},trialLength,1);
        correctTRs{iSub,rwd} = temp(:);
        
        %this is for all data
        stimTrials{iSub,rwd} = 1-nullTrials{goodSubs(iSub),rwd};
        freqTrials{goodSubs(iSub),rwd} = freqTrials{goodSubs(iSub),rwd}.*stimTrials{iSub,rwd}.*correctTrials{iSub,rwd};%set null trials to 0
        contrastTrials{goodSubs(iSub),rwd} = contrastTrials{goodSubs(iSub),rwd}.*stimTrials{iSub,rwd}.*correctTrials{iSub,rwd};%set null trials to 0
        
        %new
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
    
    %     stimPupilConcat{iSub} = [stimPupilRwd{goodSubs(iSub),1}(:,1:pupilTrialTime); stimPupilRwd{goodSubs(iSub),2}(:,1:pupilTrialTime)];
    %     nullPupilConcat{iSub} = [nullPupilRwd{goodSubs(iSub),1}(:,1:pupilTrialTime); nullPupilRwd{goodSubs(iSub),2}(:,1:pupilTrialTime)];
    
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
    subNullPupilStdMedian(iSub) = nanmedian(subStdNullPupil{iSub});
    subStimPupilBaseMedian(iSub) = nanmedian(subBaseStimPupil{iSub});
    subNullPupilBaseMedian(iSub) = nanmedian(subBaseNullPupil{iSub});
    
    %compute mean pupil for stim trials
    meanStimPupilTrial(iSub,:) = squeeze(nanmean(fillmissing(stimPupilConcat{iSub},'linear',2),1));
    
    %         %regress mean trial
    %         A = zscore(meanStimPupilTrial(iSub,:));
    %         B = fillmissing(stimPupilConcat{iSub},'linear',2);
    %         %     B = stimPupilConcat{iSub};
    %         subBetasStimPupil{iSub} = A'\B';%A'x=B
    
    
    
    
    %pulse stats
    for rwd=1:2
        stimTrials{iSub,rwd} = stimTrials{iSub,rwd}.*correctTrials{iSub,rwd};
        nullTrials{goodSubs(iSub),rwd} = nullTrials{goodSubs(iSub),rwd}.*correctTrials{iSub,rwd};
    end
    subStimPulse{iSub} = [rwdPulseTrials{goodSubs(iSub),1}(:,stimTrials{iSub,1}==1) rwdPulseTrials{goodSubs(iSub),2}(:,stimTrials{iSub,2}==1)];%concatenate pulse of null trials
    subNullPulse{iSub} = [rwdPulseTrials{goodSubs(iSub),1}(:,nullTrials{goodSubs(iSub),1}==1) rwdPulseTrials{goodSubs(iSub),2}(:,nullTrials{goodSubs(iSub),2}==1)];%concatenate pulse of null trials
    
    subStimPpg{iSub} = [rwdPpgTrials{goodSubs(iSub),1}(:,stimTrials{iSub,1}==1) rwdPpgTrials{goodSubs(iSub),2}(:,stimTrials{iSub,2}==1)];%concatenate pulse of null trials
    subNullPpg{iSub} = [rwdPpgTrials{goodSubs(iSub),1}(:,nullTrials{goodSubs(iSub),1}==1) rwdPpgTrials{goodSubs(iSub),2}(:,nullTrials{goodSubs(iSub),2}==1)];%concatenate pulse of null trials
    
    subStimResp{iSub} = [rwdRvTrials{goodSubs(iSub),1}(:,stimTrials{iSub,1}==1) rwdRvTrials{goodSubs(iSub),2}(:,stimTrials{iSub,2}==1)];%concatenate pulse of null trials
    subNullResp{iSub} = [rwdRvTrials{goodSubs(iSub),1}(:,nullTrials{goodSubs(iSub),1}==1) rwdRvTrials{goodSubs(iSub),2}(:,nullTrials{goodSubs(iSub),2}==1)];%concatenate pulse of null trials
    
    
    %subtract mean null trial
    pulseMeanNull(iSub,:) = squeeze(nanmean(subNullPulse{iSub},2));%mean across trials
    subPulseMeanNull = repmat(pulseMeanNull(iSub,:),[size(subStimPulse{iSub},2),1]);%repeat X stim trials
    if subtractMean
        subStimPulse{iSub} = subStimPulse{iSub} - subPulseMeanNull';
    end
    
    subStdStimPulse{iSub} = std(fillmissing(subStimPulse{iSub},'linear',2));
    subStimPulseStdMedian(iSub) = nanmedian(subStdStimPulse{iSub});
    subStdNullPulse{iSub} = std(fillmissing(subNullPulse{iSub},'linear',2));
    subNullPulseStdMedian(iSub) = nanmedian(subStdNullPulse{iSub});
    
    subStimRwd{iSub} = 1:length(subStdStimPulse{iSub});
    %     subStimRwdMedian(iSub) = length(freqTrials{goodSubs(iSub),1})+0.5;%median(subRwd{iSub});
    subStimRwdMedian(iSub) = sum(stimTrials{iSub,1})+0.5;
    %     subStimRwdMedian(iSub) = sum(nullTrialsTRs{goodSubs(iSub),1}==0)/trialLength+0.5;%after the last trial of high reward
    
    subNullRwd{iSub} = 1:length(subStdNullPulse{iSub});
    subNullRwdMedian(iSub) = sum(nullTrials{iSub,1})+0.5;
    %     subStimPulseStd{iSub} = std(subStimPulse{iSub});
    subMeanStimPulse{iSub} = nanmean(subStimPulse{iSub});
    
    %         f=fft(fillmissing(subStimPulse{iSub},'linear',2),[],2);
    %         %     f=fft(subStimPulse{iSub});
    %         subStimPulsePh{iSub} = angle(f(2,:));
    %         subStimPulseAmp{iSub} = abs(f(2,:));
    %compute mean pulse for stim trials
    %         meanStimPulseTrial(iSub,:) = squeeze(nanmean(fillmissing(subStimPulse{iSub},'linear',2),2));
    %regress mean trial
    %         A = zscore(meanStimPulseTrial(iSub,:));
    %         B = fillmissing(subStimPulse{iSub},'linear',2);
    %         %     B = stimPupilConcat{iSub};
    %         subBetasStimPulse{iSub} = A'\B;%A'x=B
    
    subNullPulseStd{iSub} = std(subNullPulse{iSub});
    subNullPulseMean{iSub} = mean(subNullPulse{iSub});
    %             f=fft(fillmissing(subNullPulse{iSub},'linear',2),[],2);
    %         f=fft(subNullPulse{iSub});
    %         subNullPulsePh{iSub} = angle(f(2,:));
    %         subNullPulseAmp{iSub} = abs(f(2,:));
    
    
%     subNullResp{iSub} = [rwdRvTrials{goodSubs(iSub),1}(:,nullTrials{goodSubs(iSub),1}==1) rwdRvTrials{goodSubs(iSub),2}(:,nullTrials{goodSubs(iSub),2}==1)];%concatenate pulse of null trials
    subNullRespStd{iSub} = std(subNullResp{iSub});
    subNullRespMean{iSub} = mean(subNullResp{iSub});
    %         f=fft(subNullResp{iSub});
    %         subNullRespPh{iSub} = angle(f(2,:));
    %         subNullRespAmp{iSub} = abs(f(2,:));
    
    %PPG: subtract mean null trial
    ppgMeanNull(iSub,:) = squeeze(nanmean(subNullPpg{iSub},2));%mean across trials
    subPpgMeanNull = repmat(ppgMeanNull(iSub,:),[size(subStimPpg{iSub},2),1]);%repeat X stim trials
    if subtractMean
        subStimPpg{iSub} = subStimPpg{iSub} - subPpgMeanNull';
    end
    
    %Respiration: subtract mean null trial
    respMeanNull(iSub,:) = squeeze(nanmean(subNullResp{iSub},2));%mean across trials
    subRespMeanNull = repmat(respMeanNull(iSub,:),[size(subStimResp{iSub},2),1]);%repeat X stim trials
    if subtractMean
        subStimResp{iSub} = subStimResp{iSub} - subRespMeanNull';
    end
    
    
    %get mean pupil for high and low contrast, and for frequencies
    
    stimContrastTrials{iSub} = [contrastTrials{goodSubs(iSub),1}(stimTrials{iSub,1}==1) contrastTrials{goodSubs(iSub),2}(stimTrials{iSub,2}==1)];
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
    for iRoi=1:length(controlROIs)
        roiTseriesL = [roiTC{goodSubs(iSub),controlROIs(iRoi),1}.tSeries roiTC{goodSubs(iSub),controlROIs(iRoi),2}.tSeries];
        roiTseriesR = [roiTC{goodSubs(iSub),controlROIs(iRoi)+1,1}.tSeries roiTC{goodSubs(iSub),controlROIs(iRoi)+1,2}.tSeries];
        roiTseries = [roiTseriesL; roiTseriesR];
        
        for rwd=1:2
            controlTseriesRwdL = roiTC{goodSubs(iSub),controlROIs(iRoi),rwd}.tSeries;
            controlTseriesRwdR = roiTC{goodSubs(iSub),controlROIs(iRoi),rwd}.tSeries;
            controlTseriesRwd = [controlTseriesRwdL; controlTseriesRwdR];
            controlNullTseriesRwd{iSub,iRoi,rwd} = controlTseriesRwd(:,nullTrialsTRs{goodSubs(iSub),rwd}==1 & correctTRs{iSub,rwd});
            controlStimTseriesRwd{iSub,iRoi,rwd} = controlTseriesRwd(:,nullTrialsTRs{goodSubs(iSub),rwd}==0 & correctTRs{iSub,rwd});
        end
        
        %concatenate across rwd, average across voxels
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
        
        %         controlStimStd{iSub,iRoi} = squeeze(std(controlStimTrialTseries{iSub,iRoi},0,1));
        %         f = fft(controlStimTrialTseries{iSub,iRoi});
        %         controlStimAmp{iSub,iRoi} = squeeze(abs(f(2,:)));
        controlMeanStimTrial{iSub,iRoi} = squeeze(nanmean(controlStimTrialTseries{iSub,iRoi},2));%T
        ntrials = size(controlStimTrialTseries{iSub,iRoi},2);
        A = controlMeanStimTrial{iSub,iRoi};
        A = zscore(A);
        B = controlStimTrialTseries{iSub,iRoi};
        controlStimBetas{iSub,iRoi} = A\B;%A'x=B
        controlStimBetasMedian(iSub) = median(controlStimBetas{iSub,1});
        
        controlMeanNullTrial{iSub,iRoi} = squeeze(nanmean(controlNullTrialTseries{iSub,iRoi},2));%T
        ntrials = size(controlNullTrialTseries{iSub,iRoi},2);
        A = controlMeanNullTrial{iSub,iRoi};
        A = zscore(A);
        B = controlNullTrialTseries{iSub,iRoi};
        controlNullBetas{iSub,iRoi} = A\B;%A'x=B
        controlNullBetasMedian(iSub) = median(controlNullBetas{iSub,1});
    end
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
            %STD
            subBinStdNull{iSub,iRoi,ibin} = squeeze(std(binNullTrialTseries{iSub,iRoi,ibin},0,2));%vox, trial
            subBinStdStim{iSub,iRoi,ibin} = squeeze(std(binStimTrialTseries{iSub,iRoi,ibin},0,2));%vox, trial
            %FFT amplitude
            %             f = fft(binStimTrialTseries{iSub,iRoi,ibin},[],2);
            %             subBinAmpStim{iSub,iRoi,ibin} = squeeze(abs(f(:,2,:)));
            
            %compute mean stim trial per voxel
            binMeanStimTrial{iSub,iRoi,ibin} = squeeze(nanmean(binStimTrialTseries{iSub,iRoi,ibin},3));%vox, T
            %regress mean trial
            ntrials = size(binStimTrialTseries{iSub,iRoi,ibin},3);
            subBinBetasStim{iSub,iRoi,ibin}=[];
            for ivox=1:numBinVoxels(iSub,iRoi,ibin)
                
                A = zscore(binMeanStimTrial{iSub,iRoi,ibin}(ivox,:));
                B = squeeze(binStimTrialTseries{iSub,iRoi,ibin}(ivox,:,:));
                subBinBetasStim{iSub,iRoi,ibin}(ivox,:) = A'\B;%A'x=B
                %                 Y = squeeze(binStimTrialTseries{iSub,iRoi,ibin}(ivox,:,:));
                %                 b = regress(Y,X);%Y=X*b
            end
            
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
clear classAccFreq classAccContrast permAccFreq permAccContrast
for arousalType=arousalTypes
    
    for iSub=1:length(goodSubs)
        
        switch arousalType
            case 1
                trialArousal = subBaseStimPupil{iSub}';
                medianArousal = subStimPupilBaseMedian(iSub);
                trialArousalNull = subBaseNullPupil{iSub}';
                medianArousalNull = subNullPupilBaseMedian(iSub);
            case 2
                trialArousal = subStdStimPupil{iSub}';
                medianArousal = subStimPupilStdMedian(iSub);
                trialArousalNull = subStdNullPupil{iSub}';
                medianArousalNull = subNullPupilStdMedian(iSub);
            case 3
                trialArousal = subStdStimPulse{iSub}';
                medianArousal = subStimPulseStdMedian(iSub);
                trialArousalNull = subStdNullPulse{iSub}';
                medianArousalNull = subNullPulseStdMedian(iSub);
            case 4
                trialArousal = -subStimRwd{iSub}';
                medianArousal = -subStimRwdMedian(iSub);
                trialArousalNull = -subNullRwd{iSub}';
                medianArousalNull = -subNullRwdMedian(iSub);
            case 5
                trialArousal = controlStimBetas{iSub,1}';
                medianArousal = controlStimBetasMedian(iSub);
                trialArousalNull = controlNullBetas{iSub,1}';
                medianArousalNull = controlNullBetasMedian(iSub);
        end
        
        for arousal=1:2
            if arousal==1
                arousalTrials{arousalType,iSub,arousal} = trialArousal>medianArousal;
                arousalTrialsNull{arousalType,iSub,arousal} = trialArousalNull>medianArousalNull;
            else
                arousalTrials{arousalType,iSub,arousal} = trialArousal<medianArousal;
                arousalTrialsNull{arousalType,iSub,arousal} = trialArousalNull<medianArousalNull;
            end
        end
        
        %get mean bin BOLD timeseries for each arousal
        %separately for stim and for null
        for arousal=1:2
            for ibin=1:nbins
                arousalBinStimTrial(arousalType,iSub,ibin,arousal,:) = squeeze(nanmean(nanmean(binStimTrialTseries{iSub,iRoi,ibin}(:,:,arousalTrials{arousalType,iSub,arousal}),3)));%vox, T
                arousalBinNullTrial(arousalType,iSub,ibin,arousal,:) = squeeze(nanmean(nanmean(binNullTrialTseries{iSub,iRoi,ibin}(:,:,arousalTrialsNull{arousalType,iSub,arousal}),3)));%vox, T
            end
            arousalStimPulseTrial(arousalType,iSub,arousal,:) = squeeze(nanmean(subStimPulse{iSub}(:,arousalTrials{arousalType,iSub,arousal}),2));
            arousalNullPulseTrial(arousalType,iSub,arousal,:) = squeeze(nanmean(subNullPulse{iSub}(:,arousalTrialsNull{arousalType,iSub,arousal}),2));
            arousalStimPpgTrial(arousalType,iSub,arousal,:) = squeeze(nanmean(subStimPpg{iSub}(:,arousalTrials{arousalType,iSub,arousal}),2));
            arousalNullPpgTrial(arousalType,iSub,arousal,:) = squeeze(nanmean(subNullPpg{iSub}(:,arousalTrialsNull{arousalType,iSub,arousal}),2));
            arousalStimRespTrial(arousalType,iSub,arousal,:) = squeeze(nanmean(subStimResp{iSub}(:,arousalTrials{arousalType,iSub,arousal}),2));
            arousalNullRespTrial(arousalType,iSub,arousal,:) = squeeze(nanmean(subNullResp{iSub}(:,arousalTrialsNull{arousalType,iSub,arousal}),2));

            arousalStimPupilTrial(arousalType,iSub,arousal,:) = squeeze(nanmean(stimPupilConcat{iSub}(arousalTrials{arousalType,iSub,arousal},:)));
            arousalNullPupilTrial(arousalType,iSub,arousal,:) = squeeze(nanmean(nullPupilConcat{iSub}(arousalTrialsNull{arousalType,iSub,arousal},:)));

        end
    end
end
%%
ifig=0;

for arousalType=arousalTypes
    
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

    
    ifig=ifig+1;
    rows=3;
    cols=nbins;
    figure(ifig)
    minY = min(arousalBinStimTrial(:));
    maxY = max(arousalBinStimTrial(:));
    for ibin=1:nbins
        subplot(rows,cols,ibin)
        for arousal=1:2
            plot(squeeze(nanmean(arousalBinStimTrial(arousalType,:,ibin,arousal,:))),'color',plotColors{arousal});
            hold on
        end
        if ibin==1
            title([arousalStr ': stim'])
            ylabel('std fMRI signal');
        end
        xlim([0 trialLength+1]);
        ylim([minY maxY]);
        if ibin>1
            yticks([]);
        end
    end
    for ibin=1:nbins
        subplot(rows,cols,ibin+cols)
        for arousal=1:2
            plot(squeeze(nanmean(arousalBinNullTrial(arousalType,:,ibin,arousal,:))),'color',plotColors{arousal});
            hold on
        end
        if ibin==1
            title('null')
            ylabel('std fMRI signal');
        end
        xlim([0 trialLength+1]);
        ylim([minY maxY]);
        if ibin>1
            yticks([]);
        end
    end
    for ibin=1:nbins
        subplot(rows,cols,ibin+2*cols)
        for arousal=1:2
            plot(squeeze(nanmean(arousalBinStimTrial(arousalType,:,ibin,arousal,:)-arousalBinNullTrial(arousalType,:,ibin,arousal,:))),'color',plotColors{arousal});
            hold on
        end
        if ibin==1
            title('stim-null')
            ylabel('std fMRI signal');
        end
        xlim([0 trialLength+1]);
        ylim([minY maxY]);
        if ibin>1
            yticks([]);
        end
        xlabel('time (TR)');
    end
    
    ifig=ifig+1;
    figure(ifig)
    rows=3;
    cols=4;
    subplot(rows,cols,1)
     for arousal=1:2
        plot(squeeze(nanmean(arousalStimPulseTrial(arousalType,:,arousal,:))),'color',plotColors{arousal});
        hold on
     end
     title([arousalStr ': stim pulse']);
     subplot(rows,cols,2)
     for arousal=1:2
        plot(squeeze(nanmean(arousalStimPpgTrial(arousalType,:,arousal,:))),'color',plotColors{arousal});
        hold on
     end
     title('stim ppg');
          subplot(rows,cols,3)
     for arousal=1:2
        plot(squeeze(nanmean(arousalStimRespTrial(arousalType,:,arousal,:))),'color',plotColors{arousal});
        hold on
     end
     title('stim respiration');
     subplot(rows,cols,4)
     for arousal=1:2
        plot(squeeze(nanmean(arousalStimPupilTrial(arousalType,:,arousal,:))),'color',plotColors{arousal});
        hold on
     end
     title('stim pupil');
     
     subplot(rows,cols,1+cols)
     for arousal=1:2
        plot(squeeze(nanmean(arousalNullPulseTrial(arousalType,:,arousal,:))),'color',plotColors{arousal});
        hold on
     end
     title('null pulse');
     subplot(rows,cols,2+cols)
     for arousal=1:2
        plot(squeeze(nanmean(arousalNullPpgTrial(arousalType,:,arousal,:))),'color',plotColors{arousal});
        hold on
     end
     title('null ppg');
          subplot(rows,cols,3+cols)
     for arousal=1:2
        plot(squeeze(nanmean(arousalNullRespTrial(arousalType,:,arousal,:))),'color',plotColors{arousal});
        hold on
     end
     title('null respiration');
     subplot(rows,cols,4+cols)
     for arousal=1:2
        plot(squeeze(nanmean(arousalNullPupilTrial(arousalType,:,arousal,:))),'color',plotColors{arousal});
        hold on
     end
     title('null pupil');
    
    
    %     %%
    %
    %     ifig=ifig+1; figure(ifig); clf
    %     rows=1;
    %     cols=4;
    %     subplot(rows,cols,1:2)
    %     pvalThresh=0.05;
    %     markerSize=80;
    %
    %     % legend('std','amp','regress');
    %     legend('high','low');
    %
    %     old = {' ', ','};
    %     new = '_';
    %
    %     subplot(rows,cols,3)
    %     for arousal=1:2
    %         bar((arousal-1)*nfreqs+[1:5],squeeze(nanmean(freqCondsTrue(arousalType,featureType,:,arousal,:))),'facecolor',plotColors{arousal},'linestyle',lineStyles{arousal});
    %         hold on
    %         er=errorbar((arousal-1)*nfreqs+[1:nfreqs],squeeze(nanmean(freqCondsTrue(arousalType,featureType,:,arousal,:))),squeeze(nanstd(freqCondsTrue(arousalType,featureType,:,arousal,:))));
    %         er.Color = [0 0 0];
    %         er.LineStyle = 'none';
    %         %         histogram(squeeze(mean(freqCondsTrue(:,arousal,:))),nfreqs+1,'facecolor',plotColors{arousal},'linestyle',lineStyles{arousal});
    %
    %     end
    %     subplot(rows,cols,4)
    %     for arousal=1:2
    %         bar((arousal-1)*ncontrasts+[1:ncontrasts],squeeze(nanmean(contrastCondsTrue(arousalType,featureType,:,arousal,:))),'facecolor',plotColors{arousal},'linestyle',lineStyles{arousal});
    %         hold on
    %         er=errorbar((arousal-1)*ncontrasts+[1:ncontrasts],squeeze(nanmean(contrastCondsTrue(arousalType,featureType,:,arousal,:))),squeeze(nanstd(contrastCondsTrue(arousalType,featureType,:,arousal,:))));
    %         er.Color = [0 0 0];
    %         er.LineStyle = 'none';
    %         %         histogram(squeeze(mean(contrastCondsTrue(:,arousal,:))),ncontrasts+1,'facecolor',plotColors{arousal},'linestyle',lineStyles{arousal});
    %
    %     end
    set(gcf,'position',[100 100 1100 400])
end
% savefig(replace(titleStr,old, new));


%%
figure
rows=1;
cols=4;
subplot(rows,cols,1)
plot(squeeze(nanmean(meanContrastPupil))')
subplot(rows,cols,2)
plot(squeeze(nanmean(meanContrastPulse))')
subplot(rows,cols,3)
plot(squeeze(nanmean(meanFreqPupil))')
subplot(rows,cols,4)
plot(squeeze(nanmean(meanFreqPulse))')

toc
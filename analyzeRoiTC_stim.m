close all
mrQuit
clear all
bensonROIs = 1%[1:3];
eccMin = 0.2;
eccMax = 70;
nbins = 12;
saveFigs=1;
tic
toZscore=1;%0 or 1
concatProj= 0;
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
    'nbins','binBorders','binCenters',...
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
nullPupil = origNullPupil;
stimPupil = origStimPupil;
for iSub=1:length(goodSubs)
    
    for rwd=1:2
        %remove the first trial of each run
        junkedIndices = zeros(size(eyeStimTrials{goodSubs(iSub),rwd}));
        junkedIndices(1:17:end) = ones;
        stimTrialIndices = find(eyeStimTrials{goodSubs(iSub),rwd}==1);
        junkedStimTrials = find(mod(stimTrialIndices,17)==1);
        stimPupil{goodSubs(iSub),rwd}(junkedStimTrials,:) = [];
        
        junkedIndices = zeros(size(eyeNullTrials{goodSubs(iSub),rwd}));
        junkedIndices(1:17:end) = ones;
        nullTrialIndices = find(eyeNullTrials{goodSubs(iSub),rwd}==1);
        junkedNullTrials = find(mod(nullTrialIndices,17)==1);
        nullPupil{goodSubs(iSub),rwd}(junkedNullTrials,:) = [];
        
        %get measure of pupil size per trial
        stdStimPupil{iSub,rwd} = nanstd(stimPupil{goodSubs(iSub),rwd},0,2);%per trial
        baseStimPupil{iSub,rwd} = nanmean(stimPupil{goodSubs(iSub),rwd}(:,1:baselineT),2);%per trial
        f=fft(fillmissing(stimPupil{goodSubs(iSub),rwd},'linear',2),[],2);
        phStimPupil{iSub,rwd} = angle(f(:,2));
        ampStimPupil{iSub,rwd} = abs(f(:,2));
        
        stdNullPupil{iSub,rwd} = nanstd(nullPupil{goodSubs(iSub),rwd},0,2);%per trial
        baseNullPupil{iSub,rwd} = nanmean(nullPupil{goodSubs(iSub),rwd}(:,1:baselineT),2);%per trial
        f=fft(fillmissing(nullPupil{goodSubs(iSub),rwd},'linear',2),[],2);
        phNullPupil{iSub,rwd} = angle(f(:,2));
        ampNullPupil{iSub,rwd} = abs(f(:,2));
        
        
        
        
        
        stimFreqTrials = freqTrials{goodSubs(iSub),rwd}(stimTrials{goodSubs(iSub),rwd}==1);
        for iRoi=1:length(roiNames)
            for ifreq=1:max(freqTrials{goodSubs(iSub),rwd})
                subRoiRwdFreqTseries{iSub,iRoi,rwd,ifreq} = squeeze(nanmean(stimTrialTseries{goodSubs(iSub),iRoi,rwd}(:,:,stimFreqTrials==ifreq)));%mean across voxels
            end
            
            %get measure of fMRI response amplitude
            roiStimTrialTseries{iSub,iRoi,rwd} = squeeze(nanmean(stimTrialTseries{goodSubs(iSub),iRoi,rwd}));%averaged over voxels
            rwdRoiStdStim{iSub,iRoi,rwd} = std(roiStimTrialTseries{iSub,iRoi,rwd});
            f = fft(roiStimTrialTseries{iSub,iRoi,rwd});
            rwdRoiAmpStim{iSub,iRoi,rwd} = abs(f(2,:));
            rwdRoiPhStim{iSub,iRoi,rwd} = angle(f(2,:));
            
            roiNullTrialTseries{iSub,iRoi,rwd} = squeeze(nanmean(nullTrialTseries{goodSubs(iSub),iRoi,rwd}));%averaged over voxels
            rwdRoiStdNull{iSub,iRoi,rwd} = std(roiNullTrialTseries{iSub,iRoi,rwd});
            f = fft(roiNullTrialTseries{iSub,iRoi,rwd});
            rwdRoiAmpNull{iSub,iRoi,rwd} = abs(f(2,:));
            rwdRoiPhNull{iSub,iRoi,rwd} = angle(f(2,:));
        end
    end
    
    %concatenate across rwd, separately for each subject
    subStdStimPupil{iSub} = [stdStimPupil{iSub,1}; stdStimPupil{iSub,2}];
    subBaseStimPupil{iSub} = [baseStimPupil{iSub,1}; baseStimPupil{iSub,2}];
    subPhStimPupil{iSub} = [phStimPupil{iSub,1}; phStimPupil{iSub,2}];
    subAmpStimPupil{iSub} = [ampStimPupil{iSub,1}; ampStimPupil{iSub,2}];
    subStimPupilStdMedian(iSub) = nanmedian(subStdStimPupil{iSub});
    subStimPupilBaseMedian(iSub) = nanmedian(subBaseStimPupil{iSub});
    
    subStdNullPupil{iSub} = [stdNullPupil{iSub,1}; stdNullPupil{iSub,2}];
    subBaseNullPupil{iSub} = [baseNullPupil{iSub,1}; baseNullPupil{iSub,2}];
    subPhNullPupil{iSub} = [phNullPupil{iSub,1}; phNullPupil{iSub,2}];
    subAmpNullPupil{iSub} = [ampNullPupil{iSub,1}; ampNullPupil{iSub,2}];
    
    %pulse stats
    subStimPulse{iSub} = [rwdPulseTrials{goodSubs(iSub),1}(:,stimTrials{goodSubs(iSub),1}==1) rwdPulseTrials{goodSubs(iSub),2}(:,stimTrials{goodSubs(iSub),2}==1)];%concatenate pulse of null trials
    subStimPulseStd{iSub} = std(subStimPulse{iSub});
    subStimPulseMean{iSub} = mean(subStimPulse{iSub});
    %             f=fft(fillmissing(subPulse{iSub},'linear',2),[],2);
    f=fft(subStimPulse{iSub});
    subStimPulsePh{iSub} = angle(f(2,:));
    subStimPulseAmp{iSub} = abs(f(2,:));
    
    
    subNullPulse{iSub} = [rwdPulseTrials{goodSubs(iSub),1}(:,nullTrials{goodSubs(iSub),1}==1) rwdPulseTrials{goodSubs(iSub),2}(:,nullTrials{goodSubs(iSub),2}==1)];%concatenate pulse of null trials
    subNullPulseStd{iSub} = std(subNullPulse{iSub});
    subNullPulseMean{iSub} = mean(subNullPulse{iSub});
    %             f=fft(fillmissing(subPulse{iSub},'linear',2),[],2);
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
    stimPupilConcat{iSub} = [stimPupil{goodSubs(iSub),1}(:,1:7500); stimPupil{goodSubs(iSub),2}(:,1:7500)];
    
    stimContrastTrials = [contrastTrials{goodSubs(iSub),1}(stimTrials{goodSubs(iSub),1}==1) contrastTrials{goodSubs(iSub),2}(stimTrials{goodSubs(iSub),2}==1)];
    for icontrast=1:max(stimContrastTrials)
        contrastPupil{iSub,icontrast} = stimPupilConcat{iSub}(stimContrastTrials==icontrast,:);
        meanContrastPupil(iSub,icontrast,:) = nanmean(contrastPupil{iSub,icontrast});
    end
    
    stimFreqTrials = [freqTrials{goodSubs(iSub),1}(stimTrials{goodSubs(iSub),1}==1) freqTrials{goodSubs(iSub),2}(stimTrials{goodSubs(iSub),2}==1)];
    for ifreq=1:max(stimFreqTrials)
        freqPupil{iSub,ifreq} = stimPupilConcat{iSub}(stimFreqTrials==ifreq,:);
        meanFreqPupil(iSub,ifreq,:) = nanmean(freqPupil{iSub,ifreq});
    end
    
    
    %concatenate physio design matrcies
    pulseDesignMat{iSub} = [designMatPulse{iSub,1} designMatPulse{iSub,2}];
    respDesignMat{iSub} = [designMatResp{iSub,1} designMatResp{iSub,2}];
    respPulseDesignMat{iSub} = [designMatRespPulse{iSub,1} designMatRespPulse{iSub,2}] ;
    
    for iRoi=1:length(roiNames)
        %concatenate across rwd
        for ifreq=1:size(subRoiRwdFreqTseries,4)
            subRoiFreqTseries{iSub,iRoi,ifreq} = [subRoiRwdFreqTseries{iSub,iRoi,1,ifreq} subRoiRwdFreqTseries{iSub,iRoi,2,ifreq}];
            meanFreqTrial(iSub,iRoi,ifreq,:) = mean(subRoiFreqTseries{iSub,iRoi,ifreq},2);
            subRoiFreqStd{iSub,iRoi,ifreq} =std(subRoiFreqTseries{iSub,iRoi,ifreq});
        end
        
        subStimTrialTseries{iSub,iRoi} = [roiStimTrialTseries{iSub,iRoi,1} roiStimTrialTseries{iSub,iRoi,2}];
        subRoiStdStim{iSub,iRoi} = [rwdRoiStdStim{iSub,iRoi,1}'; rwdRoiStdStim{iSub,iRoi,2}'];
        subRoiAmpStim{iSub,iRoi} = [rwdRoiAmpStim{iSub,iRoi,1}'; rwdRoiAmpStim{iSub,iRoi,2}'];
        subRoiPhStim{iSub,iRoi} = [rwdRoiPhStim{iSub,iRoi,1}'; rwdRoiPhStim{iSub,iRoi,2}'];
        
        subNullTrialTseries{iSub,iRoi} = [roiNullTrialTseries{iSub,iRoi,1} roiNullTrialTseries{iSub,iRoi,2}];
        subRoiStdNull{iSub,iRoi} = [rwdRoiStdNull{iSub,iRoi,1}'; rwdRoiStdNull{iSub,iRoi,2}'];
        subRoiAmpNull{iSub,iRoi} = [rwdRoiAmpNull{iSub,iRoi,1}'; rwdRoiAmpNull{iSub,iRoi,2}'];
        subRoiPhNull{iSub,iRoi} = [rwdRoiPhNull{iSub,iRoi,1}'; rwdRoiPhNull{iSub,iRoi,2}'];
        
        %concatenate across rwd for physio regression
        %         subTseries{iSub,iRoi} = nanmean([roiTC{iSub,iRoi,1}.tSeries roiTC{iSub,iRoi,2}.tSeries]);%mean over voxels
        
        %         pulseKernelRoi(iSub,iRoi,:) = pulseDesignMat{iSub}'\subTseries{iSub,iRoi}';
        %         respKernelRoi(iSub,iRoi,:) = respDesignMat{iSub}'\subTseries{iSub,iRoi}';
        %         respPulseKernelRoi(iSub,iRoi,:) = respPulseDesignMat{iSub}'\subTseries{iSub,iRoi}';
        
    end
    
    iRoi=1;
    %     pulsePupilBensonNull{iSub} = [subNullRespStd{iSub}', subNullRespMean{iSub}',subNullPulseStd{iSub}', subNullPulseMean{iSub}', ...
    %         subStdNullPupil{iSub}, subAmpNullPupil{iSub}, subBaseNullPupil{iSub}, subPhNullPupil{iSub}, ...
    %         subRoiAmpNull{iSub,iRoi}, subRoiStdNull{iSub,iRoi}, subRoiPhNull{iSub,iRoi}];
    %     [corrPulsePupilBensonNull(iSub,:,:) pvalPulsePupilBensonNull(iSub,:,:)] = corr(pulsePupilBensonNull{iSub},'rows','complete');
    
end
meanFreqTrialStd = squeeze(std(meanFreqTrial,0,4));
%%

binBorders = logspace(log10(eccMin),log10(eccMax),nbins+1);
nbins = length(binBorders)-1;
for i=2:length(binBorders)
    binCenters(i-1) = (binBorders(i)+binBorders(i-1))/2;
end

for iSub=1:length(goodSubs)
    for iRoi=bensonROIs%1:length(roiNames) %length(bensonROIs)
        for ibin=1:nbins
            for rwd=1:2
                binVoxels = eccen{goodSubs(iSub),iRoi}>binBorders(ibin) & eccen{goodSubs(iSub),iRoi}<=binBorders(ibin+1);
                numBinVoxels(iSub,iRoi,ibin) = sum(binVoxels);
                %                                 binVoxels = binVoxels & areas{goodSubs(iSub),iRoi}==1;%ONLY V1
                binVoxels = binVoxels & areas{goodSubs(iSub),iRoi}>0;%ONLY V1,V2,V3
                binMeanTseries = nanmean(roiTC{goodSubs(iSub),iRoi,rwd}.tSeries(binVoxels,:),1);%mean timecourse across voxels
                binNullTseries{iSub,iRoi,ibin,rwd} = binMeanTseries(nullTrialsTRs{goodSubs(iSub),rwd}==1);
                binNullTrialTseries{iSub,iRoi,ibin,rwd} = reshape(binNullTseries{iSub,iRoi,ibin,rwd},trialLength,[]);%iSub,T,trial
                %get fMRI response amplitude
                rwdBinStdNull{iSub,iRoi,ibin,rwd} = std(binNullTrialTseries{iSub,iRoi,ibin,rwd});%
                f = fft(binNullTrialTseries{iSub,iRoi,ibin,rwd});
                rwdBinAmpNull{iSub,iRoi,ibin,rwd} = abs(f(2,:));
                rwdBinPhNull{iSub,iRoi,ibin,rwd} = angle(f(2,:));
                %stim
                binStimTseries{iSub,iRoi,ibin,rwd} = binMeanTseries(nullTrialsTRs{goodSubs(iSub),rwd}==0);
                binStimTrialTseries{iSub,iRoi,ibin,rwd} = reshape(binStimTseries{iSub,iRoi,ibin,rwd},trialLength,[]);%iSub,T,trial
                
                %frequency
                %                 stimFreqTrials = freqTrials{goodSubs(iSub),rwd}(stimTrials{iSub,rwd}==1);
                %                 binMeanTrials = reshape(binMeanTseries,trialLength,[]);
                %                 for ifreq=1:max(freqTrials{goodSubs(iSub),rwd})
                %                     subBinFreqTseries{iSub,iRoi,ibin,ifreq} = squeeze(nanmean(stimTrialTseries{goodSubs(iSub),iRoi,rwd}(:,:,stimFreqTrials==ifreq)));%mean across voxels
                %                     subBinFreqMeanTrial(iSub,iRoi,ibin,ifreq,:)
                %                 end
            end
            %concatenate across rwd
            subBinNullTrialTseries{iSub,iRoi,ibin} = [binNullTrialTseries{iSub,iRoi,ibin,1} binNullTrialTseries{iSub,iRoi,ibin,2}];
            subBinStdNull{iSub,iRoi,ibin} = [rwdBinStdNull{iSub,iRoi,ibin,1}'; rwdBinStdNull{iSub,iRoi,ibin,2}'];
            %             subBinAmpNull{iSub,iRoi,ibin} = [rwdBinAmpNull{iSub,iRoi,ibin,1}'; rwdBinAmpNull{iSub,iRoi,ibin,2}'];
            %             subBinPhNull{iSub,iRoi,ibin} = [rwdBinPhNull{iSub,iRoi,ibin,1}'; rwdBinPhNull{iSub,iRoi,ibin,2}'];
            
            
            subBinStimTrialTseries{iSub,iRoi,ibin} = [binStimTrialTseries{iSub,iRoi,ibin,1} binStimTrialTseries{iSub,iRoi,ibin,2}];
            subStimTrialRwd(iSub,iRoi,1:size(binStimTrialTseries{iSub,iRoi,ibin,1},2)) = ones;
            subStimTrialRwd(iSub,iRoi,size(binStimTrialTseries{iSub,iRoi,ibin,1},2)+1:end) = 2*ones;
            %             [size(binStimTrialTseries{iSub,iRoi,ibin,1},2), size(binStimTrialTseries{iSub,iRoi,ibin,2},2)];
            
            %frequency
            stimFreqTrials = [freqTrials{goodSubs(iSub),1}(stimTrials{goodSubs(iSub),1}==1) freqTrials{goodSubs(iSub),2}(stimTrials{goodSubs(iSub),2}==1)];
            for ifreq=1:max(freqTrials{goodSubs(iSub),rwd})
                subBinFreqTseries{iSub,iRoi,ibin,ifreq} = subBinStimTrialTseries{iSub,iRoi,ibin}(:,stimFreqTrials==ifreq);
                subBinFreqMeanTrial(iSub,iRoi,ibin,ifreq,:) = mean(subBinFreqTseries{iSub,iRoi,ibin,ifreq},2);
                for rwd=1:2
                    subBinFreqTseriesRwd = binStimTrialTseries{iSub,iRoi,ibin,rwd}(:,freqTrials{goodSubs(iSub),rwd}(stimTrials{goodSubs(iSub),rwd}==1)==ifreq);
                    %                     subBinFreqTseriesRwd = binStimTrialTseries{iSub,iRoi,ibin,rwd}(:,stimTrials{goodSubs(iSub),rwd}==1);
                    subBinFreqMeanTrialRwd(iSub,iRoi,ibin,rwd,ifreq,:) = mean(subBinFreqTseriesRwd,2);
                end
                for arousal=1:2
                    if arousal==1
%                         arousalTrials = subStdStimPupil{iSub}'>subStimPupilStdMedian(iSub);
                        arousalTrials = subBaseStimPupil{iSub}'>subStimPupilBaseMedian(iSub);
                    else
%                         arousalTrials = subStdStimPupil{iSub}'<subStimPupilStdMedian(iSub);
                        arousalTrials = subBaseStimPupil{iSub}'<subStimPupilBaseMedian(iSub);
                    end
                    subBinFreqTseriesArousal = subBinStimTrialTseries{iSub,iRoi,ibin}(:,stimFreqTrials==ifreq & arousalTrials);
                    subBinFreqMeanTrialArousal(iSub,iRoi,ibin,arousal,ifreq,:) = mean(subBinFreqTseriesArousal,2);
                end
            end
            %contrast
            for icontrast=1:max(contrastTrials{goodSubs(iSub),rwd})
                stimContrastTrials =  [contrastTrials{goodSubs(iSub),1}(stimTrials{goodSubs(iSub),1}==1) contrastTrials{goodSubs(iSub),2}(stimTrials{goodSubs(iSub),2}==1)];
                subBinContrastTseries{iSub,iRoi,ibin,icontrast} = subBinStimTrialTseries{iSub,iRoi,ibin}(:,stimContrastTrials==icontrast);
                subBinContrastMeanTrial(iSub,iRoi,ibin,icontrast,:) = mean(subBinContrastTseries{iSub,iRoi,ibin,icontrast},2);
            end
            
        end
        
    end
    
end


meanBinFreqStd = std(subBinFreqMeanTrial,0,5);
meanBinContrastStd = std(subBinContrastMeanTrial,0,5);

meanBinFreqStdRwd = std(subBinFreqMeanTrialRwd,0,6);
meanBinFreqStdArousal = std(subBinFreqMeanTrialArousal,0,6);


toc

%% FIT TUNING CURVES
%% Fit SF tuning curves
options = optimset('MaxFunEvals',500, 'MaxIter',500);%used by fminsearch
% initParams = [3, 1, 1, 0];`
initParams = [0, 10];

for iRoi = bensonROIs
    for ibin=1:nbins
        d = squeeze(mean(meanBinFreqStd(:,iRoi,ibin,:)));
        d = zscore(d);
        myfunc = @(x) gaussianRms(x, d, freqs');
        [x fval exflag] = fminsearch(myfunc, initParams, options);
        sfStdTuning(iRoi,ibin,:) = x;
        sfStdExitFlag(iRoi,ibin) = exflag;
        sfStdRms(iRoi,ibin) = fval;%1 - correlation
        for rwd=1:2
            d = squeeze(mean(meanBinFreqStdRwd(:,iRoi,ibin,rwd,:)));
            d = zscore(d);
            myfunc = @(x) gaussianRms(x, d, freqs');
            [x fval exflag] = fminsearch(myfunc, initParams, options);
            sfStdTuningRwd(iRoi,ibin,rwd,:) = x;
            sfStdExitFlagRwd(iRoi,ibin,rwd) = exflag;
            sfStdRmsRwd(iRoi,ibin,rwd) = fval;%1 - correlation
        end
        for arousal=1:2
            d = squeeze(mean(meanBinFreqStdArousal(:,iRoi,ibin,arousal,:)));
            d = zscore(d);
            myfunc = @(x) gaussianRms(x, d, freqs');
            [x fval exflag] = fminsearch(myfunc, initParams, options);
            sfStdTuningArousal(iRoi,ibin,arousal,:) = x;
            sfStdExitFlagArousal(iRoi,ibin,arousal) = exflag;
            sfStdRmsArousal(iRoi,ibin,arousal) = fval;%1 - correlation
        end
    end
end

% for ibin=1:nbins
%     for rwd=1:2
%         respStd = squeeze(mean(freqDeconvStd(:,ibin,rwd,:)));%responses
%         myfunc = @(x) gaussianRms(x, respStd, (1:numFreqs)');
%         [x fval exflag] = fminsearch(myfunc, initParams, options);
%         x(3) = abs(x(3));
%         sfStdTuning(ibin,rwd,:) = x;
%         sfStdExitFlag(ibin,rwd) = exflag;
%         sfStdRms(ibin,rwd) = fval;%fit RMS
%
%         for c=1:numContrasts
%             respStd = squeeze(mean(contrastFreqDeconvStd(:,ibin,rwd,(c-1)*numFreqs+1:c*numFreqs)));%responses
%             myfunc = @(x) gaussianRms(x, respStd, (1:numFreqs)');
%             [x fval exflag] = fminsearch(myfunc, initParams, options);
%             x(3) = abs(x(3));
%             sfContrastStdTuning(ibin,rwd,c,:) = x;
%             sfContrastStdExitFlag(ibin,rwd,c) = exflag;
%             sfContrastStdRms(ibin,rwd,c) = fval;%fit RMS
%         end
%     end
% end





toc
%%
ifig=0;
plotColors = {[1 0 0], [0 0 1], [0 1 0], [1 1 0]};
rows=length(roiNames);
cols=6;
linestyle = '-';


%%
% respKernelRoi = zscore(respKernelRoi,0,3);
% pulseKernelRoi = zscore(pulseKernelRoi,0,3);
% respPulseKernelRoi = zscore(respPulseKernelRoi,0,3);

%%
iRoi=1;
nfreqs=size(meanBinFreqStd,4);
ncontrasts =size(meanBinContrastStd,4);
ifig=ifig+1; figure(ifig); clf
plotcmap = cool;
rows=2; cols=3;
%frequency
subplot(rows,cols,1)
for ifreq=1:nfreqs
    plot(squeeze(mean(meanBinFreqStd(:,iRoi,:,ifreq))),'Color',plotcmap(ceil(size(plotcmap,1)*ifreq/nfreqs),:))
    hold on
end
xlabel('eccentricity bin');
subplot(rows,cols,2)
for ibin=1:nbins
    plot(squeeze(mean(meanBinFreqStd(:,iRoi,ibin,:))),'Color',plotcmap(ceil(size(plotcmap,1)*ibin/nbins),:))
    hold on
end
xlabel('frequency')

subplot(rows,cols,3)
for ibin=1:nbins
    plot(zscore(squeeze(mean(meanBinFreqStd(:,iRoi,ibin,:)))),'Color',plotcmap(ceil(size(plotcmap,1)*ibin/nbins),:))
    hold on
end
xlabel('frequency')
colormap cool
colorbar

%contrast
subplot(rows,cols,cols+1)
for icontrast=1:ncontrasts
    plot(squeeze(mean(meanBinContrastStd(:,iRoi,:,icontrast))),'Color',plotcmap(ceil(size(plotcmap,1)*icontrast/ncontrasts),:))
    hold on
end
xlabel('eccentricity bin');
subplot(rows,cols,cols+2)
for ibin=1:nbins
    plot(squeeze(mean(meanBinContrastStd(:,iRoi,ibin,:))),'Color',plotcmap(ceil(size(plotcmap,1)*ibin/nbins),:))
    hold on
end
xlabel('contrast')

subplot(rows,cols,cols+3)
for ibin=1:nbins
    plot(zscore(squeeze(mean(meanBinContrastStd(:,iRoi,ibin,:)))),'Color',plotcmap(ceil(size(plotcmap,1)*ibin/nbins),:))
    hold on
end
xlabel('contrast')

%% plot null and stim responses per bin

rows=length(bensonROIs);
cols=nbins;
for scaled=0:1
    scaledStr = '';
    if scaled
        scaledStr = '_scaled';
    end
    ifig=ifig+1; figure(ifig); clf
    for iRoi=bensonROIs
        for ibin=1:nbins
            subplot(rows,cols,ibin + (iRoi-1)*cols)
            plot(mean(subBinStimTrialTseries{iSub,iRoi,ibin},2),'m');
            hold on
            plot(mean(subBinNullTrialTseries{iSub,iRoi,ibin},2),'g');
            title([num2str(binBorders(ibin),'%.1f') '-' num2str(binBorders(ibin+1),'%.1f')]);
            if ibin==1
                legend('stim','null')
                title({roiNames{iRoi}; [num2str(binBorders(ibin),'%.1f') '-' num2str(binBorders(ibin+1),'%.1f')]});
            end
            if ~scaled
                ylim([-0.6 1]);
                if ibin>1
                    yticks([]);
                    yticklabels('');
                end
            end
        end
    end
    set(gcf,'position',[100 100 700 400])
    if saveFigs;  print('-painters','-dpdf',['~/Documents/MATLAB/rwdRapid/figures/rwdBinsNullVsStim' scaledStr '.pdf']); end;
end

%% SF tuning for ALL TRIALS
ifig=ifig+1; figure(ifig); clf
rows=length(bensonROIs);
cols=nbins+2;
for iRoi = bensonROIs
    for ibin=1:nbins
        subplot(rows,cols,ibin + (iRoi-1)*cols)
        d = squeeze(mean(meanBinFreqStd(:,iRoi,ibin,:)));
        zi = exp(-((log(freqs)-log(sfStdTuning(iRoi,ibin,1))).^2/2/sfStdTuning(iRoi,ibin,2)^2));
        plot(freqs,zscore(d),'k')
        hold on
        errorbar(freqs,zscore(d), squeeze(std(zscore(meanBinFreqStd(:,iRoi,ibin,:),0,4))),'k');
        plot(freqs,zscore(zi),'r');
        title(num2str(1-sfStdRms(iRoi,ibin),'%.2f'));
    end
    subplot(rows,cols,cols-1 + (iRoi-1)*cols)
    plot(sfStdTuning(iRoi,:,1));
    subplot(rows,cols,cols + (iRoi-1)*cols)
    plot(sfStdTuning(iRoi,:,2));
end
set(gcf,'position',[100 100 1000,500]);
%% SF tuning for RWD
ifig=ifig+1; figure(ifig); clf
rows=2*length(bensonROIs);
cols=nbins;
temp = squeeze(mean(meanBinFreqStdRwd));%mean over subjects
maxY= max(temp(:));
minY = 0;%min(temp(:));
for iRoi = bensonROIs
    for rwd=1:2
        for ibin=1:nbins
            subplot(rows,cols,ibin + (iRoi-1)*2*cols + (rwd-1)*cols)
            d = squeeze(mean(meanBinFreqStdRwd(:,iRoi,ibin,rwd,:)));
            zi = exp(-((log(freqs)-log(sfStdTuningRwd(iRoi,ibin,rwd,1))).^2/2/sfStdTuningRwd(iRoi,ibin,rwd,2)^2));
            plot(freqs,d,'k')
            hold on
%             errorbar(freqs,zscore(d), squeeze(std(zscore(meanBinFreqStdRwd(:,iRoi,ibin,rwd,:),0,5))),'k');
            plot(freqs,mean(d) + zscore(zi)*std(d),'r');
            hline(std(mean(subBinNullTrialTseries{iSub,iRoi,ibin},2)));
            title(num2str(1-sfStdRmsRwd(iRoi,ibin,rwd),'%.2f'));
            ylim([minY maxY]);
        end
%         subplot(rows,cols,cols-1 + (iRoi-1)*2*cols + (rwd-1)*cols)
%         plot(sfStdTuningRwd(iRoi,:,rwd,1));
%         subplot(rows,cols,cols + (iRoi-1)*2*cols + (rwd-1)*cols)
%         plot(sfStdTuningRwd(iRoi,:,rwd,2));
    end
end
set(gcf,'position',[100 100 1000,500]);

%% SF tuning for AROUSAL
ifig=ifig+1; figure(ifig); clf
rows=2*length(bensonROIs);
cols=nbins;%+2;
temp = squeeze(mean(meanBinFreqStdArousal));%mean over subjects
maxY= max(temp(:));
minY = 0;%min(temp(:));
for iRoi = bensonROIs
    for rwd=1:2
        for ibin=1:nbins
            subplot(rows,cols,ibin + (iRoi-1)*2*cols + (rwd-1)*cols)
            d = squeeze(mean(meanBinFreqStdArousal(:,iRoi,ibin,rwd,:)));
            zi = exp(-((log(freqs)-log(sfStdTuningArousal(iRoi,ibin,rwd,1))).^2/2/sfStdTuningArousal(iRoi,ibin,rwd,2)^2));
            plot(freqs,d,'k')
            hold on
            %             errorbar(freqs,d, squeeze(std(meanBinFreqStdArousal(:,iRoi,ibin,rwd,:))),'k');
            plot(freqs,mean(d) + zscore(zi)*std(d),'r');
            hline(std(mean(subBinNullTrialTseries{iSub,iRoi,ibin},2)));
            title(num2str(1-sfStdRmsArousal(iRoi,ibin,rwd),'%.2f'));
            ylim([minY maxY]);
        end
        %         subplot(rows,cols,cols-1 + (iRoi-1)*2*cols + (rwd-1)*cols)
        %         plot(sfStdTuningArousal(iRoi,:,rwd,1));
        %         subplot(rows,cols,cols + (iRoi-1)*2*cols + (rwd-1)*cols)
        %         plot(sfStdTuningArousal(iRoi,:,rwd,2));
    end
end
set(gcf,'position',[100 100 1000,500]);

%% MEAN pupil for contrasts and frequencies
ifig=ifig+1; figure(ifig); clf
rows=1;
cols=4;
subplot(rows,cols,1)
for icontrast=1:2
    plot(squeeze(meanContrastPupil(:,icontrast,:))','color',plotColors{icontrast})
    hold on
end
legend('high contrast');
title('single subjects');
subplot(rows,cols,2)
plot(squeeze(meanContrastPupil(:,1,:) - meanContrastPupil(:,2,:))')
hold all
title('high - low contrast');
subplot(rows,cols,3)
for icontrast=1:2
    plot(squeeze(mean(meanContrastPupil(:,icontrast,:)))','color',plotColors{icontrast})
    hold on
end
title('group average');
subplot(rows,cols,4)
plot(squeeze(mean(meanFreqPupil))')
title('spatial frequencies');
% plot(squeeze(meanFreqPupil),'g')
set(gcf,'position',[100 100 1000,300]);
%%
%%%%%
function [c] = gaussianRms(x, orig, xi)
% mean = x(1);
% sigma = x(2);
% amp = x(3);
% baseline = x(4);
% x(3) = abs(x(3));%so that we don't have inverted tuning curves
% zi = x(4) + x(3)*exp(-((xi-x(1)).^2/2/x(2)^2));
% c = sqrt(sum((orig-zi).^2));
% if x(1)<0; x(1) = 0; end;
zi = exp(-((log(xi)-log(x(1))).^2/2/x(2)^2));
c = 1-corr(zi,orig);


end
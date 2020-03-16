close all
mrQuit
clear all
tic
toZscore=1;%0 or 1
concatProj= 1;
% dataFolder = '/Volumes/MH02086153MACDT-Drobo/rwdFmri/';
dataFolder = '/Users/rothzn/rwdFmri/';

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
goodSubs = 1:2;
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
        

        for iRoi=1:length(roiNames)       
            stimFreqTrials = freqTrials{goodSubs(iSub),rwd}(stimTrials{iSub,rwd}==1);
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
        subTseries{iSub,iRoi} = nanmean([roiTC{iSub,iRoi,1}.tSeries roiTC{iSub,iRoi,2}.tSeries]);%mean over voxels

        pulseKernelRoi(iSub,iRoi,:) = pulseDesignMat{iSub}'\subTseries{iSub,iRoi}';
        respKernelRoi(iSub,iRoi,:) = respDesignMat{iSub}'\subTseries{iSub,iRoi}';
        respPulseKernelRoi(iSub,iRoi,:) = respPulseDesignMat{iSub}'\subTseries{iSub,iRoi}';

    end
    
    iRoi=1;
    pulsePupilBensonNull{iSub} = [subNullRespStd{iSub}', subNullRespMean{iSub}',subNullPulseStd{iSub}', subNullPulseMean{iSub}', ...
        subStdNullPupil{iSub}, subAmpNullPupil{iSub}, subBaseNullPupil{iSub}, subPhNullPupil{iSub}, ...
        subRoiAmpNull{iSub,iRoi}, subRoiStdNull{iSub,iRoi}, subRoiPhNull{iSub,iRoi}];
    [corrPulsePupilBensonNull(iSub,:,:) pvalPulsePupilBensonNull(iSub,:,:)] = corr(pulsePupilBensonNull{iSub},'rows','complete');
 
end
meanFreqTrialStd = squeeze(std(meanFreqTrial,0,4));
%%
eccMin = 0.2;
eccMax = 70;
nbins = 12;
binBorders = logspace(log10(eccMin),log10(eccMax),nbins+1);
nbins = length(binBorders)-1;
for i=2:length(binBorders)
    binCenters(i-1) = (binBorders(i)+binBorders(i-1))/2;
end
for iSub=1:length(goodSubs)
    for iRoi=1:length(roiNames)
        for ibin=1:nbins
            for rwd=1:2
                binVoxels = eccen{goodSubs(iSub),iRoi}>binBorders(ibin) & eccen{goodSubs(iSub),iRoi}<=binBorders(ibin+1);
                numBinVoxels(iSub,iRoi,ibin) = sum(binVoxels);
                %                 binVoxels = binVoxels & areas{iSub,iRoi}==1;%ONLY V1
                binVoxels = binVoxels & areas{goodSubs(iSub),iRoi}>0;%ONLY V1,V2,V3
                binMeanTseries = nanmean(roiTC{goodSubs(iSub),iRoi,rwd}.tSeries(binVoxels,:));%mean timecourse across voxels
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
            subBinAmpNull{iSub,iRoi,ibin} = [rwdBinAmpNull{iSub,iRoi,ibin,1}'; rwdBinAmpNull{iSub,iRoi,ibin,2}'];
            subBinPhNull{iSub,iRoi,ibin} = [rwdBinPhNull{iSub,iRoi,ibin,1}'; rwdBinPhNull{iSub,iRoi,ibin,2}'];
            
            meanResp = mean(subBinNullTrialTseries{iSub,iRoi,ibin},2);
            f=fft(meanResp);
            meanPh = angle(f(2));
            %convert phase values to be raletive to mean response
            subBinPhNull{iSub,iRoi,ibin} = mod(subBinPhNull{iSub,iRoi,ibin} - meanPh + pi, 2*pi)-pi;
            
            subBinStimTrialTseries{iSub,iRoi,ibin} = [binStimTrialTseries{iSub,iRoi,ibin,1} binStimTrialTseries{iSub,iRoi,ibin,2}];

            %frequency
            for ifreq=1:max(freqTrials{goodSubs(iSub),rwd})           
                stimFreqTrials = [freqTrials{goodSubs(iSub),1}(stimTrials{iSub,1}==1) freqTrials{goodSubs(iSub),2}(stimTrials{iSub,2}==1)];
                subBinFreqTseries{iSub,iRoi,ibin,ifreq} = subBinStimTrialTseries{iSub,iRoi,ibin}(:,stimFreqTrials==ifreq);
                subBinFreqMeanTrial(iSub,iRoi,ibin,ifreq,:) = mean(subBinFreqTseries{iSub,iRoi,ibin,ifreq},2);
            end
            %contrast
            for icontrast=1:max(contrastTrials{goodSubs(iSub),rwd})
                stimContrastTrials =  [contrastTrials{goodSubs(iSub),1}(stimTrials{iSub,1}==1) contrastTrials{goodSubs(iSub),2}(stimTrials{iSub,2}==1)];
                subBinContrastTseries{iSub,iRoi,ibin,icontrast} = subBinStimTrialTseries{iSub,iRoi,ibin}(:,stimContrastTrials==icontrast);
                subBinContrastMeanTrial(iSub,iRoi,ibin,icontrast,:) = mean(subBinContrastTseries{iSub,iRoi,ibin,icontrast},2);
            end
            
        end
        % model null trials as linear combination of control ROI amplitude
        ampPredictors = [subBinStdNull{iSub,iRoi,nbins} subBaseNullPupil{iSub} subStdNullPupil{iSub} subNullPulseStd{iSub}' ];
%         ampPredictors = [subBinStdNull{iSub,iRoi,nbins}];

%         predictors = [subBinAmpNull{iSub,iRoi,nbins} subNullPulseStd{iSub}' ];
        phPredictors = [subBinPhNull{iSub,iRoi,nbins} subBaseNullPupil{iSub}];
        for ibin=1:nbins
            [subBinAmpNullCoef(iSub,iRoi,ibin,:),BINT,subBinAmpNullResid{iSub,iRoi,ibin}] = regress(subBinStdNull{iSub,iRoi,ibin},ampPredictors);
            [subBinPhNullCoef(iSub,iRoi,ibin,:),BINT,subBinPhNullResid{iSub,iRoi,ibin}] = regress(subBinPhNull{iSub,iRoi,ibin},phPredictors);
            
%             subBinAmpNullCoef(iSub,iRoi,ibin,:) = predictors\subBinAmpNull{iSub,iRoi,ibin};
%             subBinAmpNullResid{iSub,iRoi,ibin} = subBinAmpNull{iSub,iRoi,ibin} - predictors*subBinAmpNullCoef(iSub,iRoi,ibin,:);
%             pulseKernelRoi(iSub,iRoi,:) = pulseDesignMat{iSub}'\subTseries{iSub,iRoi}';
            subBinMeanResp(iSub,iRoi,ibin,:) = mean(subBinNullTrialTseries{iSub,iRoi,ibin},2);
            estAmp = repmat(ampPredictors*squeeze(subBinAmpNullCoef(iSub,iRoi,ibin,:)),1,10);
            meanResp = repmat(zscore(squeeze(subBinMeanResp(iSub,iRoi,ibin,:))'),size(estAmp,1),1);
%             subBinNullEst{iSub,iRoi,ibin} = estAmp.*meanResp;
            subBinNullAmpEst{iSub,iRoi,ibin} = ampPredictors*squeeze(subBinAmpNullCoef(iSub,iRoi,ibin,:)).*zscore(squeeze(subBinMeanResp(iSub,iRoi,ibin,:))');
            phShifts = -floor(phPredictors*squeeze(subBinPhNullCoef(iSub,iRoi,ibin,:))*trialLength/(2*pi));%Fourier phase inverse of circshift
            phShifts(isnan(phShifts)) = zeros;
            for itrial=1:size(subBinNullAmpEst{iSub,iRoi,ibin},1)
                subBinNullEst{iSub,iRoi,ibin}(itrial,:) = circshift(subBinNullAmpEst{iSub,iRoi,ibin}(itrial,:),phShifts(itrial));
            end
 
%             subBinNullEst{iSub,iRoi,ibin} = predictors*squeeze(subBinAmpNullCoef(iSub,iRoi,ibin,:)).*squeeze(subBinMeanResp(iSub,iRoi,ibin,:))';
            subBinNullAmpResid{iSub,iRoi,ibin} = subBinNullTrialTseries{iSub,iRoi,ibin}' - subBinNullAmpEst{iSub,iRoi,ibin};
            subBinNullResid{iSub,iRoi,ibin} = subBinNullTrialTseries{iSub,iRoi,ibin} - subBinNullEst{iSub,iRoi,ibin}';
            origVar(iSub,iRoi,ibin) = var(subBinNullTrialTseries{iSub,iRoi,ibin}(:));
            residVar(iSub,iRoi,ibin) = nanvar(subBinNullResid{iSub,iRoi,ibin}(:));%some trials are NaNs!
            ampResidVar(iSub,iRoi,ibin) = nanvar(subBinNullAmpResid{iSub,iRoi,ibin}(:));%some trials are NaNs!
        end
    end
    
    iRoi=1;

    for ihemi=1:2
        hemiBinsNull{iSub,ihemi} = [subNullRespStd{iSub}',subNullRespAmp{iSub}', subNullRespMean{iSub}', subNullRespPh{iSub}'...
            subNullPulseStd{iSub}',subNullPulseAmp{iSub}', subNullPulseMean{iSub}', subNullPulsePh{iSub}'...
            subStdNullPupil{iSub}, subAmpNullPupil{iSub}, subBaseNullPupil{iSub}, subPhNullPupil{iSub}];%, ...
        for ibin=1:nbins
            hemiBinsNull{iSub,ihemi} =  [hemiBinsNull{iSub,ihemi}, subBinStdNull{iSub,1+ihemi,ibin}];
        end
        for ibin=1:nbins
            hemiBinsNull{iSub,ihemi} =  [hemiBinsNull{iSub,ihemi}, subBinAmpNull{iSub,1+ihemi,ibin}];
        end
        for ibin=1:nbins
            hemiBinsNull{iSub,ihemi} =  [hemiBinsNull{iSub,ihemi}, subBinPhNull{iSub,1+ihemi,ibin}];
        end
    end

    [corrHemiBinsNull(iSub,:,:), pvalHemiBinsNull(iSub,:,:)] = corr(hemiBinsNull{iSub,1}, hemiBinsNull{iSub,2},'rows','complete');
    
end
labels = {'resp std','resp amp', 'resp mean','resp phase',...
    'pulse std','pulse amp', 'pulse mean', 'pulse phase',...
    'pupil std', 'pupil amp', 'pupil baseline', 'pupil phase',...
    'roi amp', 'roi std', 'roi phase'};

meanBinFreqStd = std(subBinFreqMeanTrial,0,5);
meanBinContrastStd = std(subBinContrastMeanTrial,0,5);
toc
%%
ifig=0;
plotColors = {[1 0 0], [0 0 1], [0 1 0], [1 1 0]};
rows=length(roiNames);
cols=6;
ifig=ifig+1; figure(ifig); clf
linestyle = '-';
%%
ifig=ifig+1; figure(ifig); clf
% imagesc(squeeze(mean(corrPulsePupilBenson)));
imagesc(squeeze(mean(corrHemiBinsNull)));

title('trial-by-trial correlation');
yticklabels(labels);
yticks(1:length(labels));
set(gcf,'position',[100 100 800 700])
% caxis([-1 1]);

%%
respKernelRoi = zscore(respKernelRoi,0,3);
pulseKernelRoi = zscore(pulseKernelRoi,0,3);
% respPulseKernelRoi = zscore(respPulseKernelRoi,0,3);

%%
iRoi=1;
nfreqs=size(meanBinFreqStd,4);
ncontrasts =size(meanBinContrastStd,4);
ifig=ifig+1; figure(ifig); clf
plotcmap = cool;
rows=2; cols=2;
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
%contrast
subplot(rows,cols,3)
for icontrast=1:ncontrasts
    plot(squeeze(mean(meanBinContrastStd(:,iRoi,:,icontrast))),'Color',plotcmap(ceil(size(plotcmap,1)*icontrast/ncontrasts),:))
   hold on
end
xlabel('eccentricity bin');
subplot(rows,cols,4)
for ibin=1:nbins
    plot(squeeze(mean(meanBinContrastStd(:,iRoi,ibin,:))),'Color',plotcmap(ceil(size(plotcmap,1)*ibin/nbins),:))
   hold on
end
xlabel('contrast')
%%
ifig=ifig+1; figure(ifig); clf
iSub=1;
rows=5;
cols=nbins;
for ibin=1:nbins
   subplot(rows,cols,ibin)
   plot(subBinNullTrialTseries{iSub,iRoi,ibin});
   subplot(rows,cols,ibin+cols)
   plot(subBinNullAmpEst{iSub,iRoi,ibin}');
    subplot(rows,cols,ibin+2*cols)
   plot(subBinNullAmpResid{iSub,iRoi,ibin}');
   subplot(rows,cols,ibin+3*cols)
   plot(subBinNullEst{iSub,iRoi,ibin}');
   subplot(rows,cols,ibin+4*cols)
   plot(subBinNullResid{iSub,iRoi,ibin});
%    plot(subBinNullTrialTseries{iSub,iRoi,ibin},'color',plotColors{1});
%    subplot(rows,cols,ibin+cols)
%    plot(subBinNullEst{iSub,iRoi,ibin}','color',plotColors{2});
%     subplot(rows,cols,ibin+2*cols)
%    plot(subBinNullResid{iSub,iRoi,ibin}','color',plotColors{3});
end
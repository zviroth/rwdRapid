close all
mrQuit
clear all
saveFigs=1;
arousalType = 5; %1=pupil std, 2=pupil baseline, 3=pulse std, 4=pulse mean, 5=ppg std, 6=ppg mean
tic
toZscore=0;%0 or 1
concatProj= 0;
% dataFolder = '/Volumes/MH02086153MACDT-Drobo/rwdFmri/';
dataFolder = '/Users/rothzn/rwdFmri/';
% dataFolder = '/Volumes/MH02086153MACDT-Drobo/rwdFmri/';
%Load eyetracking data
saveFolder = '~/rwdFmri/';
load([saveFolder 'rwdRapidEyeData.mat'], 'subFolders', 'samplerate',  ...
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
load([saveFolder 'roiTC_' zScoreString concatProjStr '.mat'], 'subFolders', 'roiNames', ...
    'numRuns','numTRs','concatInfo',...
    'frames', 'junkedFrames', 'TR', 'trialsPerRun', 'trialLength', 'nVolumes',...
    'nbins','binBorders','binCenters',...
    'eccen','ang','areas',...
    'numVox','roiTC','nullTrialsTRs',...
    'nullTrialsRun','nullTrials','contrastTrialsRun','freqTrialsRun','contrastTrials','freqTrials',...
    'nullTseries','nullTrialTseries','stimTseries','stimTrialTseries');

goodSubs = [1:12 14];
goodSubs = [1:14];
bensonROIs = 1:3;
%%
nullPupil = origNullPupil;
for iSub=1:length(goodSubs)
    
    for rwd=1:2
        %remove the first trial of each run
        junkedIndices = zeros(size(eyeNullTrials{goodSubs(iSub),rwd}));
        junkedIndices(1:17:end) = ones;
        nullTrialIndices = find(eyeNullTrials{goodSubs(iSub),rwd}==1);
        junkedNullTrials = find(mod(nullTrialIndices,17)==1);
        nullPupil{goodSubs(iSub),rwd}(junkedNullTrials,:) = [];
        
        %get measure of pupil size per trial
        stdNullPupil{iSub,rwd} = nanstd(nullPupil{goodSubs(iSub),rwd},0,2);%per trial
        baseNullPupil{iSub,rwd} = nanmean(nullPupil{goodSubs(iSub),rwd}(:,1:baselineT),2);%per trial
        f=fft(fillmissing(nullPupil{goodSubs(iSub),rwd},'linear',2),[],2);
        phNullPupil{iSub,rwd} = angle(f(:,2));
        ampNullPupil{iSub,rwd} = abs(f(:,2));
        
        
        
        
        for iRoi=1:length(roiNames)
            %get measure of fMRI response amplitude
            roiNullTrialTseries{iSub,iRoi,rwd} = squeeze(nanmean(nullTrialTseries{goodSubs(iSub),iRoi,rwd}));%averaged over voxels
            rwdRoiStdNull{iSub,iRoi,rwd} = std(roiNullTrialTseries{iSub,iRoi,rwd});
            f = fft(roiNullTrialTseries{iSub,iRoi,rwd});
            rwdRoiAmpNull{iSub,iRoi,rwd} = abs(f(2,:));
            rwdRoiPhNull{iSub,iRoi,rwd} = angle(f(2,:));
            %correlate pupil amplitude/baseline with fMRI amplitude
            %             [corrStdPupilRoiNullRwd(iSub,iRoi,rwd) pvalsStdRwd(iSub,iRoi,rwd)] = corr(rwdRoiStdNull{iSub,iRoi,rwd}',stdNullPupil{iSub,rwd}, 'rows','complete');
            %             [corrBasePupilRoiNullRwd(iSub,iRoi,rwd) pvalsBaseRwd(iSub,iRoi,rwd)] = corr(rwdRoiStdNull{iSub,iRoi,rwd}',baseNullPupil{iSub,rwd}, 'rows','complete');
            %             [corrStdPupilAmpRoiNullRwd(iSub,iRoi,rwd) pvalsStdAmpRwd(iSub,iRoi,rwd)] = corr(rwdRoiAmpNull{iSub,iRoi,rwd}',stdNullPupil{iSub,rwd}, 'rows','complete');
            %             [corrBasePupilAmpRoiNullRwd(iSub,iRoi,rwd) pvalsBaseAmpRwd(iSub,iRoi,rwd)] = corr(rwdRoiAmpNull{iSub,iRoi,rwd}',baseNullPupil{iSub,rwd}, 'rows','complete');
        end
    end
    
    %concatenate across rwd, separately for each subject
    
    subStdNullPupil{iSub} = [stdNullPupil{iSub,1}; stdNullPupil{iSub,2}];
    subBaseNullPupil{iSub} = [baseNullPupil{iSub,1}; baseNullPupil{iSub,2}];
    subPhNullPupil{iSub} = [phNullPupil{iSub,1}; phNullPupil{iSub,2}];
    subAmpNullPupil{iSub} = [ampNullPupil{iSub,1}; ampNullPupil{iSub,2}];
    
    
    subPulse{iSub} = [rwdPulseTrials{goodSubs(iSub),1}(:,nullTrials{goodSubs(iSub),1}==1) rwdPulseTrials{goodSubs(iSub),2}(:,nullTrials{goodSubs(iSub),2}==1)];%concatenate pulse of null trials
    subPulseStd{iSub} = std(subPulse{iSub});
    subPulseMean{iSub} = mean(subPulse{iSub});
    %             f=fft(fillmissing(subPulse{iSub},'linear',2),[],2);
    f=fft(subPulse{iSub});
    subPulsePh{iSub} = angle(f(2,:));
    subPulseAmp{iSub} = abs(f(2,:));
    
    subPpg{iSub} = [rwdPpgTrials{goodSubs(iSub),1}(:,nullTrials{goodSubs(iSub),1}==1) rwdPpgTrials{goodSubs(iSub),2}(:,nullTrials{goodSubs(iSub),2}==1)];%concatenate ppg-amp of null trials
    subPpgStd{iSub} = std(subPpg{iSub});
    subPpgMean{iSub} = mean(subPpg{iSub});
    f=fft(subPpg{iSub});
    subPpgPh{iSub} = angle(f(2,:));
    subPpgAmp{iSub} = abs(f(2,:));
    
    %divide trials according to pupil response
    switch arousalType
        case 1
            subMedianStdPupil(iSub) = median(subStdNullPupil{iSub},'omitnan');
            pupilArousal{iSub,1} = subStdNullPupil{iSub}>subMedianStdPupil(iSub);
            pupilArousal{iSub,2} = subStdNullPupil{iSub}<subMedianStdPupil(iSub);
        case 2
            subMedianBasePupil(iSub) = median(subBaseNullPupil{iSub},'omitnan');
            pupilArousal{iSub,1} = subBaseNullPupil{iSub}>subMedianBasePupil(iSub);
            pupilArousal{iSub,2} = subBaseNullPupil{iSub}<subMedianBasePupil(iSub);
        case 3
            subMedianStdPulse(iSub) = median(subPulseStd{iSub},'omitnan');
            pupilArousal{iSub,1} = subPulseStd{iSub}>subMedianStdPulse(iSub);
            pupilArousal{iSub,2} = subPulseStd{iSub}<subMedianStdPulse(iSub);
        case 4
            subMedianMeanPulse(iSub) = median(subPulseMean{iSub},'omitnan');
            pupilArousal{iSub,1} = subPulseMean{iSub}>subMedianMeanPulse(iSub);
            pupilArousal{iSub,2} = subPulseMean{iSub}<subMedianMeanPulse(iSub);
        case 5
            subMedianStdPpg(iSub) = median(subPpgStd{iSub},'omitnan');
            pupilArousal{iSub,1} = subPpgStd{iSub}>subMedianStdPpg(iSub);
            pupilArousal{iSub,2} = subPpgStd{iSub}<subMedianStdPpg(iSub);
        case 6
            subMedianMeanPpg(iSub) = median(subPpgMean{iSub},'omitnan');
            pupilArousal{iSub,1} = subPpgMean{iSub}>subMedianMeanPpg(iSub);
            pupilArousal{iSub,2} = subPpgMean{iSub}<subMedianMeanPpg(iSub);
    end
    
    
    for arousal=1:2
        subArousalPulse(iSub,arousal,:) = mean(subPulse{iSub}(:,pupilArousal{iSub,arousal}),2);
        subArousalPulseStd(iSub,arousal) = nanstd(squeeze(subArousalPulse(iSub,arousal,:)));
        subArousalPulseMean(iSub,arousal) = nanmean(squeeze(subArousalPulse(iSub,arousal,:)));
    end
    
    for arousal=1:2
        subArousalPpg(iSub,arousal,:) = mean(subPpg{iSub}(:,pupilArousal{iSub,arousal}),2);
        subArousalPpgStd(iSub,arousal) = nanstd(squeeze(subArousalPpg(iSub,arousal,:)));
        subArousalPpgMean(iSub,arousal) = nanmean(squeeze(subArousalPpg(iSub,arousal,:)));
    end
    
    subResp{iSub} = [rwdRvTrials{goodSubs(iSub),1}(:,nullTrials{goodSubs(iSub),1}==1) rwdRvTrials{goodSubs(iSub),2}(:,nullTrials{goodSubs(iSub),2}==1)];%concatenate pulse of null trials
    subRespStd{iSub} = std(subResp{iSub});
    subRespMean{iSub} = mean(subResp{iSub});
    f=fft(subResp{iSub});
    subRespPh{iSub} = angle(f(2,:));
    subRespAmp{iSub} = abs(f(2,:));
    
    
    for arousal=1:2
        subArousalResp(iSub,arousal,:) = mean(subResp{iSub}(:,pupilArousal{iSub,arousal}),2);
        subArousalRespStd(iSub,arousal) = nanstd(squeeze(subArousalResp(iSub,arousal,:)));
        subArousalRespMean(iSub,arousal) = nanmean(squeeze(subArousalResp(iSub,arousal,:)));
    end
    
    %     %correlate pulse STD with pupil STD
    %     [corrStdPupilStdPulseNull(iSub) pvalsStdPupilStdPulse(iSub)] = corr(subPulseStd{iSub}',subStdNullPupil{iSub},'rows','complete');
    %     %correlate mean pulse with baseline pupil
    %     [corrBasePupilMeanPulseNull(iSub) pvalsBasePupilMeanPulse(iSub)] = corr(subPulseMean{iSub}',subBaseNullPupil{iSub},'rows','complete');
    %     %correlate std pulse with baseline pupil
    %     [corrBasePupilStdPulseNull(iSub) pvalsBasePupilStdPulse(iSub)] = corr(subPulseStd{iSub}',subBaseNullPupil{iSub},'rows','complete');
    %     %correlate mean pulse with std pupil
    %     [corrStdPupilMeanPulseNull(iSub) pvalsStdPupilMeanPulse(iSub)] = corr(subPulseMean{iSub}',subStdNullPupil{iSub},'rows','complete');
    
    %concatenate physio design matrcies
    pulseDesignMat{iSub} = [designMatPulse{iSub,1} designMatPulse{iSub,2}];
    ppgDesignMat{iSub} = [designMatPpg{iSub,1} designMatPpg{iSub,2}];
    respDesignMat{iSub} = [designMatResp{iSub,1} designMatResp{iSub,2}];
    respPulseDesignMat{iSub} = [designMatRespPulse{iSub,1} designMatRespPulse{iSub,2}] ;
    respPulsePpgDesignMat{iSub} = [designMatRespPulsePpg{iSub,1} designMatRespPulsePpg{iSub,2}] ;
    
    for iRoi=1:length(roiNames)
        
        subNullTrialTseries{iSub,iRoi} = [roiNullTrialTseries{iSub,iRoi,1} roiNullTrialTseries{iSub,iRoi,2}];
        subRoiStdNull{iSub,iRoi} = [rwdRoiStdNull{iSub,iRoi,1}'; rwdRoiStdNull{iSub,iRoi,2}'];
        subRoiAmpNull{iSub,iRoi} = [rwdRoiAmpNull{iSub,iRoi,1}'; rwdRoiAmpNull{iSub,iRoi,2}'];
        subRoiPhNull{iSub,iRoi} = [rwdRoiPhNull{iSub,iRoi,1}'; rwdRoiPhNull{iSub,iRoi,2}'];
        %         %correlate std measure with pupil
        %         [corrStdPupilStdRoiNull(iSub,iRoi) pvalsStd(iSub,iRoi)]= corr(subRoiStdNull{iSub,iRoi},subStdNullPupil{iSub}, 'rows','complete');
        %         [corrBasePupilStdRoiNull(iSub,iRoi) pvalsBase(iSub,iRoi)] = corr(subRoiStdNull{iSub,iRoi},subBaseNullPupil{iSub}, 'rows','complete');
        %         %correlate fft amp measure with pupil
        %         [corrStdPupilAmpRoiNull(iSub,iRoi) pvalsAmpStd(iSub,iRoi)]= corr(subRoiAmpNull{iSub,iRoi},subStdNullPupil{iSub}, 'rows','complete');
        %         [corrBasePupilAmpRoiNull(iSub,iRoi) pvalsAmpBase(iSub,iRoi)] = corr(subRoiAmpNull{iSub,iRoi},subBaseNullPupil{iSub}, 'rows','complete');
        
        %         %correlate std measure with pulse
        %         [corrStdPulseStdRoiNull(iSub,iRoi) pvalsStdPulse(iSub,iRoi)]= corr(subRoiStdNull{iSub,iRoi},subPulseStd{iSub}', 'rows','complete');
        %         [corrMeanPulseStdRoiNull(iSub,iRoi) pvalsMeanPulse(iSub,iRoi)]= corr(subRoiStdNull{iSub,iRoi},subPulseMean{iSub}', 'rows','complete');
        %
        %         %correlate fft amp measure with pulse
        %         [corrStdPulseAmpRoiNull(iSub,iRoi) pvalsAmpStdPulse(iSub,iRoi)]= corr(subRoiAmpNull{iSub,iRoi},subPulseStd{iSub}', 'rows','complete');
        
        %concatenate across rwd for physio regression
        subTseries{iSub,iRoi} = nanmean([roiTC{iSub,iRoi,1}.tSeries roiTC{iSub,iRoi,2}.tSeries]);%mean over voxels

        pulseKernelRoi(iSub,iRoi,:) = pulseDesignMat{iSub}'\subTseries{iSub,iRoi}';
        ppgKernelRoi(iSub,iRoi,:) = ppgDesignMat{iSub}'\subTseries{iSub,iRoi}';
        respKernelRoi(iSub,iRoi,:) = respDesignMat{iSub}'\subTseries{iSub,iRoi}';
        respPulseKernelRoi(iSub,iRoi,:) = respPulseDesignMat{iSub}'\subTseries{iSub,iRoi}';
        respPulsePpgKernelRoi(iSub,iRoi,:) = respPulsePpgDesignMat{iSub}'\subTseries{iSub,iRoi}';
        
        %divide into 2 arousal bins
        for arousal=1:2
            subArousalRoiMeanTrial(iSub,iRoi,arousal,:) = mean(subNullTrialTseries{iSub,iRoi}(:,pupilArousal{iSub,arousal}),2);
            subArousalRoiMeanStd(iSub,iRoi,arousal) = mean(subRoiStdNull{iSub,iRoi}(pupilArousal{iSub,arousal}));
            subArousalRoiMeanAmp(iSub,iRoi,arousal) = mean(subRoiAmpNull{iSub,iRoi}(pupilArousal{iSub,arousal}));
        end
    end
    
    iRoi=1;
    pulsePupilBenson{iSub} = [subRespStd{iSub}', subRespMean{iSub}',...
        subPulseStd{iSub}', subPulseMean{iSub}', ...
        subPpgStd{iSub}', subPpgMean{iSub}', ...
        subStdNullPupil{iSub}, subAmpNullPupil{iSub}, subBaseNullPupil{iSub}, subPhNullPupil{iSub}, ...
        subRoiAmpNull{iSub,iRoi}, subRoiStdNull{iSub,iRoi}, subRoiPhNull{iSub,iRoi}];
    [corrPulsePupilBenson(iSub,:,:) pvalPulsePupilBenson(iSub,:,:)] = corr(pulsePupilBenson{iSub},'rows','complete');
    
    
    
end

%%
eccMin = 0.2;
eccMax = 70;
nbins = 8;
binBorders = logspace(log10(eccMin),log10(eccMax),nbins+1);
nbins = length(binBorders)-1;
for i=2:length(binBorders)
    binCenters(i-1) = (binBorders(i)+binBorders(i-1))/2;
end
for iSub=1:length(goodSubs)
    for iRoi=bensonROIs%1:length(roiNames)
        for ibin=1:nbins
            for rwd=1:2
                binVoxels = eccen{goodSubs(iSub),iRoi}>binBorders(ibin) & eccen{goodSubs(iSub),iRoi}<=binBorders(ibin+1);
                numBinVoxels(iSub,iRoi,ibin) = sum(binVoxels);
                %                 binVoxels = binVoxels & areas{iSub,iRoi}==1;%ONLY V1
                binVoxels = binVoxels & areas{goodSubs(iSub),iRoi}>0;%ONLY V1,V2,V3
                binMeanTseries = nanmean(roiTC{goodSubs(iSub),iRoi,rwd}.tSeries(binVoxels,:),1);%mean timecourse across voxels
                binNullTseries{iSub,iRoi,ibin,rwd} = binMeanTseries(nullTrialsTRs{goodSubs(iSub),rwd}==1);
                binNullTrialTseries{iSub,iRoi,ibin,rwd} = reshape(binNullTseries{iSub,iRoi,ibin,rwd},trialLength,[]);%iSub,T,trial
                %get fMRI response amplitude
                rwdBinStdNull{iSub,iRoi,ibin,rwd} = std(binNullTrialTseries{iSub,iRoi,ibin,rwd});%
                f = fft(binNullTrialTseries{iSub,iRoi,ibin,rwd});
                rwdBinAmpNull{iSub,iRoi,ibin,rwd} = abs(f(2,:));
                rwdBinPhNull{iSub,iRoi,ibin,rwd} = angle(f(2,:));
                
                %mean across trials
                subRwdBinMeanTrial(iSub,iRoi,ibin,rwd,:) = mean(binNullTrialTseries{iSub,iRoi,ibin,rwd},2);
                subRwdBinMeanStd(iSub,iRoi,ibin,rwd) = mean(rwdBinStdNull{iSub,iRoi,ibin,rwd},2);
                subRwdBinMeanAmp(iSub,iRoi,ibin,rwd) = mean(rwdBinAmpNull{iSub,iRoi,ibin,rwd},2);
                subRwdBinVarStd(iSub,iRoi,ibin,rwd) = std(rwdBinStdNull{iSub,iRoi,ibin,rwd},0,2);
                subRwdBinVarAmp(iSub,iRoi,ibin,rwd) = std(rwdBinAmpNull{iSub,iRoi,ibin,rwd},0,2);
                subRwdBinVarPh(iSub,iRoi,ibin,rwd) = circ_std(rwdBinPhNull{iSub,iRoi,ibin,rwd},[],[],2);
                subRwdBinMeanVar(iSub,iRoi,ibin,rwd) = mean(std(binNullTrialTseries{iSub,iRoi,ibin,rwd},0,1));
            end
            %concatenate across rwd
            subBinNullTrialTseries{iSub,iRoi,ibin} = [binNullTrialTseries{iSub,iRoi,ibin,1} binNullTrialTseries{iSub,iRoi,ibin,2}];
            subBinStdNull{iSub,iRoi,ibin} = [rwdBinStdNull{iSub,iRoi,ibin,1}'; rwdBinStdNull{iSub,iRoi,ibin,2}'];
            subBinAmpNull{iSub,iRoi,ibin} = [rwdBinAmpNull{iSub,iRoi,ibin,1}'; rwdBinAmpNull{iSub,iRoi,ibin,2}'];
            subBinPhNull{iSub,iRoi,ibin} = [rwdBinPhNull{iSub,iRoi,ibin,1}'; rwdBinPhNull{iSub,iRoi,ibin,2}'];
            
            

            %             %correlate with pupil
            %             [corrStdPupilStdBinNull(iSub,iRoi,ibin) pvalsBinStd(iSub,iRoi,ibin)]= corr(subBinStdNull{iSub,iRoi,ibin},subStdNullPupil{iSub}, 'rows','complete');
            %             [corrStdPupilAmpBinNull(iSub,iRoi,ibin) pvalsBinAmp(iSub,iRoi,ibin)]= corr(subBinAmpNull{iSub,iRoi,ibin},subStdNullPupil{iSub}, 'rows','complete');
            %
            %             %correlate with pulse
            %             [corrStdPulseStdBinNull(iSub,iRoi,ibin) pvalsBinStdPulse(iSub,iRoi,ibin)]= corr(subBinStdNull{iSub,iRoi,ibin},subPulseStd{iSub}', 'rows','complete');
            %             [corrStdPulseAmpBinNull(iSub,iRoi,ibin) pvalsBinAmpPulse(iSub,iRoi,ibin)]= corr(subBinAmpNull{iSub,iRoi,ibin},subPulseStd{iSub}', 'rows','complete');
            
            
            %divide into 2 arousal bins
            for arousal=1:2
                subArousalBinMeanTrial(iSub,iRoi,ibin,arousal,:) = mean(subBinNullTrialTseries{iSub,iRoi,ibin}(:,pupilArousal{iSub,arousal}),2);
                subArousalBinMeanStd(iSub,iRoi,ibin,arousal) = mean(subBinStdNull{iSub,iRoi,ibin}(pupilArousal{iSub,arousal}));
                subArousalBinMeanAmp(iSub,iRoi,ibin,arousal) = mean(subBinAmpNull{iSub,iRoi,ibin}(pupilArousal{iSub,arousal}));
                subArousalBinVarStd(iSub,iRoi,ibin,arousal) = std(subBinStdNull{iSub,iRoi,ibin}(pupilArousal{iSub,arousal}));
                %                 subArousalBinVarAmp(iSub,iRoi,ibin,arousal) = std(subBinAmpNull{iSub,iRoi,ibin}(pupilArousal{iSub,arousal}));
                f = fft(subBinNullTrialTseries{iSub,iRoi,ibin}(:,pupilArousal{iSub,arousal}));
                [s subArousalBinVarPh(iSub,iRoi,ibin,arousal)] = circ_std(angle(f(2,:)'));
                subArousalBinVarAmp(iSub,iRoi,ibin,arousal) = std(abs(f(2,:)));
                subArousalBinVarTime(iSub,iRoi,ibin,arousal) = mean(std(subBinNullTrialTseries{iSub,iRoi,ibin}(:,pupilArousal{iSub,arousal}),0,2));
            end
            
        end
    end
    
    iRoi=1;
    %     pulsePupilBins{iSub} = [subRespStd{iSub}', subRespMean{iSub}',subPulseStd{iSub}', subPulseMean{iSub}', ...
    %         subStdNullPupil{iSub}, subBaseNullPupil{iSub}];%, ...
    %subRoiAmpNull{iSub,iRoi}, subRoiStdNull{iSub,iRoi}, subRoiPhNull{iSub,iRoi}];
    
    for ihemi=1:2
        hemiBins{iSub,ihemi} = [subRespStd{iSub}',subRespAmp{iSub}', subRespMean{iSub}', subRespPh{iSub}'...
            subPulseStd{iSub}',subPulseAmp{iSub}', subPulseMean{iSub}', subPulsePh{iSub}'...
            subPpgStd{iSub}', subPpgAmp{iSub}', subPpgMean{iSub}', subPpgPh{iSub}'...
            subStdNullPupil{iSub}, subAmpNullPupil{iSub}, subBaseNullPupil{iSub}, subPhNullPupil{iSub}];%, ...
        for ibin=1:nbins
            hemiBins{iSub,ihemi} =  [hemiBins{iSub,ihemi}, subBinStdNull{iSub,1+ihemi,ibin}];
        end
        for ibin=1:nbins
            hemiBins{iSub,ihemi} =  [hemiBins{iSub,ihemi}, subBinAmpNull{iSub,1+ihemi,ibin}];
        end
        for ibin=1:nbins
            hemiBins{iSub,ihemi} =  [hemiBins{iSub,ihemi}, subBinPhNull{iSub,1+ihemi,ibin}];
        end
        %DMN & cerebellum
        for iRoi=1:2
            curRoi = bensonROIs(end)+2*(iRoi-1)+ihemi;
            hemiBins{iSub,ihemi} = [hemiBins{iSub,ihemi} subRoiAmpNull{iSub,curRoi} subRoiStdNull{iSub,curRoi} subRoiPhNull{iSub,curRoi}];
        end
    end
    %     for ibin=1:nbins
    %         pulsePupilBins{iSub} = [pulsePupilBins{iSub}, subBinStdNull{iSub,iRoi,ibin}];
    %     end
    %     for ibin=1:nbins
    %         pulsePupilBins{iSub} = [pulsePupilBins{iSub}, subBinAmpNull{iSub,iRoi,ibin}];
    %     end
    [corrHemiBins(iSub,:,:), pvalHemiBins(iSub,:,:)] = corr(hemiBins{iSub,1}, hemiBins{iSub,2},'rows','complete');
    
end
labels = {'resp std','resp amp', 'resp mean','resp phase',...
    'pulse std','pulse amp', 'pulse mean', 'pulse phase',...
    'ppg std','ppg amp', 'ppg mean', 'ppg phase',...
    'pupil std', 'pupil amp', 'pupil baseline', 'pupil phase',...
    'benson amp', 'benson std', 'benson phase','DMN amp','DMN std','DMN phase','cerebellum amp','cerebellum std','cerebellum phase'};

toc
% %%
% for iSub=1:length(goodSubs)
%     for ihemi=1:2
%         for ibin=1:nbins
%             subBinAmpNull{iSub,1+ihemi,ibin};
%             subBinPhNull{iSub,1+ihemi,ibin};
%             complexTraj{iSub}(ihemi,ibin,:) = subBinAmpNull{iSub,1+ihemi,ibin}.*exp(subBinPhNull{iSub,1+ihemi,ibin}*1i);
%         end
%         %DMN & cerebellum
%         for iRoi=1:2
%             curRoi = bensonROIs(end)+2*(iRoi-1)+ihemi;
%             complexTraj{iSub}(ihemi,nbins+iRoi,:) = subRoiAmpNull{iSub,curRoi}.*exp(subRoiPhNull{iSub,curRoi}*1i);
%         end
%         %get cartesion coordinates
%         for iRoi=1:size(complexTraj{iSub},2)
%             traj{iSub}(ihemi,iRoi,:,:) = [real(squeeze(complexTraj{iSub}(ihemi,iRoi,:))) imag(squeeze(complexTraj{iSub}(ihemi,iRoi,:)))];
%         end
%     end
% 
%     for iRoi1=1:size(complexTraj{iSub},2)
%         for iRoi2=1:size(complexTraj{iSub},2)
%             [procDist(iSub,iRoi1,iRoi2) procTrans{iSub}(iRoi1,iRoi2,:,:) procTransform{iSub,iRoi1,iRoi2}] = procrustes(squeeze(traj{iSub}(1,iRoi1,:,:)),squeeze(traj{iSub}(2,iRoi2,:,:)));%transforms the 2nd to fit the 1st
%         end
%     end
% end
%%
ifig=0;
plotColors = {[1 0 0], [0 0 1], [0 1 0], [1 1 0]};

%%
rows=length(bensonROIs);
cols=6;
ifig=ifig+1; figure(ifig); clf
linestyle = '.-';
for iRoi=bensonROIs%1:length(roiNames)
    subplot(rows,cols,1+ (iRoi-1)*cols)
    for rwd=1:2
        plot(squeeze(mean(subRwdBinMeanStd(:,iRoi,:,rwd))),linestyle,'color',plotColors{rwd});
        hold on
    end
    title({roiNames{iRoi}; 'std'});
    if iRoi==1
        legend('high','low');
    end
    
    subplot(rows,cols,2+ (iRoi-1)*cols)
    for rwd=1:2
        plot(squeeze(mean(subRwdBinMeanAmp(:,iRoi,:,rwd))),linestyle,'color',plotColors{rwd});
        hold on
    end
    title('fft amp');
    
    subplot(rows,cols,3+ (iRoi-1)*cols)
    for rwd=1:2
        plot(squeeze(mean(subRwdBinVarStd(:,iRoi,:,rwd))),linestyle,'color',plotColors{rwd});
        hold on
    end
    title('std variability');
    
    subplot(rows,cols,4+ (iRoi-1)*cols)
    for rwd=1:2
        plot(squeeze(mean(subRwdBinVarAmp(:,iRoi,:,rwd))),linestyle,'color',plotColors{rwd});
        hold on
    end
    title('amp variability');
    
    subplot(rows,cols,5+ (iRoi-1)*cols)
    for rwd=1:2
        plot(squeeze(mean(subRwdBinVarPh(:,iRoi,:,rwd))),linestyle,'color',plotColors{rwd});
        hold on
    end
    title('ph variability');
    
    subplot(rows,cols,6+ (iRoi-1)*cols)
    for rwd=1:2
        plot(squeeze(mean(subRwdBinMeanVar(:,iRoi,:,rwd))),linestyle,'color',plotColors{rwd});
        hold on
    end
    title('timepoint variability');
    

end
for isubplot=1:rows*cols
    subplot(rows,cols,isubplot)
   xlim([0 nbins+1]); 
end
set(gcf,'position',[100 100 700 400])
if saveFigs
print('-painters','-dpdf',['~/Documents/MATLAB/rwdRapid/figures/rwdBinsNull.pdf']);
end
%%
rows=length(bensonROIs);%roiNames);
cols=6;
ifig=ifig+1; figure(ifig); clf
linestyle = '.-';
for iRoi=bensonROIs%1:length(roiNames)
    subplot(rows,cols,1+ (iRoi-1)*cols)
    for arousal=1:2
        plot(squeeze(mean(subArousalBinMeanStd(:,iRoi,:,arousal))),linestyle,'color',plotColors{arousal});
        hold on
    end
    title({roiNames{iRoi}; 'std'});
    if iRoi==1
        legend('high','low');
    end
    
    subplot(rows,cols,2+ (iRoi-1)*cols)
    for arousal=1:2
        plot(squeeze(mean(subArousalBinMeanAmp(:,iRoi,:,arousal))),linestyle,'color',plotColors{arousal});
        hold on
    end
    title('fft amp');
    
    subplot(rows,cols,3+ (iRoi-1)*cols)
    for arousal=1:2
        plot(squeeze(mean(subArousalBinVarStd(:,iRoi,:,arousal))),linestyle,'color',plotColors{arousal});
        hold on
    end
    title('std variability');
    
    subplot(rows,cols,4+ (iRoi-1)*cols)
    for arousal=1:2
        plot(squeeze(mean(subArousalBinVarAmp(:,iRoi,:,arousal))),linestyle,'color',plotColors{arousal});
        hold on
    end
    title('amp variability');
    
    subplot(rows,cols,5+ (iRoi-1)*cols)
    for arousal=1:2
        plot(squeeze(mean(subArousalBinVarPh(:,iRoi,:,arousal))),linestyle,'color',plotColors{arousal});
        hold on
    end
    title('phase variability');
    
    subplot(rows,cols,6+ (iRoi-1)*cols)
    for arousal=1:2
        plot(squeeze(mean(subArousalBinVarTime(:,iRoi,:,arousal))),linestyle,'color',plotColors{arousal});
        hold on
    end
    title('timepoint variability');
end
for isubplot=1:rows*cols
    subplot(rows,cols,isubplot)
   xlim([0 nbins+1]); 
end
set(gcf,'position',[100 100 700 400])
if saveFigs
    switch arousalType
        case 1
            arousalTypeStr='stdPupil';
        case 2
            arousalTypeStr='baselinePupil';
        case 3
            arousalTypeStr='stdPulse';
        case 4
            arousalTypeStr='meanPulse';
        case 5
            arousalTypeStr='stdPpg';
        case 6
            arousalTypeStr='meanPpg';
    end
    print('-painters','-dpdf',['~/Documents/MATLAB/rwdRapid/figures/arousalBinsNull_' arousalTypeStr '.pdf']);
end
%%
ifig=ifig+1; figure(ifig); clf
% imagesc(squeeze(mean(corrPulsePupilBenson)));
imagesc(squeeze(nanmean(corrHemiBins)));

title('trial-by-trial correlation');
yticklabels(labels);
yticks([1:length(labels)-9 length(labels)-9+round(nbins/2) length(labels)-9+round(nbins/2)+nbins ...
    length(labels)-9+round(nbins/2)+2*nbins 1+length(labels)-9+3*nbins:length(labels)+3*(nbins-1) ] );
xticklabels('');

set(gcf,'position',[100 100 600 500])
% caxis([-1 1]);
if saveFigs
print('-painters','-dpdf',['~/Documents/MATLAB/rwdRapid/figures/binNullTrialCorr' concatProjStr '.pdf']);
end
%%
respKernelRoi = zscore(respKernelRoi,0,3);
pulseKernelRoi = zscore(pulseKernelRoi,0,3);
% respPulseKernelRoi = zscore(respPulseKernelRoi,0,3);

ifig=ifig+1; figure(ifig); clf
rows=1;
cols=6;
subplot(rows,cols,1)
plot(squeeze(respKernelRoi(:,1,:))')
hold on
plot(squeeze(mean(respKernelRoi(:,1,:))),'k','linewidth',2)
title('respiration');
subplot(rows,cols,2)
plot(squeeze(pulseKernelRoi(:,1,:))')
hold on
plot(squeeze(mean(pulseKernelRoi(:,1,:))),'k','linewidth',2)
title('pulse')

subplot(rows,cols,3)
plot(squeeze(ppgKernelRoi(:,1,:))')
hold on
plot(squeeze(mean(ppgKernelRoi(:,1,:))),'k','linewidth',2)
title('ppg')

% subplot(rows,cols,4)
% plot(squeeze(respPulseKernelRoi(:,1,1:10))')
% hold on
% plot(squeeze(mean(respPulseKernelRoi(:,1,1:10))),'k','linewidth',2)
% subplot(rows,cols,5)
% plot(squeeze(respPulseKernelRoi(:,1,11:20))')
% hold on
% plot(squeeze(mean(respPulseKernelRoi(:,1,11:20))),'k','linewidth',2)

subplot(rows,cols,4)
plot(squeeze(respPulsePpgKernelRoi(:,1,1:10))')
hold on
plot(squeeze(mean(respPulsePpgKernelRoi(:,1,1:10))),'k','linewidth',2)
subplot(rows,cols,5)
plot(squeeze(respPulsePpgKernelRoi(:,1,11:20))')
hold on
plot(squeeze(mean(respPulsePpgKernelRoi(:,1,11:20))),'k','linewidth',2)
subplot(rows,cols,6)
plot(squeeze(respPulsePpgKernelRoi(:,1,21:30))')
hold on
plot(squeeze(mean(respPulsePpgKernelRoi(:,1,21:30))),'k','linewidth',2)
%%
ifig=ifig+1; figure(ifig); clf
rows=2;
cols=2;
%get mean pupil response for high and low reward
subplot(rows,cols,1)
for rwd=1:2
    allRwdNullPupil{rwd} = [];
    for iSub=1:length(goodSubs)
        allRwdNullPupil{rwd} = [allRwdNullPupil{rwd}; nullPupil{goodSubs(iSub),rwd}(:,1:samplerate*trialLength*TR)];

%         plot(nanmean(nullPupil{goodSubs(iSub),rwd}),'color',plotColors{rwd});
        hold on
        
%         pupilArousal{iSub,iarousal}(1:size(nullPupil{goodSubs(iSub),1}))
    end
    plot(nanmean(allRwdNullPupil{rwd}),'color',plotColors{rwd},'linewidth',2);
end
legend('high','low');
title('pupil high/low reward');
subplot(rows,cols,2)
plot(nanmean(allRwdNullPupil{1})-nanmean(allRwdNullPupil{2}),'k','linewidth',2);
hold on
hline(0)


%arousal
subplot(rows,cols,3)
for arousal=1:2
    allArousalNullPupil{arousal} = [];
    for iSub=1:length(goodSubs)
       subPupil = [nullPupil{goodSubs(iSub),1}(:,1:samplerate*trialLength*TR); nullPupil{goodSubs(iSub),2}(:,1:samplerate*trialLength*TR)];
       allArousalNullPupil{arousal} = [allArousalNullPupil{arousal}; subPupil(pupilArousal{iSub,arousal}>0,:)];
    end
    plot(nanmean(allArousalNullPupil{arousal}),'color',plotColors{arousal},'linewidth',2);
    hold on
end
title('pupil high/low arousal');
subplot(rows,cols,4)
plot(nanmean(allArousalNullPupil{1})-nanmean(allArousalNullPupil{2}),'k','linewidth',2);
hold on
hline(0)
if saveFigs
    print('-painters','-dpdf',['~/Documents/MATLAB/rwdRapid/figures/pupilRwdArousal_' arousalTypeStr '.pdf']);
end
%%
%%
ifig=ifig+1; figure(ifig); clf
rows=2;
cols=2;
%get mean pupil response for high and low reward
subplot(rows,cols,1)
for rwd=1:2
    allRwdNullPulse{rwd} = [];
    for iSub=1:length(goodSubs)        
        allRwdNullPulse{rwd} = [allRwdNullPulse{rwd}; rwdPulseTrials{goodSubs(iSub),rwd}(:,nullTrials{goodSubs(iSub),1}==1)'];
    end
    plot(nanmean(allRwdNullPulse{rwd}),'color',plotColors{rwd},'linewidth',2);
    hold on
end
legend('high','low');
title('pulse high/low reward');
subplot(rows,cols,2)
plot(nanmean(allRwdNullPulse{1})-nanmean(allRwdNullPulse{2}),'k','linewidth',2);
hold on
hline(0)

%arousal
subplot(rows,cols,3)
for arousal=1:2
    allArousalNullPulse{arousal} = [];
    for iSub=1:length(goodSubs)
       allArousalNullPulse{arousal} = [allArousalNullPulse{arousal}; subPulse{iSub}(:,pupilArousal{iSub,arousal}>0)'];
    end
    plot(nanmean(allArousalNullPulse{arousal}),'color',plotColors{arousal},'linewidth',2);
    hold on
end
title('pulse high/low arousal');
subplot(rows,cols,4)
plot(nanmean(allArousalNullPulse{1})-nanmean(allArousalNullPulse{2}),'k','linewidth',2);
hold on
hline(0)
if saveFigs
    print('-painters','-dpdf',['~/Documents/MATLAB/rwdRapid/figures/pulseRwdArousal_' arousalTypeStr '.pdf']);
end

%%
ifig=ifig+1; figure(ifig); clf
rows=2;
cols=2;
%get mean ppg response for high and low reward
subplot(rows,cols,1)
for rwd=1:2
    allRwdNullPpg{rwd} = [];
    for iSub=1:length(goodSubs)        
        allRwdNullPpg{rwd} = [allRwdNullPpg{rwd}; rwdPpgTrials{goodSubs(iSub),rwd}(:,nullTrials{goodSubs(iSub),1}==1)'];
    end
    plot(nanmean(allRwdNullPpg{rwd}),'color',plotColors{rwd},'linewidth',2);
    hold on
end
legend('high','low');
title('ppg high/low reward');
subplot(rows,cols,2)
plot(nanmean(allRwdNullPpg{1})-nanmean(allRwdNullPpg{2}),'k','linewidth',2);
hold on
hline(0)

%arousal
subplot(rows,cols,3)
for arousal=1:2
    allArousalNullPpg{arousal} = [];
    for iSub=1:length(goodSubs)
       allArousalNullPpg{arousal} = [allArousalNullPpg{arousal}; subPpg{iSub}(:,pupilArousal{iSub,arousal}>0)'];
    end
    plot(nanmean(allArousalNullPpg{arousal}),'color',plotColors{arousal},'linewidth',2);
    hold on
end
title('ppg high/low arousal');
subplot(rows,cols,4)
plot(nanmean(allArousalNullPpg{1})-nanmean(allArousalNullPpg{2}),'k','linewidth',2);
hold on
hline(0)
if saveFigs
    print('-painters','-dpdf',['~/Documents/MATLAB/rwdRapid/figures/ppgRwdArousal_' arousalTypeStr '.pdf']);
end

% %%
% ifig=ifig+1; figure(ifig); clf
% imagesc(squeeze(mean(procDist)));
% hold on
% xticklabels('')
% yticklabels({'benson','DMN','cerebellum'})
% yticks([1 nbins+1 nbins+2]);
% colorbar
% if saveFigs
%     print('-painters','-dpdf',['~/Documents/MATLAB/rwdRapid/figures/procrustesDistMatrix_2D.pdf']);
% end
% %%
% ifig=ifig+1; figure(ifig); clf
% rows = length(goodSubs);
% cols = size(traj{1},2)+1;
% % iSub=1; 
% iRoi1=6; 
% 
% 
% markerSize = 20;
% iRoi1=6;
% for iSub=1:length(goodSubs)
%    subplot(rows,cols,(iSub-1)*cols+1);
%    targetTraj = squeeze(traj{iSub}(1,iRoi1,:,:));%target trajectory
%    colors = cool(size(targetTraj,1));
%    scatter(targetTraj(:,1),targetTraj(:,2),markerSize,colors,'filled');
%    for  iRoi2=1:size(traj{iSub},2)
%        subplot(rows,cols,(iSub-1)*cols+iRoi2+1);
%        scatter(targetTraj(:,1),targetTraj(:,2),markerSize,colors);
%        hold on
%        fitTraj = squeeze(procTrans{iSub}(iRoi1,iRoi2,:,:));
%        % plot(temp(:,1),temp(:,2),'r.-');
%        scatter(fitTraj(:,1),fitTraj(:,2),markerSize,colors,'filled');
% 
%    end
% end
% for isubplot=1:cols*rows
%     subplot(rows,cols,isubplot)
%            xticks([]);
%        yticks([]);
% end
% %%
% ifig=ifig+1; figure(ifig); clf
% iRoi1=6;
% iRoi2 = [6 9];
% rows=length(iRoi2);
% cols=4;
% markerSize=30;
% iSub=1;
% 
% 
% for ind2=1:length(iRoi2)
%     subplot(rows,cols,1 + (ind2-1)*cols)
%     targetTraj = squeeze(traj{iSub}(1,iRoi1,:,:));%target trajectory
%     colors = cool(size(targetTraj,1));
%     scatter(targetTraj(:,1),targetTraj(:,2),markerSize,colors,'filled');
%     title('target');
%     subplot(rows,cols,2 + (ind2-1)*cols)
%     origTraj = squeeze(traj{iSub}(2,iRoi2(ind2),:,:));%target trajectory
%     scatter(origTraj(:,1),origTraj(:,2),markerSize,colors,'filled');
%     title('original');
%     subplot(rows,cols,3 + (ind2-1)*cols)
%     scatter(targetTraj(:,1),targetTraj(:,2),markerSize/2,colors);
%     hold on
%     fitTraj = squeeze(procTrans{iSub}(iRoi1,iRoi2(ind2),:,:));
%     % plot(temp(:,1),temp(:,2),'r.-');
%     scatter(fitTraj(:,1),fitTraj(:,2),markerSize,colors,'filled');
%     title('transformed');
%     %
%     % procTransform{iSub,iRoi1,iRoi2};
%     c = procTransform{iSub,iRoi1,iRoi2(ind2)}.c;
%     T = procTransform{iSub,iRoi1,iRoi2(ind2)}.T;
%     b = procTransform{iSub,iRoi1,iRoi2(ind2)}.b;
%     % Z = b*Y*T + c;
%     invb = 1/b; %=inv(b);
%     invT = inv(T);
%     estTraj = invb*(targetTraj - c)*invT;
%     subplot(rows,cols,4 + (ind2-1)*cols)
%     scatter(origTraj(:,1),origTraj(:,2),markerSize/2,colors);
%     hold on
%     scatter(estTraj(:,1),estTraj(:,2),markerSize,colors,'filled');
%     title('estimated');
% end
% for isubplot=1:cols*rows
%     subplot(rows,cols,isubplot)
%     %            xticks([]);
%     %        yticks([]);
%     axis square
% end
% if saveFigs
%     print('-painters','-dpdf',['~/Documents/MATLAB/rwdRapid/figures/procrustesExamples.pdf']);
% end
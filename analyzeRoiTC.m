close all
mrQuit
clear all
tic
toZscore=0;%0 or 1
concatProj= 0;
dataFolder = '/Volumes/MH02086153MACDT-Drobo/rwdFmri/';

%Load eyetracking data
load([dataFolder 'rwdRapidEyeData.mat'], 'dataFolder', 'subFolders', 'samplerate',  ...
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
load([dataFolder 'rwdRapid_physio' ConcatProjStr '.mat'], 'concatInfo',  ...
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

goodSubs = [1:12 14];
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
        
        
        for iRoi=1:length(roiNames)
            %get measure of fMRI response amplitude
            roiNullTrialTseries{iSub,iRoi,rwd} = squeeze(nanmean(nullTrialTseries{goodSubs(iSub),iRoi,rwd}));%averaged over voxels
            rwdRoiStdNull{iSub,iRoi,rwd} = std(roiNullTrialTseries{iSub,iRoi,rwd});
            f = fft(roiNullTrialTseries{iSub,iRoi,rwd});
            rwdRoiAmpNull{iSub,iRoi,rwd} = abs(f(2,:));
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
    
    %divide trials according to pupil response
    subMedianStdPupil(iSub) = median(subStdNullPupil{iSub},'omitnan');
    subMedianBasePupil(iSub) = median(subBaseNullPupil{iSub},'omitnan');
    pupilArousal{iSub,1} = subStdNullPupil{iSub}>subMedianStdPupil(iSub);
    pupilArousal{iSub,2} = subStdNullPupil{iSub}<subMedianStdPupil(iSub);
    
    subPulse{iSub} = [rwdPulseTrials{iSub,1} rwdPulseTrials{iSub,2}];
    
    for arousal=1:2
        subArousalPulse(iSub,arousal,:) = mean(subPulse{iSub}(:,pupilArousal{iSub,arousal}),2);
        subArousalPulseStd(iSub,arousal) = nanstd(squeeze(subArousalPulse(iSub,arousal,:)));
        subArousalPulseMean(iSub,arousal) = nanmean(squeeze(subArousalPulse(iSub,arousal,:)));
    end
    
    for iRoi=1:length(roiNames)
        %concatenate across rwd
        subNullTrialTseries{iSub,iRoi} = [roiNullTrialTseries{iSub,iRoi,1} roiNullTrialTseries{iSub,iRoi,2}];
        subRoiStdNull{iSub,iRoi} = [rwdRoiStdNull{iSub,iRoi,1}'; rwdRoiStdNull{iSub,iRoi,2}'];
        subRoiAmpNull{iSub,iRoi} = [rwdRoiAmpNull{iSub,iRoi,1}'; rwdRoiAmpNull{iSub,iRoi,2}'];
        %using std measure
        [corrStdPupilStdRoiNull(iSub,iRoi) pvalsStd(iSub,iRoi)]= corr(subRoiStdNull{iSub,iRoi},subStdNullPupil{iSub}, 'rows','complete');
        [corrBasePupilStdRoiNull(iSub,iRoi) pvalsBase(iSub,iRoi)] = corr(subRoiStdNull{iSub,iRoi},subBaseNullPupil{iSub}, 'rows','complete');
        %using fft amp measure
        [corrStdPupilAmpRoiNull(iSub,iRoi) pvalsAmpStd(iSub,iRoi)]= corr(subRoiAmpNull{iSub,iRoi},subStdNullPupil{iSub}, 'rows','complete');
        [corrBasePupilAmpRoiNull(iSub,iRoi) pvalsAmpBase(iSub,iRoi)] = corr(subRoiAmpNull{iSub,iRoi},subBaseNullPupil{iSub}, 'rows','complete');
        %divide into 2 arousal bins
        for arousal=1:2
            subArousalRoiMeanTrial(iSub,iRoi,arousal,:) = mean(subNullTrialTseries{iSub,iRoi}(:,pupilArousal{iSub,arousal}),2);
            subArousalRoiMeanStd(iSub,iRoi,arousal) = mean(subRoiStdNull{iSub,iRoi}(pupilArousal{iSub,arousal}));
            subArousalRoiMeanAmp(iSub,iRoi,arousal) = mean(subRoiAmpNull{iSub,iRoi}(pupilArousal{iSub,arousal}));
        end
    end
end

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
                
            end
            %concatenate across rwd
            subBinNullTrialTseries{iSub,iRoi,ibin} = [binNullTrialTseries{iSub,iRoi,ibin,1} binNullTrialTseries{iSub,iRoi,ibin,2}];
            subBinStdNull{iSub,iRoi,ibin} = [rwdBinStdNull{iSub,iRoi,ibin,1}'; rwdBinStdNull{iSub,iRoi,ibin,2}'];
            subBinAmpNull{iSub,iRoi,ibin} = [rwdBinAmpNull{iSub,iRoi,ibin,1}'; rwdBinAmpNull{iSub,iRoi,ibin,2}'];
            %correlate with pupil
            [corrStdPupilStdBinNull(iSub,iRoi,ibin) pvalsBinStd(iSub,iRoi,ibin)]= corr(subBinStdNull{iSub,iRoi,ibin},subStdNullPupil{iSub}, 'rows','complete');
            [corrStdPupilAmpBinNull(iSub,iRoi,ibin) pvalsBinAmp(iSub,iRoi,ibin)]= corr(subBinAmpNull{iSub,iRoi,ibin},subStdNullPupil{iSub}, 'rows','complete');
            
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
end

toc
%%
ifig=0;
plotColors = {[1 0 0], [0 0 1], [0 1 0], [1 1 0]};
rows=length(roiNames);
cols=6;
ifig=ifig+1; figure(ifig); clf
linestyle = '-';
for iRoi=1:length(roiNames)
    subplot(rows,cols,1+ (iRoi-1)*cols)
    for arousal=1:2
        plot(squeeze(mean(subArousalBinMeanStd(:,iRoi,:,arousal))),'color',plotColors{arousal},'linestyle',linestyle);
        hold on
    end
    title('std');
    
    subplot(rows,cols,2+ (iRoi-1)*cols)
    for arousal=1:2
        plot(squeeze(mean(subArousalBinMeanAmp(:,iRoi,:,arousal))),'color',plotColors{arousal},'linestyle',linestyle);
        hold on
    end
    title('fft amp');
    
    subplot(rows,cols,3+ (iRoi-1)*cols)
    for arousal=1:2
        plot(squeeze(mean(subArousalBinVarStd(:,iRoi,:,arousal))),'color',plotColors{arousal},'linestyle',linestyle);
        hold on
    end
    title('std variability');
    
    subplot(rows,cols,4+ (iRoi-1)*cols)
    for arousal=1:2
        plot(squeeze(mean(subArousalBinVarAmp(:,iRoi,:,arousal))),'color',plotColors{arousal},'linestyle',linestyle);
        hold on
    end
    title('amp variability');
    
    subplot(rows,cols,5+ (iRoi-1)*cols)
    for arousal=1:2
        plot(squeeze(mean(subArousalBinVarPh(:,iRoi,:,arousal))),'color',plotColors{arousal},'linestyle',linestyle);
        hold on
    end
    title('phase variability');
    
    subplot(rows,cols,6+ (iRoi-1)*cols)
    for arousal=1:2
        plot(squeeze(mean(subArousalBinVarTime(:,iRoi,:,arousal))),'color',plotColors{arousal},'linestyle',linestyle);
        hold on
    end
    title('timepoint variability');
end
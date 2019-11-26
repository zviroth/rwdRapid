close all
clear all
tic
dataFolder = '~/data/rwdFmri/';
% subFolders = {'001920190214/'};
% roiNames = {'lh_V1_exvivo','rh_V1_exvivo'};
samplerate=500;
subFolders = {'007420190221/'};%nice pupil response
% subFolders = {'007520190318/'};%almost no pupil difference
% subFolders = {'008320190320/'};%nice pupil response
% subFolders = {'009020190325/'};%no pupil size difference
% subFolders = {'001920190403/'};%very bad eye tracking
% subFolders = {'007420190221/', '007520190318/', '008320190320/', '009020190325/', '001920190403/'};
% subFolders = {'005520190411', '008020190411', '007420190221/', '007520190318/', '008320190320/', '009020190325/', '001920190403/'};
subFolders = {'008620190410', '005520190411', '008020190411', '007420190221/', '007520190318/', '008320190320/', '009020190325/'};
% subFolders = {'008620190410', '005520190411'};
% subFolders = {'008620190410', '005520190411'};
subFolders = {'010420190805','009020190515','009720190603','001920190515', '007420190529', '008320190523',  '008620190410', '005520190411', '007520190715/', '008020190509'};
subFolders = {'010420190805','009020190515','009720190603', '007420190529', '008320190523',  '008620190410', '005520190411', '007520190715/'};

% subFolders = {'007420190221_OC/'};
% physioNames = {[
% roiNames = {'left','right'};
% roiNames = {'l_v1','r_v1'};
% roiNames = {'l_v1_ecc4','r_v1_ecc4'};
% roiNames = {'l_v1','r_v1','l_v1_ecc4','r_v1_ecc4'};
roiNames = {'lV1_eccen4','rV1_eccen4'};
% roiNames = {'lV1_eccen4'};
% roiNames = {'lV1_eccen0to1','rV1_eccen0to1','lV1_eccen1to4','rV1_eccen1to4','lV1_eccen4to10','rV1_eccen4to10','lV1_eccen10','rV1_eccen10'};
frames=170;
junkedFrames = 10;
TR=1.5;

physioSampleRate = 50; %Hz
numSubs = length(subFolders);
numRois = length(roiNames);
curFolder = pwd;
onlyCorrect=0;
cd(dataFolder);
mrQuit;
if exist('v')
    deleteView(v);
end
rows = numRois;
cols=numSubs;
trialsPerRun=16;
trialLength = 10;
nVolumes = trialsPerRun*trialLength;
numContrasts = 2;
numFreqs = 5;
plotColors = {[0 0 1], [1 0 0], [0 1 0], [0.5 1 0.2]};
plotStyles = {'-','--',':','-.','-','--',':','-.'};
clear d concatInfo sumResponse subResponse roiMeanTseries meanResponse roiTC subRunResponse trialCorrectness trialResponse trialRT propCorrect
clear contrast spatFreq nullTrial contrastConcat spatFreqConcat nullTrialConcat designMatContrast designMatFreq nullDeconv deconvMatNull meanDeconvNull meanDeconvNullConstant
clear taskResp taskPrediction meanPupil
clear ecg ecgPeaks ecgPeaksDiff ecgPeaksAmp resp respPeaks respPeaksDiff respPeaksAmp
clear designMatContrastFreq designMatNull taskBaseline taskTimes rwdPupil
clear sacRun pupilRun startTimes endTimes eyetrackerTime

L = 100;
filter = ones(L,1);
% filter = [1:L/2 L/2:-1:1];
filter = filter./sum(filter);


for iSub = 1:numSubs
    clear e runEyeSize eyeVel sac
    cd(subFolders{iSub});
    v=newView;
    % switch to the concatenation group
    v = viewSet(v, 'curGroup', 'Concatenation');
    nScans = viewGet(v, 'nscans');
    
    for iScan = 1:nScans%2 concatenations, 1 for each reward type
        d{iScan} = viewGet(v, 'd', iScan);
        %             concatInfo{iSub,iScan} = viewGet(v, 'concatInfo', iScan);
        s = viewGet(v, 'stimfile', iScan);
        
        rwdType = s{1}.myscreen.stimulus.rewardType;
        if strcmp(rwdType, 'H')
            rwd = 1;
        elseif strcmp(rwdType, 'L')
            rwd = 2;
        else
            disp('wtf');
            keyboard
        end
        concatRwdTypeNum(iScan) = rwd;
        for r=1:length(s)
            trialCorrectness{iSub,rwd}(r,:) = s{r}.task{1}{1}.correctness;
            trialResponse{iSub,rwd,r} = s{r}.task{1}{1}.response;
            propCorrect{iSub,rwd}(r) = sum(trialCorrectness{iSub,rwd}(r,:)) / length(trialCorrectness{iSub,rwd}(r,:));
            stimfile = s{r}.filename;
            e{r} = myGetTaskEyeTraces(stimfile, 'removeBlink=3');
            startTimes{iSub,rwd,r} = e{r}.complete.startTime;
            endTimes{iSub,rwd,r} = e{r}.complete.endTime;
            eyetrackerTime{iSub,rwd,r} = e{r}.complete.time;
            runEyeSize{rwd}(r,:) = size(e{r}.eye.pupil);
            nullTrialsRun{iSub,rwd,r} = s{r}.task{2}{1}.randVars.nullTrial(1:trialsPerRun+1);
        end
        trialLengthEye(rwd) = max(runEyeSize{rwd}(:,2));
        numTrials(rwd) = sum(runEyeSize{rwd}(:,1));
        sacRwd{iSub,rwd} = zeros(numTrials(rwd), trialLengthEye(rwd));
        trialCounter=0;
        numRuns(iSub,rwd) = length(s);
        for r=1:length(s)
            
            eyePos = [e{r}.complete.xPos' e{r}.complete.yPos'];% time X 2, x and y
            if sum(~isnan(eyePos(:)))~=0
                eyevel{r} = vecvel(eyePos,samplerate, 2);
                sacOutput = microsacc(eyePos, eyevel{r}, 1, 5);
                sacOnset = sacOutput(:,1);
                sacRun{iSub,rwd,r} = sacOnset;%saccade onset
                sacRunTC = zeros(size(eyePos,1),1);
                sacRunTC(sacOnset) = ones;
                sacRunSmooth{iSub,rwd,r} = conv(sacRunTC, filter, 'same');
                pupilRun{iSub,rwd,r} = e{r}.complete.pupil;%e{r}.eye.pupil is longer because of zero padding trials;
                
                for itrial=1:length(e{r}.complete.startTime)
                    trialCounter=trialCounter+1;
                    sac{r,itrial} = sacOnset(sacOnset >= e{r}.complete.startTime(itrial) & sacOnset<=e{r}.complete.endTime(itrial));%saccade onsets for this trial
                    sac{r,itrial} = sac{r,itrial} - e{r}.complete.startTime(itrial) + 1;%time relative to trial
                    sacTrial = zeros(1,trialLengthEye(rwd));
                    sacTrial(sac{r,itrial}) = ones;
                    sacRwd{iSub,rwd}(trialCounter,:) = sacTrial;
                    
                    sacRunSmoothTrial{r,itrial} = sacRunSmooth{iSub,rwd,r}(e{r}.complete.startTime(itrial):e{r}.complete.endTime(itrial));
                    %                     sacTrial = zeros(1,trialLengthEye(rwd));
                    sacTrial = NaN(1,trialLengthEye(rwd));
                    sacTrial(1:length(sacRunSmoothTrial{r,itrial})) = sacRunSmoothTrial{r,itrial};
                    sacRunSmoothTrialZeroFilled{iSub,rwd}(trialCounter,:) = sacTrial;
                    
                end
                
            end
            
        end
        
        binnedSac{iSub,rwd} = sum(sacRwd{iSub,rwd});
        rwdPupil{iSub,rwd} = NaN(numTrials(rwd), trialLengthEye(rwd));
        nullTrials{iSub,rwd} = NaN(numTrials(rwd),1);
        trialCounter=0;
        for r=1:length(s)
            rwdPupil{iSub,rwd}(trialCounter+1:trialCounter+runEyeSize{rwd}(r,1), 1:runEyeSize{rwd}(r,2)) = e{r}.eye.pupil;
            nullTrials{iSub,rwd}(trialCounter+1:trialCounter+runEyeSize{rwd}(r,1)) = nullTrialsRun{iSub,rwd,r};
            trialCounter = trialCounter + runEyeSize{rwd}(r,1);
        end
        %plot mean pupil size
        meanPupil{iSub,rwd} = nanmean(rwdPupil{iSub,rwd})';
        meanNullPupil{iSub,rwd} = nanmean(rwdPupil{iSub,rwd}(nullTrials{iSub,rwd}==1,:))';
        meanStimPupil{iSub,rwd} = nanmean(rwdPupil{iSub,rwd}(nullTrials{iSub,rwd}==0,:))';
        
        cueTime = s{r}.fixStimulus.cueTime * samplerate;
        interTime = s{r}.fixStimulus.interTime * samplerate;
        responseTime = s{r}.fixStimulus.responseTime * samplerate;
        taskTimes(iSub,:) = [cueTime interTime cueTime responseTime];
    end
    for rwd=1:2
        smoothSac{iSub,rwd} = nanmean(sacRunSmoothTrialZeroFilled{iSub,rwd});
        
    end
    deleteView(v);
    
    
    %     %smooth out the microsaccade data
    %     for rwd=1:2
    %         smoothSac{iSub,rwd} = conv(binnedSac{iSub,rwd}, filter, 'same');
    %     end
    
    
    %% return to home directory
    cd(dataFolder);
end

save([dataFolder 'rwdRapidEyeData.mat'], 'dataFolder', 'subFolders', 'samplerate',  ...
    'physioSampleRate', 'trialsPerRun', 'trialLength', ...
    'numSubs', 'onlyCorrect', ...
    'sacRwd','binnedSac', 'smoothSac', ...
    'rwdPupil','meanPupil','taskTimes',...
    'sacRun','pupilRun', 'numRuns','startTimes','endTimes','eyetrackerTime',...
    'nullTrials','meanNullPupil','meanStimPupil',...
    'sacRunSmooth', 'sacRunSmoothTrialZeroFilled','sacRunSmoothTrial','L');



%%
toc
%%
figure(1); clf;
rows=2;
cols=numSubs;
for iSub=1:numSubs
    subplot(rows,cols,iSub)
    for rwd=1:2
        plot(smoothSac{iSub,rwd})
        hold on
    end
end

for iSub=1:numSubs
    subplot(rows,cols,iSub+cols)
    for rwd=1:2
        plot(binnedSac{iSub,rwd})
        hold on
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = plotRwdStimSegments(meanPupil, allTimes)

itime = 1;
startT=1;
endT = startT+allTimes(itime);
plot(startT:endT-1,meanPupil(startT:endT-1),'lineWidth',4,'color',[0 1 1]);
itime = itime+1;
startT=endT;
endT = startT+allTimes(itime);
plot(startT:endT-1,meanPupil(startT:endT-1),'lineWidth',4,'color',[0.1 0.7 0.6]);
itime = itime+1;
startT=endT;
endT = startT+allTimes(itime);
plot(startT:endT-1,meanPupil(startT:endT-1),'lineWidth',4,'color',[0 1 1]);
itime = itime+1;
startT=endT;
endT = startT+allTimes(itime);
plot(startT:endT-1,meanPupil(startT:endT-1),'lineWidth',4,'color',[1 1 0.0]);

end
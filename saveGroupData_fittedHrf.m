dataFolder = '~/data/rwdFmri/';
% subFolders = {'001920190214/'};
% roiNames = {'lh_V1_exvivo','rh_V1_exvivo'};
samplerate=500;
subFolders = {'007420190221/'};%nice pupil response
% subFolders = {'007520190318/'};%almost no pupil difference
% subFolders = {'008320190320/'};%nice pupil response
% subFolders = {'009020190325/'};%no pupil size difference
% subFolders = {'001920190403/'};%very bad eye tracking
subFolders = {'008620190410', '005520190411', '008020190411', '007420190221/', '007520190318/', '008320190320/', '009020190325/', '001920190403/'};
% subFolders = {'008620190410', '005520190411'};


% subFolders = {'007420190221_OC/'};
% physioNames = {[
% roiNames = {'left','right'};
% roiNames = {'l_v1','r_v1'};
% roiNames = {'l_v1_ecc4','r_v1_ecc4'};
% roiNames = {'l_v1','r_v1','l_v1_ecc4','r_v1_ecc4'};
roiNames = {'lV1_eccen4','rV1_eccen4','lV1_eccen8','rV1_eccen8'};
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
clear contrast spatFreq nullTrial contrastConcat spatFreqConcat nullTrialConcat designMatContrast designMatFreq nullDeconv deconvMatNull meanDeconvNull
clear taskResp taskPrediction meanPupil
clear ecg ecgPeaks ecgPeaksDiff ecgPeaksAmp resp respPeaks respPeaksDiff respPeaksAmp
clear designMatContrastFreq designMatNull taskBaseline meanDeconvFreq meanDeconvContrastFreq

deconvLength = 10;

curSubplot=1;

for iSub = 1:numSubs
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
            rwdTypeNum = 1;
        elseif strcmp(rwdType, 'L')
            rwdTypeNum = 2;
        else
            disp('wtf');
            keyboard
        end
        
        for iRoi=1:numRois
            roiTC{iSub,iRoi,rwdTypeNum} = loadROITSeries(v, roiNames{iRoi}, iScan, [], 'keepNAN',true);
            concatInfo{iSub,rwdTypeNum} = viewGet(v, 'concatInfo', iScan);
            nvox = roiTC{iSub,iRoi,rwdTypeNum}.n;
            temp = reshape(roiTC{iSub,iRoi,rwdTypeNum}.tSeries,nvox,trialLength,[]);
            voxMeanRespRwd{iSub,iRoi,rwdTypeNum} = mean(temp,3);%[nvox,10]
            %ROI average trial
            trialTC{iSub,iRoi}(rwdTypeNum,:) = nanmean(reshape(nanmean(roiTC{iSub,iRoi,rwdTypeNum}.tSeries), trialLength, []),2);%average across voxels
            
        end

        
        %extract conditions for all scans within concatenation
        deconvMatNull{rwdTypeNum} = [];
        deconvMatContrast{rwdTypeNum} = [];
        deconvMatFreq{rwdTypeNum} = [];
        deconvMatContrastFreq{rwdTypeNum} = [];
        
                %create  vector with a 1 for each trial in a single run, to be convolved with the
        %task response (null trials). This is identical for all runs!
        taskSmallMat = ones(1,trialsPerRun);
        taskMat = upsample(taskSmallMat,trialLength);%add zeros in between
        clear contrastMat sfMat contrastSfMat nullTrialMat
        for r=1:length(s)%create design matrices per run
            nullTrials = s{r}.task{2}{1}.randVars.nullTrial(2:trialsPerRun+1);%skipping the first because of junked frames
            runContrastEvents = s{r}.task{2}{1}.randVars.contrast(2:trialsPerRun+1);
            runFreqEvents = s{r}.task{2}{1}.randVars.spatFreq(2:trialsPerRun+1);
            runContrastEvents = runContrastEvents.*(1-nullTrials);
            runFreqEvents = runFreqEvents.*(1-nullTrials);
            
            %contrast
            smallContrastMat = zeros(numContrasts,trialsPerRun);
            for c=1:numContrasts
                smallContrastMat(c,runContrastEvents==c) = ones;%1=low, 2=high
            end
            contrastMat{r} = upsample(smallContrastMat,trialLength);%number per volume

            %spatial frequency
            smallFreqMat = zeros(numFreqs,trialsPerRun);%number per trial
            for c=1:numFreqs
                smallFreqMat(c,runFreqEvents==c) = ones;
            end
            freqMat{r} = upsample(smallFreqMat,trialLength);%number per volume

            %spatial frequency & contrast!
            smallContrastSfMat = zeros(numFreqs*numContrasts,trialsPerRun);%number per trial
            for c=1:numFreqs*numContrasts
                curCond = (runContrastEvents-1)*numFreqs+runFreqEvents;
                smallContrastSfMat(c,curCond==c) = ones;
            end
            contrastSfMat{r} = upsample(smallContrastSfMat,trialLength);%number per volume

            %null vs stimulated trials
            smallNullMat = zeros(2,trialsPerRun);%number per trial
            for c=0:1
                smallNullMat(c+1,nullTrials==c) = ones; %0=stim, 1=null
            end
            nullTrialMat{r} = upsample(smallNullMat,trialLength);%number per volume
            
            %nullTrial deconvolution matrix
            deconvMatNullRun = zeros((size(nullTrialMat{r},1)-1)*deconvLength,size(nullTrialMat{r},2));
            for c=1:2
                for j=1:deconvLength
                    deconvMatNullRun((c-1)*deconvLength+j,:) = [zeros(1,j-1) nullTrialMat{r}(c,1:end-j+1)];%shift and pad with zeros
                end
            end
            deconvMatNullRun(end+1,:) = ones;
            deconvMatNull{rwdTypeNum} = [deconvMatNull{rwdTypeNum} deconvMatNullRun];
            
            %contrast deconvolution matrix
            deconvMatContrastRun = zeros((size(contrastMat{r},1)-1)*deconvLength,size(contrastMat{r},2));
            for c=1:numContrasts
                for j=1:deconvLength
                    deconvMatContrastRun((c-1)*deconvLength+j,:) = [zeros(1,j-1) contrastMat{r}(c,1:end-j+1)];%shift and pad with zeros
                end
            end
            deconvMatContrastRun(end+1,:) = ones;
            deconvMatContrast{rwdTypeNum} = [deconvMatContrast{rwdTypeNum} deconvMatContrastRun];
            
            %sf deconvolution matrix
            deconvMatFreqRun = zeros((size(freqMat{r},1)-1)*deconvLength,size(freqMat{r},2));
            for c=1:numFreqs
                for j=1:deconvLength
                    deconvMatFreqRun((c-1)*deconvLength+j,:) = [zeros(1,j-1) freqMat{r}(c,1:end-j+1)];%shift and pad with zeros
                end
            end
            deconvMatFreqRun(end+1,:) = ones;
            deconvMatFreq{rwdTypeNum} = [deconvMatFreq{rwdTypeNum} deconvMatFreqRun];
            
            %sf & contrast deconvolution matrix
            deconvMatContrastFreqRun = zeros((size(contrastSfMat{r},1)-1)*deconvLength,size(contrastSfMat{r},2));
            for c=1:numFreqs*numContrasts
                for j=1:deconvLength
                    deconvMatContrastFreqRun((c-1)*deconvLength+j,:) = [zeros(1,j-1) contrastSfMat{r}(c,1:end-j+1)];%shift and pad with zeros
                end
            end
            deconvMatContrastFreqRun(end+1,:) = ones;
            deconvMatContrastFreq{rwdTypeNum} = [deconvMatContrastFreq{rwdTypeNum} deconvMatContrastFreqRun];

            
        end
    end
    
    
    %get mean response per voxel - to be used instead of canonical HRF
    for iRoi=1:numRois
        voxMeanResp{iSub,iRoi} = (voxMeanRespRwd{iSub,iRoi,1}+voxMeanRespRwd{iSub,iRoi,2})/2;
    end
            
    for rwdTypeNum=1:2
        for iRoi=1:numRois
            %create design matrices for each voxel, using mean response
            %in place of canonical hrf
            nvox = roiTC{iSub,iRoi,rwdTypeNum}.n;
            for i=1:nvox
                tc = roiTC{iSub,iRoi,rwdTypeNum}.tSeries(i,:);
                hrf = voxMeanResp{iSub,iRoi}(i,:);
                hrf = zscore(hrf);
                
                designMatContrastVox=[];
                designMatFreqVox=[];
                designMatContrastFreqVox=[];
                designMatNullVox=[];
                for r=1:length(s)%create design matrices per run. at the beginning of each run the predictors must be zeroed
                    
                    %contrast
                    for c=1:numContrasts
                        temp = conv(contrastMat{r}(c,:),hrf);
                        designMatContrastVoxRun(c,:) =  temp(1:trialsPerRun*trialLength);
                    end
                    designMatContrastVoxRun(numContrasts+1,:) = ones;
                    designMatContrastVox = [designMatContrastVox designMatContrastVoxRun];
                    
                    
                    %sf
                    for c=1:numFreqs
                        temp = conv(freqMat{r}(c,:),hrf);
                        designMatFreqVoxRun(c,:) =  temp(1:trialsPerRun*trialLength);
                    end
                    designMatFreqVoxRun(numFreqs+1,:) = ones;
                    designMatFreqVox = [designMatFreqVox designMatFreqVoxRun];
                    
                    
                    %contrast & sf
                    for c=1:numFreqs*numContrasts
                        temp = conv(contrastSfMat{r}(c,:),hrf);
                        designMatContrastFreqVoxRun(c,:) =  temp(1:trialsPerRun*trialLength);
                    end
                    designMatContrastFreqVoxRun(numFreqs*numContrasts+1,:) = ones;
                    designMatContrastFreqVox = [designMatContrastFreqVox designMatContrastFreqVoxRun];
                    
                    
                    %null vs stim
                    for c=1:2
                        temp = conv(nullTrialMat{r}(c,:),hrf);
                        designMatNullVoxRun(c,:) =  temp(1:trialsPerRun*trialLength);
                    end
                    designMatNullVoxRun(2+1,:) = ones;
                    designMatNullVox = [designMatNullVox designMatNullVoxRun];
                    
                end
                
                %perform GLM to get betas
                contrastBetasWithTask{iSub,iRoi,rwdTypeNum}(:,i) = designMatContrastVox'\tc';
                sfBetasWithTask{iSub,iRoi,rwdTypeNum}(:,i) = designMatFreqVox'\tc';
                sfContrastBetasWithTask{iSub,iRoi,rwdTypeNum}(:,i)= designMatContrastFreqVox'\tc';
                nullBetas{iSub,iRoi,rwdTypeNum}(:,i) = designMatNullVox'\tc';

%                 %deconvolution with stim/null predictors
%                 nullDeconv{iRoi,rwdTypeNum}(:,i) = deconvMatNull{iScan}'\tc';
                
            end
            %deconvolution with stim/null predictors
            nullDeconv{iSub,iRoi,rwdTypeNum} = deconvMatNull{rwdTypeNum}'\roiTC{iSub,iRoi,rwdTypeNum}.tSeries';
                
            %deconvolution with contrast predictors
            contrastDeconv{iSub,iRoi,rwdTypeNum} = deconvMatContrast{rwdTypeNum}'\roiTC{iSub,iRoi,rwdTypeNum}.tSeries';
            
            %deconvolution with frequency predictors
            freqDeconv{iSub,iRoi,rwdTypeNum} = deconvMatFreq{rwdTypeNum}'\roiTC{iSub,iRoi,rwdTypeNum}.tSeries';
            
            %deconvolution with frequency & contrast predictors
            contrastFreqDeconv{iSub,iRoi,rwdTypeNum} = deconvMatContrastFreq{rwdTypeNum}'\roiTC{iSub,iRoi,rwdTypeNum}.tSeries';

            
            %get the null response
            taskResp{iSub,iRoi}(rwdTypeNum,:,:) = nullDeconv{iSub,iRoi,rwdTypeNum}(deconvLength+1:2*deconvLength,:); %for each voxel
            taskBaseline{iSub,iRoi}(rwdTypeNum,:) = nullDeconv{iSub,iRoi,rwdTypeNum}(2*deconvLength+1,:); %for each voxel
            
            %get mean deconvolved null/stim responses
            for c=1:2
                for t=1:deconvLength
                    meanDeconvNull(iSub,iRoi,rwdTypeNum,c,t) = nanmean(nullDeconv{iSub,iRoi,rwdTypeNum}((c-1)*deconvLength+t,:));%mean over voxels
                end
            end
            meanDeconvNullConstant(iSub,iRoi,rwdTypeNum) = nanmean(nullDeconv{iSub,iRoi,rwdTypeNum}(2*deconvLength+1,:));
            

            %get mean deconvolved contrast responses
            for c=1:numContrasts
                for t=1:deconvLength
                    meanDeconvContrast(iSub,iRoi,rwdTypeNum,c,t) = nanmean(contrastDeconv{iSub,iRoi,rwdTypeNum}((c-1)*deconvLength+t,:));%mean over voxels
                end
            end
            meanDeconvContrastConstant(iSub,iRoi,rwdTypeNum) = nanmean(contrastDeconv{iSub,iRoi,rwdTypeNum}(numContrasts*deconvLength+1,:));

            %get mean deconvolved frequency responses
            for c=1:numFreqs
                for t=1:deconvLength
                    meanDeconvFreq(iSub,iRoi,rwdTypeNum,c,t) = nanmean(freqDeconv{iSub,iRoi,rwdTypeNum}((c-1)*deconvLength+t,:));%mean over voxels
                end
            end
            meanDeconvFreqConstant(iSub,iRoi,rwdTypeNum) = nanmean(freqDeconv{iSub,iRoi,rwdTypeNum}(numFreqs*deconvLength+1,:));

            %get mean deconvolved contrast & frequency responses
            for c=1:numFreqs*numContrasts
                for t=1:deconvLength
                    meanDeconvContrastFreq(iSub,iRoi,rwdTypeNum,c,t) = nanmean(contrastFreqDeconv{iSub,iRoi,rwdTypeNum}((c-1)*deconvLength+t,:));%mean over voxels
                end
            end
            meanDeconvContrastFreqConstant(iSub,iRoi,rwdTypeNum) = nanmean(contrastFreqDeconv{iSub,iRoi,rwdTypeNum}(numFreqs*numContrasts*deconvLength+1,:));

            
            
            %convole with task events
            taskPrediction{iSub,iRoi,rwdTypeNum}=[];
            clear taskPredictionRun
            for r=1:length(s)%per run
                for i=1:nvox
                    temp = conv(taskMat, taskResp{iSub,iRoi}(rwdTypeNum,:,i));
                    taskPredictionRun(i,:) = temp(1:trialsPerRun*trialLength);
                end
                taskPrediction{iSub,iRoi,rwdTypeNum} = [taskPrediction{iSub,iRoi,rwdTypeNum} taskPredictionRun];
            end
            
            %subtract the task response
            tcMinusTask{iSub,iRoi,rwdTypeNum} = roiTC{iSub,iRoi,rwdTypeNum}.tSeries - taskPrediction{iSub,iRoi,rwdTypeNum};
            %             tcMinusTask{iRoi,rwdTypeNum} = zscore(tcMinusTask{iRoi,rwdTypeNum},0,2);
            
            
            %recompute all betas, after null response is removed
            for i=1:nvox
                tc = tcMinusTask{iSub,iRoi,rwdTypeNum}(1,:);
                hrf = voxMeanResp{iSub,iRoi}(i,:);%should we recompute the hrf??
                hrf = zscore(hrf);
                
                clear designMatContrastRun
                designMatContrastVox=[];
                designMatFreqVox=[];
                designMatContrastFreqVox=[];
                designMatNullVox=[];
                for r=1:length(s)%create design matrices per run
                    %contrast
                    for c=1:numContrasts
                        temp = conv(contrastMat{r}(c,:),hrf);
                        designMatContrastVoxRun(c,:) =  temp(1:trialsPerRun*trialLength);
                    end
                    designMatContrastVoxRun(numContrasts+1,:) = ones;
                    designMatContrastVox = [designMatContrastVox designMatContrastVoxRun];
                    %sf
                    for c=1:numFreqs
                        temp = conv(freqMat{r}(c,:),hrf);
                        designMatFreqVoxRun(c,:) =  temp(1:trialsPerRun*trialLength);
                    end
                    designMatFreqVoxRun(numFreqs+1,:) = ones;
                    designMatFreqVox = [designMatFreqVox designMatFreqVoxRun];
                    %contrast & sf
                    for c=1:numFreqs*numContrasts
                        temp = conv(contrastSfMat{r}(c,:),hrf);
                        designMatContrastFreqVoxRun(c,:) =  temp(1:trialsPerRun*trialLength);
                    end
                    designMatContrastFreqVoxRun(numFreqs*numContrasts+1,:) = ones;
                    designMatContrastFreqVox = [designMatContrastFreqVox designMatContrastFreqVoxRun];
                    %null vs stim
                    for c=1:2
                        temp = conv(nullTrialMat{r}(c,:),hrf);
                        designMatNullVoxRun(c,:) =  temp(1:trialsPerRun*trialLength);
                    end
                    designMatNullVoxRun(2+1,:) = ones;
                    designMatNullVox = [designMatNullVox designMatNullVoxRun];
                end
                %perform GLM to get betas
                contrastBetas{iSub,iRoi,rwdTypeNum}(:,i) = designMatContrastVox'\tc';
                sfBetas{iSub,iRoi,rwdTypeNum}(:,i) = designMatFreqVox'\tc';
                sfContrastBetas{iSub,iRoi,rwdTypeNum}(:,i)= designMatContrastFreqVox'\tc';
                nullDeconvBetas{iSub,iRoi,rwdTypeNum}(:,i) = designMatNullVox'\tc';
            end
            
        end
    end
    deleteView(v);
    

    %% return to home directory
    cd('..');
end
%%
save([dataFolder 'rwdRapidData_fittedHrf.mat'], 'dataFolder', 'subFolders', 'samplerate', 'roiNames', ...
    'frames', 'junkedFrames', 'TR', 'trialsPerRun', 'trialLength', 'nVolumes',...
    'numContrasts', 'numFreqs', 'numSubs', 'numRois', 'onlyCorrect', ...
    'concatInfo', 'contrastBetasWithTask', 'sfBetasWithTask', 'sfContrastBetasWithTask', ...
    'taskResp','taskBaseline', 'nullBetas', 'contrastBetas','sfBetas','sfContrastBetas', ...
    'trialTC','meanDeconvNull','meanDeconvNullConstant','meanDeconvFreqConstant','nullDeconvBetas',...
    'meanDeconvContrastConstant', 'meanDeconvContrast', ...
    'meanDeconvFreqConstant', 'meanDeconvFreq',...
    'meanDeconvContrastFreq', 'meanDeconvContrastFreqConstant');


%%
%%%%%%%%%%%%%%%%%%%%%%%%%
% increases the sample rate of x by inserting n - 1 zeros between samples. If x is a matrix, the function treats each column as a separate sequence.
% inserts additional n-1 zeros after the last sample.
%%%%%%%%%%%%%%%%%%%%%%%%%%
function bigMatrix = upsample(smallMatrix,volumesPerTrial)
numTrials = size(smallMatrix,2);
numConditions = size(smallMatrix,1);
bigMatrix = zeros(numConditions, volumesPerTrial*numTrials);
bigMatrix(:,1:volumesPerTrial:end) = smallMatrix;
end

%%
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
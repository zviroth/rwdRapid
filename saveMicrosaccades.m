close all;
clear all;
dataFolder = '/Volumes/MH02086153MACDT-Drobo/rwdFmri/';
plotColors = {[0 0 1], [1 0 0], [0 1 0], [0.5 1 0.2]};
plotStyles = {'-','--',':','-.','-','--',':','-.'};

load([dataFolder 'rwdRapidEyeData.mat'], 'dataFolder', 'subFolders', 'samplerate',  ...
    'trialsPerRun', 'trialLength', ...
    'numSubs', 'onlyCorrect', ...
    'sacRwd','binnedSac', 'smoothSac', ...
    'rwdPupil','meanPupil','taskTimes',...
    'sacRun','pupilRun', 'numRuns','startTimes','endTimes','eyetrackerTime',...
    'nullTrials','nullPupil','stimPupil','meanNullPupil','meanStimPupil',...
    'sacRunSmooth', 'sacRunSmoothTrialZeroFilled','sacRunSmoothTrial','L');

clear sacPupilFunc sacPupilFuncRun resid residRun subMeanPupilLength sacMat sacPupilRun
TR=1.5;
% numSubs=2;

tic
% sacRun{iSub,rwdTypeNum,r}; %saccade onset
% pupilRun{iSub,rwdTypeNum,r}; %pupil size
deconvLength= 500; %in eyetracker samples
for iSub=1:numSubs
%     numRuns(iSub) = size(pupilRun,3);%what if there is a different number of high and low runs?
    concatSac=[];
    concatPupil=[];
    concatSacMat=[];
    for rwd=1:2

        for r=1:numRuns(iSub,rwd)

%             %detrend, convert to % signal change
%             goodtimepoints = ~isnan(pupilRun{iSub,rwd,r});
%             temp = pupilRun{iSub,rwd,r}(goodtimepoints);
%             temp = percentTSeries(temp')';
%             pupilRun{iSub,rwd,r}(goodtimepoints) = temp;


            runLength = length(pupilRun{iSub,rwd,r}(:));
            if runLength==0
                runLength = (trialsPerRun+1)*samplerate*trialLength*TR;
                pupilRun{iSub,rwd,r} = NaN(1,runLength);
            end
            sacVec = zeros(1,runLength);
            sacVec(sacRun{iSub,rwd,r}) = ones;%one for each microsaccade onset
            
            sacMat{iSub,rwd,r} = zeros(deconvLength,runLength);
            for j=1:deconvLength
                sacMat{iSub,rwd,r}(j,:) = [zeros(1,j-1) sacVec(1:end-j+1)];%shift and pad with zeros
            end
            sacMat{iSub,rwd,r}(deconvLength+1,:) = ones;
            
            goodTimepoints = ~isnan(pupilRun{iSub,rwd,r});
            %single run pupil function, for illustration purposes...
            sacPupilFuncRun{iSub}(rwd,r,:) = sacMat{iSub,rwd,r}(:,goodTimepoints)'\pupilRun{iSub,rwd,r}(goodTimepoints)';%deconvolution
            %             residRun{iSub,rwd,r}(:) = pupilRun{iSub,rwd,r}' - sacMat'*squeeze(sacPupilFuncRun{iSub}(rwd,r,:));%residual timecourse
            
            concatSacMat = [concatSacMat sacMat{iSub,rwd,r}];
            concatSac = [concatSac; sacRun{iSub,rwd,r}+length(concatPupil)];
            concatPupil = [concatPupil pupilRun{iSub,rwd,r}];
            
            
        end
    end
    goodTimepoints = ~isnan(concatPupil);
    sacPupilFunc(iSub,:) = concatSacMat(:,goodTimepoints)'\concatPupil(goodTimepoints)';
    
    %set the function integral to zero
    offset = sum(sacPupilFunc(iSub,1:end-1))./length(sacPupilFunc(iSub,1:end-1));
    sacPupilFunc(iSub,1:end-1) = sacPupilFunc(iSub,1:end-1)-offset;
    
%      sacPupilFunc(iSub,end) =  sacPupilFunc(iSub,end) + offset;
    
%     sacPupilFunc(iSub,:) = sacPupilFunc(iSub,:) - sacPupilFunc(iSub,1);%subtract the first time point of the pupil-ms function from the entire function, so that it starts at zero
    
    
%     resid{iSub} = concatPupil' - concatSacMat(1:end-1,:)'*sacPupilFunc(iSub,1:end-1)';
%     resid{iSub} = resid{iSub} + sacPupilFunc(iSub,end);%add the constant
    
    %perform regression for each run, using the concatenated data's ms
    %pupil function
    for rwd=1:2
        for r=1:numRuns(iSub,rwd)
            sacPupilRun{iSub,rwd,r} = sacMat{iSub,rwd,r}(1:end-1,:)'*squeeze(sacPupilFunc(iSub,1:end-1))';%ms pupil function convolved with all microsaccade onsets
            residRun{iSub,rwd,r}(:) = pupilRun{iSub,rwd,r}' - sacMat{iSub,rwd,r}(1:end-1,:)'*squeeze(sacPupilFunc(iSub,1:end-1))';%residual timecourse. 
%             residRun{iSub,rwd,r}(:) = residRun{iSub,rwd,r}(:) + sacPupilFunc(iSub,end); %add the constant 
        end
    end
    
    %cut residuals into trials, and average
    for rwd=1:2
        residRwdPupil{iSub,rwd} = NaN(size(rwdPupil{iSub,rwd}));%numTrials X maxTrialLength
        regressRwdPupil{iSub,rwd} = NaN(size(rwdPupil{iSub,rwd}));
        trialCounter=0;%trial counter across runs
       %get residual for each run
       for r=1:numRuns(iSub,rwd)
           for itrial=1:size(startTimes{iSub,rwd,r},2)
               
               thisTrialTimes = startTimes{iSub,rwd,r}(itrial):endTimes{iSub,rwd,r}(itrial);
               [~, t1] = min(abs(eyetrackerTime{iSub,rwd,r}(1:length(pupilRun{iSub,rwd,r}))-thisTrialTimes(1)));
               t1=t1(1);
               [~, t2] = min(abs(eyetrackerTime{iSub,rwd,r}(1:length(pupilRun{iSub,rwd,r}))-thisTrialTimes(end)));
               t2=t2(1);
               
               trialResidual = residRun{iSub,rwd,r}(t1:t2);
%                trialResidual = pupilRun{iSub,rwd,r}(t1:t2); % TESTING

               
               
%                trialPupil = residRun{iSub,rwd,r}(startTimes{iSub,rwd,r}(itrial):endTimes{iSub,rwd,r}(itrial)-1);%endTime is startTime of next trial
               trialCounter = trialCounter+1;
               residRwdPupil{iSub,rwd}(trialCounter,1:length(trialResidual)) = trialResidual;
               
               trialRegressed = sacPupilRun{iSub,rwd,r}(t1:t2);
               regressRwdPupil{iSub,rwd}(trialCounter,1:length(trialResidual)) = trialRegressed;%this is regressed out, the pupil modulation from microsaccades
           end
       end
       
       %average across trials
       residMeanPupil{iSub,rwd} = nanmean(residRwdPupil{iSub,rwd});
       meanResidNullPupil{iSub,rwd} = nanmean(residRwdPupil{iSub,rwd}(nullTrials{iSub,rwd}==1,:))';
       meanResidStimPupil{iSub,rwd} = nanmean(residRwdPupil{iSub,rwd}(nullTrials{iSub,rwd}==0,:))';
       
       regressMeanPupil{iSub,rwd} = nanmean(regressRwdPupil{iSub,rwd});
    end
    


    l = min(length(meanPupil{iSub,1}), length(meanPupil{iSub,2}));
    subMeanPupilLength(iSub) = l;
end

%%
allMeanPupil = zeros(2,min(subMeanPupilLength));
allResidMeanPupil = zeros(2,min(subMeanPupilLength));
sumBinnedSac = zeros(2,min(subMeanPupilLength));
sumSmoothSac = zeros(2,min(subMeanPupilLength));
allMeanStimPupil = zeros(2,min(subMeanPupilLength));
allResidStimPupil = zeros(2,min(subMeanPupilLength));
allMeanNullPupil = zeros(2,min(subMeanPupilLength));
allResidNullPupil = zeros(2,min(subMeanPupilLength));
allRegressMeanPupil = zeros(2,min(subMeanPupilLength));

for iSub=1:numSubs
    for rwd=1:2
        allMeanPupil(rwd,:) = allMeanPupil(rwd,:) + meanPupil{iSub,rwd}(1:min(subMeanPupilLength))';
        allResidMeanPupil(rwd,:) = allResidMeanPupil(rwd,:) + residMeanPupil{iSub,rwd}(1:min(subMeanPupilLength));
        
        sumBinnedSac(rwd,:) =  sumBinnedSac(rwd,:) + binnedSac{iSub,rwd}(1:min(subMeanPupilLength));
        if length(smoothSac{iSub,rwd})>1
            sumSmoothSac(rwd,:) =  sumSmoothSac(rwd,:) + smoothSac{iSub,rwd}(1:min(subMeanPupilLength));
        end
        
        allMeanStimPupil(rwd,:) = allMeanStimPupil(rwd,:) + meanStimPupil{iSub,rwd}(1:min(subMeanPupilLength))';
        allResidStimPupil(rwd,:) = allResidStimPupil(rwd,:) + meanResidStimPupil{iSub,rwd}(1:min(subMeanPupilLength))';
        
        allMeanNullPupil(rwd,:) = allMeanNullPupil(rwd,:) + meanNullPupil{iSub,rwd}(1:min(subMeanPupilLength))';
        allResidNullPupil(rwd,:) = allResidNullPupil(rwd,:) + meanResidNullPupil{iSub,rwd}(1:min(subMeanPupilLength))';
        
        
        allRegressMeanPupil(rwd,:) = allRegressMeanPupil(rwd,:) + regressMeanPupil{iSub,rwd}(1:min(subMeanPupilLength));
    end
end
allMeanPupil = allMeanPupil./numSubs;
allResidMeanPupil = allResidMeanPupil./numSubs;
sumBinnedSac = sumBinnedSac./numSubs;
sumSmoothSac = sumSmoothSac./numSubs;
allMeanStimPupil = allMeanStimPupil./numSubs;
allResidStimPupil = allResidStimPupil./numSubs;
allMeanNullPupil = allMeanNullPupil./numSubs;
allResidNullPupil = allResidNullPupil./numSubs;
        
allRegressMeanPupil = allRegressMeanPupil./numSubs;


%%
save([dataFolder 'rwdRapid_microsaccade.mat'], 'dataFolder', 'subFolders', 'samplerate',  ...
    'trialsPerRun', 'trialLength', ...
    'numSubs', 'onlyCorrect', ...
    'sacRwd','binnedSac', 'smoothSac', ...
    'rwdPupil','meanPupil','taskTimes',...
    'sacRun','pupilRun', 'numRuns','startTimes','endTimes','eyetrackerTime',...
    'nullTrials','nullPupil','stimPupil','meanNullPupil','meanStimPupil',...
    'sacRunSmooth', 'sacRunSmoothTrialZeroFilled','sacRunSmoothTrial','L',...
    'smoothSac', 'sacPupilFuncRun', 'sacPupilFunc','sumBinnedSac','sumSmoothSac','allMeanPupil', ...
    'allResidMeanPupil', 'allMeanStimPupil','allResidStimPupil','allMeanNullPupil','allResidNullPupil','allRegressMeanPupil',...
    'concatSac','residRwdPupil','regressRwdPupil','residRun','sacPupilRun','residMeanPupil','meanResidNullPupil','meanResidNullPupil');
toc

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
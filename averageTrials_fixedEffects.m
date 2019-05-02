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

% subFolders = {'007420190221_OC/'};
% physioNames = {[
% roiNames = {'left','right'};
% roiNames = {'l_v1','r_v1'};
% roiNames = {'l_v1_ecc4','r_v1_ecc4'};
% roiNames = {'l_v1','r_v1','l_v1_ecc4','r_v1_ecc4'};
% roiNames = {'lV1_eccen4','rV1_eccen4','lV2_eccen4','rV2_eccen4','lV3_eccen4','rV3_eccen4'};
roiNames = {'lV1_eccen4','rV1_eccen4'};
roiNames = {'lV1_eccen8','rV1_eccen8'};
% roiNames = {'lV1_eccen0to1','rV1_eccen0to1','lV1_eccen1to4','rV1_eccen1to4','lV1_eccen4to10','rV1_eccen4to10','lV1_eccen10','rV1_eccen10'};
% roiNames = {'rV1_eccen0to1','rV1_eccen1to4','rV1_eccen4to10','rV1_eccen10'};
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
clear contrastTrials contrastAmp freqTrials freqAmp allTrials nullTrials

clear contrastMeanTrial contrastMeanAmp freqMeanTrial freqMeanAmp nullMeanTrial nullMeanAmp contrastFreqMeanTrial contrastFreqMeanAmp
clear contrastFreqRdm  
% figure(1)
% clf;


for iSub = 1:numSubs
    cd(subFolders{iSub});
    v=newView;
    % switch to the concatenation group
    v = viewSet(v, 'curGroup', 'Concatenation');
    nScans = viewGet(v, 'nscans');
    clear allTrials contrastTrials contrastAmp  freqTrials freqAmp  nullTrials nullAmp contrastFreqTrials contrastFreqAmp
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
%         sacRwd{rwdTypeNum}=[];
%         for r=1:length(s)
%             trialCorrectness{iSub,rwdTypeNum}(r,:) = s{r}.task{1}{1}.correctness;
%             trialResponse{iSub,rwdTypeNum,r} = s{r}.task{1}{1}.response;
%             propCorrect{iSub,rwdTypeNum}(r) = sum(trialCorrectness{iSub,rwdTypeNum}(r,:)) / length(trialCorrectness{iSub,rwdTypeNum}(r,:));
%             stimfile = s{r}.filename;
%             e{r} = myGetTaskEyeTraces(stimfile, 'removeBlink=3');
%             runEyeSize{rwdTypeNum}(r,:) = size(e{r}.eye.pupil);
%             sacRun{r} = [];
%             for itrial=1:e{r}.nTrials
%                 eyePos = [e{r}.eye.xPos(itrial,:)'  e{r}.eye.yPos(itrial,:)'];
%                 eyevel{r,itrial} = vecvel(eyePos,samplerate, 2);
%                 sac{r,itrial} = microsacc(eyePos, eyevel{r,itrial}, 1, 5);%saccade onset
%                 sacRun{r} = [sacRun{r}; sac{r,itrial}];
%             end
%             sacRwd{rwdTypeNum} = [sacRwd{rwdTypeNum}; sacRun{r}];
%         end
%         [sacRwd{rwdTypeNum}(:,1), ind] = sort(sacRwd{rwdTypeNum}(:,1));%saccade onset
%         sacRwd{rwdTypeNum}(:,2:7) = sacRwd{rwdTypeNum}(ind,2:7);%saccade end, velocity, horiz/vert component, horiz/vert amplitude
%         [binnedSac{rwdTypeNum},edges] = histcounts(sacRwd{rwdTypeNum}(:,1),max(sacRwd{rwdTypeNum}(:,1)));
%         
%         trialLengthEye(rwdTypeNum) = max(runEyeSize{rwdTypeNum}(:,2));
%         numTrials(rwdTypeNum) = sum(runEyeSize{rwdTypeNum}(:,1));
%         rwdPupil{rwdTypeNum} = NaN(numTrials(rwdTypeNum), trialLengthEye(rwdTypeNum));
%         trialCounter=0;
%         for r=1:length(s)
%             rwdPupil{rwdTypeNum}(trialCounter+1:trialCounter+runEyeSize{rwdTypeNum}(r,1), 1:runEyeSize{rwdTypeNum}(r,2)) = e{r}.eye.pupil;
%             trialCounter = trialCounter + runEyeSize{rwdTypeNum}(r,1);
%         end
%         %plot mean pupil size
%         meanPupil{rwdTypeNum} = nanmean(rwdPupil{rwdTypeNum})';
        
        
        
        %extract conditions for all scans within concatenation
        % doing this manually because the volNum is messed up, missed
        % triggers, instead of using getTaskParameters
        % names_: {'contrast'  'spatFreq'  'nullTrial'}
        contrastConcat{iScan} = [];
        freqConcat{iScan} = [];
        nullTrialConcat{iScan} = [];
        for r=1:length(s)%per run
            
            contrastConcat{iScan} = [contrastConcat{iScan} s{r}.task{2}{1}.randVars.contrast(2:trialsPerRun+1)];%skipping the first because of junked frames
            freqConcat{iScan} = [freqConcat{iScan} s{r}.task{2}{1}.randVars.spatFreq(2:trialsPerRun+1)];
            nullTrialConcat{iScan} = [nullTrialConcat{iScan} s{r}.task{2}{1}.randVars.nullTrial(2:trialsPerRun+1)];
        end
        
        contrastConcat{iScan} = contrastConcat{iScan}.*(1-nullTrialConcat{iScan});
        freqConcat{iScan} = freqConcat{iScan}.*(1-nullTrialConcat{iScan});
        
        totalTrials = length(contrastConcat{iScan});
        %contrast
        smallContrastMat = zeros(numContrasts,totalTrials);%number per trial
        for c=1:numContrasts
            smallContrastMat(c,contrastConcat{iScan}==c) = ones;
        end
        
        %spatial frequency
        smallFreqMat = zeros(numFreqs,totalTrials);%number per trial
        for c=1:numFreqs
            smallFreqMat(c,freqConcat{iScan}==c) = ones;
        end
        
        
        %spatial frequency & contrast!
        smallContrastSfMat = zeros(numFreqs*numContrasts,totalTrials);%number per trial
        for c=1:numFreqs*numContrasts
            curCond = (contrastConcat{iScan}-1)*numFreqs+freqConcat{iScan};
            smallContrastSfMat(c,curCond==c) = ones;
        end
        
        
        %null vs stimulated trials
        smallNullMat = zeros(2,totalTrials);%number per trial
        for c=0:1
            smallNullMat(c+1,nullTrialConcat{iScan}==c) = ones; %1=stim, 2=null
        end
        
        
        
        
%         %create  vector with a 1 for each trial, to be convolved with the
%         %task response (null trials)
%         taskSmallMat = ones(1,totalTrials);
%         taskMat = upsample(taskSmallMat,trialLength);%add zeros in between
        for iRoi=1:numRois
            roiTC{iRoi,rwd} = loadROITSeries(v, roiNames{iRoi}, iScan, [], 'keepNAN',true);
            concatInfo{iSub,rwd} = viewGet(v, 'concatInfo', iScan);
            allTrials{iRoi}(rwd,:,:,:) = reshape(roiTC{iRoi,rwd}.tSeries,[],trialLength,totalTrials);%voxels,trialLength, totalTrials
            %average trial
            trialTC(iSub,iRoi,rwd,:) = nanmean(reshape(nanmean(roiTC{iRoi,rwd}.tSeries), trialLength, totalTrials),2);
            
            %contrast
            for c=1:numContrasts
                temp = squeeze(allTrials{iRoi}(rwd,:,:,contrastConcat{iScan}==c));%voxels,time,trials
                contrastTrials{iRoi}(rwd,c,:,:) = squeeze(nanmean(temp,3));%average over trials
                contrastAmp{iRoi}(rwd,c,:) = std(contrastTrials{iRoi}(rwd,c,:,:),0,4);%3rd dim is voxels, 4th dim is time
                contrastMeanTrial(iSub,iRoi,rwd,c,:) = nanmean(contrastTrials{iRoi}(rwd,c,:,:),3);%average over voxels
                contrastMeanAmp(iSub,iRoi,rwd,c) = squeeze(nanmean(contrastAmp{iRoi}(rwd,c,:)));
            end
            
            %frequency
            for c=1:numFreqs
                temp = squeeze(allTrials{iRoi}(rwd,:,:,freqConcat{iScan}==c));
                freqTrials{iRoi}(rwd,c,:,:) = squeeze(nanmean(temp,3));%average over trials
                freqAmp{iRoi}(rwd,c,:) = std(freqTrials{iRoi}(rwd,c,:,:),0,4);%3rd dim is voxels, 4th dim is time
                freqMeanTrial(iSub,iRoi,rwd,c,:) = nanmean(freqTrials{iRoi}(rwd,c,:,:),3);%average over voxels
                freqMeanAmp(iSub,iRoi,rwd,c) = squeeze(nanmean(freqAmp{iRoi}(rwd,c,:)));
            end
            goodVox = ~isnan(sum(freqAmp{iRoi}(rwd,:,:)));
            corMat = corr(squeeze(freqAmp{iRoi}(rwd,:,goodVox))');
            freqRdm(iSub,iRoi,rwd,:) = squareform(1-corMat);
            
            
            for c=1:numFreqs*numContrasts
                temp = squeeze(allTrials{iRoi}(rwd,:,:,(contrastConcat{iScan}-1)*numFreqs+freqConcat{iScan}==c));
                contrastFreqTrials{iRoi}(rwd,c,:,:) = squeeze(nanmean(temp,3));%average over trials
                contrastFreqAmp{iRoi}(rwd,c,:) = std(contrastFreqTrials{iRoi}(rwd,c,:,:),0,4);%3rd dim is voxels, 4th dim is time
                contrastFreqMeanTrial(iSub,iRoi,rwd,c,:) = nanmean(contrastFreqTrials{iRoi}(rwd,c,:,:),3);%average over voxels
                contrastFreqMeanAmp(iSub,iRoi,rwd,c) = squeeze(nanmean(contrastFreqAmp{iRoi}(rwd,c,:)));
            end
            goodVox = ~isnan(sum(contrastFreqAmp{iRoi}(rwd,:,:)));
            corMat = corr(squeeze(contrastFreqAmp{iRoi}(rwd,:,goodVox))');
            contrastFreqRdm(iSub,iRoi,rwd,:) = squareform(1-corMat);
            
            %null/stim
            for c=0:1
                temp = squeeze(allTrials{iRoi}(rwd,:,:,nullTrialConcat{iScan}==c));
                nullTrials{iRoi}(rwd,c+1,:,:) = squeeze(nanmean(temp,3));%average over trials
                nullAmp{iRoi}(rwd,c+1,:) = std(nullTrials{iRoi}(rwd,c+1,:,:),0,4);%3rd dim is voxels, 4th dim is time
                nullMeanTrial(iSub,iRoi,rwd,c+1,:) = nanmean(nullTrials{iRoi}(rwd,c+1,:,:),3);%average over voxels
                nullMeanAmp(iSub,iRoi,rwd,c+1) = squeeze(nanmean(nullAmp{iRoi}(rwd,c+1,:)));
            end
            
        end
        
        
        
    end
    deleteView(v);
    %% return to home directory
    cd('..');
end

%% average trial for high/low reward
widelinewidth = 3;
narrowlinewidth = 0.5;
i=1;
figure(i); clf
rows=1;
cols=numRois;
rwdTypeStr = {'H','L'};
for iRoi=1:numRois
    subplot(rows,cols,iRoi)
    for rwd=1:2
        plot(squeeze(mean(trialTC(:,iRoi,rwd,:)))','linewidth',widelinewidth,'color',plotColors{rwd});
        hold all
        plot(squeeze(trialTC(:,iRoi,rwd,:))','linewidth',narrowlinewidth,'color',plotColors{rwd});
%         trialTC(iSub,iRoi,rwd,:)
    end
    xlabel('time(TRs)');
    ylabel('% signal change');
    title(roiNames{iRoi});
end

legend('high', 'low');
%% NULL/STIM figures

% %separate subplots for null/stim reward
% i=i+1; figure(i); clf
% rows=2;
% cols=numRois;
% nullStr = {'stim', 'null'};
% for iRoi=1:numRois
%        for rwd=1:2
%           subplot(rows,cols,iRoi + (rwd-1)*numRois)
%           plot(squeeze(mean(nullMeanTrial(:,iRoi,rwd,:,:)))','linewidth',widelinewidth,'color',plotColors{rwd});
%           hold all
%           for c=1:2
%             plot(squeeze(nullMeanTrial(:,iRoi,rwd,c,:))','linewidth',narrowlinewidth,'color',plotColors{rwd});
%           end
%           title([roiNames{iRoi} ' '  rwdTypeStr{rwd}]);
%        end
% end
% legend('stim, average', 'stim, individuals');


%separate subplots for high/low reward
i=i+1; figure(i); clf
rows=2;
cols=numRois;
nullStr = {'stim', 'null'};
for iRoi=1:numRois
       for rwd=1:2
          subplot(rows,cols,iRoi + (rwd-1)*numRois)
          
          for c=1:2
            plot(squeeze(nullMeanTrial(:,iRoi,rwd,c,:))','linewidth',narrowlinewidth,'color',plotColors{c});
            hold all
            plot(squeeze(mean(nullMeanTrial(:,iRoi,rwd,c,:)))','linewidth',widelinewidth,'color',plotColors{c});
          end
          title([roiNames{iRoi} ' '  rwdTypeStr{rwd}]);
       end
end
legend('stim', 'null');



i=i+1; figure(i); clf
rows=2;
for iRoi=1:numRois
    
    for rwd=1:2
        subplot(rows,cols,iRoi)
        plot(squeeze(mean(nullMeanAmp(:,iRoi,rwd,:))),'linewidth',widelinewidth);
        hold all
        if rwd==2; title(['null/stim mean amplitude,  ' roiNames{iRoi}]); end
        xticks([1 2]);
        xticklabels({'null','stim'});
        subplot(rows,cols,iRoi+cols)
        for c=1:2
            temp = mean(nullMeanTrial(:,iRoi,rwd,c,:));
            groupMeanTrialNull(iRoi,rwd,c,:) = temp;
            groupMeanAmpNull(iRoi,rwd,c) = std(temp);
        end
        temp = squeeze(groupMeanAmpNull(iRoi,rwd,:));
        plot(temp,'linewidth',widelinewidth);
        hold all
        if rwd==2; title(['amplitude of mean timecourse,  ' roiNames{iRoi}]); end
        xticks([1 2]);
        xticklabels({'null','stim'});
    end
end
legend(rwdTypeStr);

%% CONTRAST figures

%separate subplots for high/low reward
i=i+1; figure(i); clf
rows=2;
cols=numRois;
rwdTypeStr = {'H','L'};
contrastStr = {'low contrast', 'high contrast'};
for iRoi=1:numRois
       for rwd=1:2
          subplot(rows,cols,iRoi + (rwd-1)*numRois)
            
          for c=1:numContrasts
              plot(squeeze(contrastMeanTrial(:,iRoi,rwd,c,:))','linewidth',narrowlinewidth,'color',plotColors{c});
              hold all
              plot(squeeze(mean(contrastMeanTrial(:,iRoi,rwd,c,:)))','linewidth',widelinewidth,'color',plotColors{c});
          end
        title([roiNames{iRoi} ' '  rwdTypeStr{rwd}]);
       end
end
legend(contrastStr);

i=i+1; figure(i); clf
rows=2;
for iRoi=1:numRois
    
    for rwd=1:2
        subplot(rows,cols,iRoi)
        plot(squeeze(mean(contrastMeanAmp(:,iRoi,rwd,:))),'linewidth',widelinewidth);
        hold all
        if rwd==2; title(['contrast mean amplitude,  ' roiNames{iRoi}]); end
        subplot(rows,cols,iRoi+cols)
        for c=1:numContrasts
            temp = mean(contrastMeanTrial(:,iRoi,rwd,c,:));
            groupMeanTrialContrast(iRoi,rwd,c,:) = temp;
            groupMeanAmpContrast(iRoi,rwd,c) = std(temp);
        end
        temp = squeeze(groupMeanAmpContrast(iRoi,rwd,:));
        plot(temp,'linewidth',widelinewidth);
        hold all
        if rwd==2; title(['amplitude of mean timecourse,  ' roiNames{iRoi}]); end
    end
end
legend(rwdTypeStr);

%% FREQUENCY figures
i=i+1; figure(i); clf
freqStr = {'freq 1', 'freq 2', 'freq 3', 'freq 4', 'freq 5'};
rows=2;
cols=numRois;
for iRoi=1:numRois
       for rwd=1:2
          subplot(rows,cols,iRoi + (rwd-1)*numRois)
%           for c=1:numFreqs
%               plot(squeeze(freqMeanTrial(:,iRoi,rwd,c,:))','linewidth',narrowlinewidth,'color',plotColors{c});
%               hold all
%               plot(squeeze(mean(freqMeanTrial(:,iRoi,rwd,c,:)))','linewidth',widelinewidth,'color',plotColors{c});
%           end
            plot(squeeze(mean(freqMeanTrial(:,iRoi,rwd,:,:)))','linewidth',widelinewidth)
          hold all
          title([roiNames{iRoi} ' '  rwdTypeStr{rwd}]);
       end
end
legend(freqStr);

i=i+1; figure(i); clf
rows=2;
for iRoi=1:numRois
    for rwd=1:2
        subplot(rows,cols,iRoi)
        plot(squeeze(mean(freqMeanAmp(:,iRoi,rwd,:))),'linewidth',widelinewidth);
        hold all
        if rwd==2; title(['mean amplitude,  ' roiNames{iRoi}]); end
        subplot(rows,cols,iRoi + cols)
        
        for c=1:numFreqs
            temp = mean(freqMeanTrial(:,iRoi,rwd,c,:));
            groupMeanTrialFreq(iRoi,rwd,c,:) = temp;
            groupMeanAmpFreq(iRoi,rwd,c) = std(temp);
        end
        temp = squeeze(groupMeanAmpFreq(iRoi,rwd,:));
        plot(temp,'linewidth',widelinewidth);
        hold all
        if rwd==2; title(['amplitude of mean timecourse,  ' roiNames{iRoi}]); end
    end
end
legend(rwdTypeStr);

%% CONTRAST & FREQUENCY figures
i=i+1; figure(i); clf
freqStr = {'freq 1', 'freq 2', 'freq 3', 'freq 4', 'freq 5'};
rows=2;
cols=numRois;

rows=2;
for iRoi=1:numRois
    
    for rwd=1:2
        subplot(rows,cols,iRoi)
        temp = squeeze(mean(contrastFreqMeanAmp(:,iRoi,rwd,:)));
        plot([temp(1:5); NaN; temp(6:end)],'linewidth',widelinewidth);
        hold all
        if rwd==2; title(['mean amplitude,  ' roiNames{iRoi}]); end
        subplot(rows,cols,iRoi+cols)
        for c=1:numContrasts*numFreqs
            temp = mean(contrastFreqMeanTrial(:,iRoi,rwd,c,:));
            groupMeanTrialContrastFreq(iRoi,rwd,c,:) = temp;
            groupMeanAmpContrastFreq(iRoi,rwd,c) = std(temp);
        end
        temp = squeeze(groupMeanAmpContrastFreq(iRoi,rwd,:));
        plot([temp(1:5); NaN; temp(6:end)],'linewidth',widelinewidth);
        hold all
        if rwd==2; title(['amplitude of mean timecourse,  ' roiNames{iRoi}]); end
        
    end
end
legend(rwdTypeStr);

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
close all
clear all
tic
dataFolder = '/Volumes/MH02086153MACDT-Drobo/rwdFmri/';
concatProj=0;
subFolders = {'010620190823','001920190515','008020190509','007420190529',...
    '007520190715','008320190523','009020190515','009720190603','010420190805',...
    '010720190911',...
    '005520190411','008620190410','008720190503','009920190627'};
numSubs=length(subFolders);

junkedFrames = 10;
trialLength=10;
TR=1.5;
ecgselect=0.009;
respselect = 0.1;
display=0;
ecgSampleRate = 50;
ecgTrial = ecgSampleRate*TR*trialLength;
trialsPerRun=16;
ecgRunLength = ecgTrial*(trialsPerRun+1);
ecgInterpMethod = 'linear';
respInterpMethod = 'linear';
respStdWindow = 6*ecgSampleRate;
deconvLength = 10;
cd(dataFolder);
mrQuit;
if exist('v')
    deleteView(v);
end

clear d concatInfo sumResponse subResponse roiMeanTseries meanResponse roiTC subRunResponse trialCorrectness trialResponse trialRT propCorrect subMeanRunTC
regressBetasGlobal={};
plotColors = {[1 0 0], [0 0 1], [0 1 0], [0.5 1 0.2],[0 0 1], [1 0 0], [0 1 0], [0.5 1 0.2]};
plotStyles = {'-','--',':','-.','-','--',':','-.','-','--',':','-.','-','--',':','-.'};
rows=2;
cols=numSubs;
for iSub = 1:numSubs
    
    cd(subFolders{iSub});
    v=newView;
    
    % Find the correct Concatenation group
    if concatProj
        concatGroupNum = viewGet(v,'groupNum','concatProj'); %concatenation of multiple sessions
    else
        concatGroupNum = viewGet(v,'groupNum','concatSessions'); %concatenation of multiple sessions
        if isempty(concatGroupNum)%single session
            concatGroupNum = viewGet(v,'groupNum','Concatenation');
        end
    end
    
    % switch to the concatenation group
    v = viewSet(v, 'curGroup', concatGroupNum);
    nScans = viewGet(v, 'nscans');
    clear rois
    
    for iScan = 1:2%nScans%2 concatenations, 1 for each reward type
        s = viewGet(v, 'stimfile', iScan);
%         if ~isfield(s{1}, 'stimulus')
%             s=s{1};
%         end
        rwdType = s{1}.myscreen.stimulus.rewardType;
        rwdValue = s{1}.myscreen.stimulus.rewardValue;
        if strcmp(rwdType, 'H')
            rwd = 1;
        elseif strcmp(rwdType, 'L')
            rwd = 2;
        else
            disp('wtf');
            keyboard
        end
        numRuns(iSub,rwd) = length(s);
        concatInfo{iSub,rwd} = viewGet(v, 'concatInfo', iScan);

        %get physio data
        for r=1:numRuns(iSub,rwd)
            physioRunRwd{iSub,iScan}(r) = s{r}.myscreen.stimulus.rewardType;
            physioRunRwdVal{iSub,iScan}(r) = s{r}.myscreen.stimulus.rewardValue;
            
            %RESPIRATION
            respfilename = fullfile(s{r}.myscreen.originalmlrdir,s{r}.myscreen.physiodir,s{r}.myscreen.respirationfilename);
            if ~isfile(respfilename)
                %if the directory has been moved, but this is the original
                %session this should work.
                respfilename = fullfile(s{r}.myscreen.physiodir,s{r}.myscreen.respirationfilename);
                if ~isfile(respfilename)%concatenated session, and original directory has been moved?
                    keyboard
                end
            end
            resp{iSub,rwd,r}=load(respfilename);
            [respPeaks{iSub,rwd,r},criterion] = pickpeaks(resp{iSub,rwd,r},respselect,display);
            respPeaksDiff{iSub,rwd,r} = diff(respPeaks{iSub,rwd,r});
            respPeaksAmp{iSub,rwd,r} = resp{iSub,rwd,r}(respPeaks{iSub,rwd,r});
            interpPeaksAmp{iSub,rwd,r} = interp1(respPeaks{iSub,rwd,r}, respPeaksAmp{iSub,rwd,r}, 1:length(resp{iSub,rwd,r}), respInterpMethod, 'extrap');
            interpPeaksAmp{iSub,rwd,r}(end:ecgRunLength) = NaN;
            rv{iSub,rwd,r} = zeros(size(resp{iSub,rwd,r}));
            
            %             [respTroughs{iSub,rwd,r},criterion] = pickpeaks(-resp{iSub,rwd,r},respselect,display);%gives the times
            %             respTroughsAmp{iSub,rwd,r} = resp{iSub,rwd,r}(respTroughs{iSub,rwd,r});
            %             interpTroughsAmp{iSub,rwd,r} = interp1(respTroughs{iSub,rwd,r}, respTroughsAmp{iSub,rwd,r}, 1:length(resp{iSub,rwd,r}), respInterpMethod, 'extrap');
            %             interpTroughsAmp{iSub,rwd,r}(end:ecgRunLength) = NaN;
            %
            %             respRateTime{iSub,rwd,r} = respPeaks{iSub,rwd,r}(1:end-1) + 0.5*diff(respPeaks{iSub,rwd,r});%timepoints between the peaks
            %             interpRespPeaksDiff{iSub,rwd,r} = interp1(respRateTime{iSub,rwd,r}, respPeaksDiff{iSub,rwd,r}, 1:length(resp{iSub,rwd,r}), respInterpMethod, 'extrap');
            %             interpRespPeaksDiff{iSub,rwd,r}(end:ecgRunLength) = NaN;
            %
            %             rvt{iSub,rwd,r} = (interpPeaksAmp{iSub,rwd,r} - interpTroughsAmp{iSub,rwd,r})./interpRespPeaksDiff{iSub,rwd,r};
            %
            %we will junk the first 10 TRs of BOLD anyway, so starting from
            %6 sec. BUT, this measure is lagging.
            for t=respStdWindow:length(resp{iSub,rwd,r})
                rv{iSub,rwd,r}(t) = nanstd(resp{iSub,rwd,r}(t-respStdWindow+1:t));
            end
            rv{iSub,rwd,r}(1:respStdWindow) = rv{iSub,rwd,r}(t);
            rv{iSub,rwd,r}(end:ecgRunLength) = NaN;
%             rvt{iSub,rwd,r}(end:ecgRunLength) = NaN;
%             temp = reshape(rv{iSub,rwd,r}, ecgTrial,[]);
            %             rv{iSub,rwd,r}(1:respStdWindow) = squeeze(mean(temp(1:respStdWindow,:),2));%extrapolate to first trial from mean across trials
%             trialsRV{iSub,rwd}(r,:,:) = reshape(rv{iSub,rwd,r}, ecgTrial,[]);
            %             trialsRVT{iSub,rwd}(r,:,:) = reshape(rvt{iSub,rwd,r}, ecgTrial,[]);
%             meanTrialRvRun{iSub,rwd}(r,:) = squeeze(mean(trialsRV{iSub,rwd}(r,:,2:end),3));%mean across trials, exclude first trial
            %             meanTrialRvtRun{iSub,rwd}(r,:) = squeeze(mean(trialsRVT{iSub,rwd}(r,:,2:end),3));%mean across trials, exclude first trial
            
            %ECG
            %first load the ECG file
            ecgfilename = fullfile(s{r}.myscreen.originalmlrdir,s{r}.myscreen.physiodir,s{r}.myscreen.ecgfilename);
            if ~isfile(ecgfilename)
                %if the directory has been moved, but this is the original
                %session this should work.
                ecgfilename = fullfile(s{r}.myscreen.physiodir,s{r}.myscreen.ecgfilename);
                if ~isfile(ecgfilename)%concatenated session, and original directory has been moved?
                    keyboard
                end
            end
            ecg{iSub,rwd,r}=load(ecgfilename);
            [ecgPeaks{iSub,rwd,r},criterion] = pickpeaks(ecg{iSub,rwd,r},ecgselect,display);
            ecgPeaksDiff{iSub,rwd,r} = diff(ecgPeaks{iSub,rwd,r});%units are timepoints = 1/ecgSampleRate sec.
            ecgPeaksAmp{iSub,rwd,r} = ecg{iSub,rwd,r}(ecgPeaks{iSub,rwd,r});
            ecgRateTime{iSub,rwd,r} = ecgPeaks{iSub,rwd,r}(1:end-1) + 0.5*diff(ecgPeaks{iSub,rwd,r});%timepoints between the peaks
            scaleFactor = mean(ecgPeaksAmp{iSub,rwd,r})/mean(ecgPeaksDiff{iSub,rwd,r});%for visualization purposes
            interpPeaksDiff{iSub,rwd,r} = interp1(ecgRateTime{iSub,rwd,r}, ecgPeaksDiff{iSub,rwd,r}, 1:length(ecg{iSub,rwd,r}), ecgInterpMethod, 'extrap');
            interpPeaksDiff{iSub,rwd,r}(end:ecgRunLength) = NaN;
            
            trialsPeakDiff{iSub,rwd}(r,:,:) = reshape(interpPeaksDiff{iSub,rwd,r}, ecgTrial,[]);
            meanTrialPeaksDiffRun{iSub,rwd}(r,:) = squeeze(mean(trialsPeakDiff{iSub,rwd}(r,:,:),3));%mean across trials
            
            ecgPulseRate{iSub,rwd,r} = 1./ecgPeaksDiff{iSub,rwd,r};
            interpPulseRate{iSub,rwd,r} = interp1(ecgRateTime{iSub,rwd,r}, ecgPulseRate{iSub,rwd,r}, 1:length(ecg{iSub,rwd,r}), ecgInterpMethod, 'extrap');
            interpPulseRate{iSub,rwd,r}(end:ecgRunLength) = NaN;
            
%             trialsPulseRate{iSub,rwd}(r,:,:) = reshape(interpPulseRate{iSub,rwd,r}, ecgTrial,[]);
%             meanTrialPulseRateRun{iSub,rwd}(r,:) = squeeze(mean(trialsPulseRate{iSub,rwd}(r,:,:),3));%mean across trials
            
        end
        
        %RESPIRATION
%         rwdMeanRVT(iSub,rwd,:) = mean(meanTrialRvtRun{iSub,rwd});
%         rwdMeanRV(iSub,rwd,:) = mean(meanTrialRvRun{iSub,rwd});
        rwdRvTC{iSub,rwd} = [];
%         rwdRvtTC{iSub,rwd} = [];
        for r=1:numRuns(iSub,rwd)
            %removing junked frames & concatenate
            rwdRvTC{iSub,rwd} = horzcat(rwdRvTC{iSub,rwd}, rv{iSub,rwd,r}(ecgSampleRate*1.5*junkedFrames+1:end)');
%             rwdRvtTC{iSub,rwd} = horzcat(rwdRvtTC{iSub,rwd}, rvt{iSub,rwd,r}(ecgSampleRate*1.5*junkedFrames+1:end));
        end
        %downsample pulse using median
        trRV = reshape(rwdRvTC{iSub,rwd}, ecgSampleRate*1.5,[]);
%         rwdSubMeanPulse(iSub,rwd,:) = mean(trRV,2);
%         trRVT = reshape(rwdRvtTC{iSub,rwd}, ecgSampleRate*1.5,[]);
        downsampledRV{iSub,rwd} = nanmedian(trRV);
%         downsampledRVT{iSub,rwd} = nanmedian(trRVT);

        %create design matrix for deconvolution
        designMatRv{iSub,rwd}(1,:) = downsampledRV{iSub,rwd};
%         designMatRvt{iSub,rwd}(1,:) = downsampledRVT{iSub,rwd};
        for i=2:deconvLength
            designMatRv{iSub,rwd}(i,:) = circshift(designMatRv{iSub,rwd}(i-1,:),1);
%             designMatRvt{iSub,rwd}(i,:) = circshift(designMatRvt{iSub,rwd}(i-1,:),1);
        end
%         designMatResp{iSub,rwd} = [designMatRv{iSub,rwd}; designMatRvt{iSub,rwd}];
        designMatResp{iSub,rwd} = [designMatRv{iSub,rwd}];%ONLY RV
        
        %ECG
%         rwdMeanPeaksDiff(iSub,rwd,:) = mean(meanTrialPeaksDiffRun{iSub,rwd});
%         rwdMeanPulseRate(iSub,rwd,:) = mean(meanTrialPulseRateRun{iSub,rwd});
        %concatenate all pulse traces for this rwd
        rwdPulseTC{iSub,rwd} = [];
        for r=1:numRuns(iSub,rwd)
            %removing junked frames & concatenate
            rwdPulseTC{iSub,rwd} = horzcat(rwdPulseTC{iSub,rwd}, interpPulseRate{iSub,rwd,r}(ecgSampleRate*TR*junkedFrames+1:end));
        end
        %downsample pulse using median
        trPulse = reshape(rwdPulseTC{iSub,rwd}, ecgSampleRate*1.5,[]);
        downsampledPulse{iSub,rwd} = nanmedian(trPulse);
        
        %create design matrix for deconvolution
        designMatPulse{iSub,rwd}(1,:) = downsampledPulse{iSub,rwd};
        for i=2:deconvLength
            designMatPulse{iSub,rwd}(i,:) = circshift(designMatPulse{iSub,rwd}(i-1,:),1);
        end
        designMatRespPulse{iSub,rwd} = [designMatResp{iSub,rwd}; designMatPulse{iSub,rwd}];
        
%         %KEEP ONLY GOOD TRIALS
%         designMatPulse{iSub,rwd} = designMatPulse{iSub,rwd}(:,goodTRs);
%         designMatResp{iSub,rwd} = designMatResp{iSub,rwd}(:,goodTRs);
%         designMatRespPulse{iSub,rwd} = designMatRespPulse{iSub,rwd}(:,goodTRs);
%         
        rwdPulseTrials{iSub,rwd} = reshape(rwdPulseTC{iSub,rwd},ecgTrial,[]);
        rwdRvTrials{iSub,rwd} = reshape(rwdRvTC{iSub,rwd},ecgTrial,[]);

        meanPulse(iSub,rwd,:) = mean(rwdPulseTrials{iSub,rwd},2);
        meanRV(iSub,rwd,:) = mean(rwdRvTrials{iSub,rwd},2);
    end
    
    subplot(rows,cols,iSub)
    for rwd=1:2
       plot( mean(rwdPulseTrials{iSub,rwd},2),'color',plotColors{rwd});
       hold on
    end
    subplot(rows,cols,iSub+cols)  
    for rwd=1:2
       plot( mean(rwdRvTrials{iSub,rwd},2),'color',plotColors{rwd});
       hold on
    end
    title(getLastDir(pwd));
    deleteView(v);
    cd ..
end
set(gcf,'position',[100 100 1000 500]);
%%
figure
subplot(1,2,1)
for rwd=1:2
   plot(squeeze(mean(meanPulse(:,rwd,:))),'color',plotColors{rwd});
   hold on
end
subplot(1,2,2)
for rwd=1:2
   plot(squeeze(mean(meanRV(:,rwd,:))),'color',plotColors{rwd});
   hold on
end
%%

ConcatProjStr = '';
if concatProj
    ConcatProjStr = 'concatProj';
end
%%
save([dataFolder 'rwdRapid_physio.mat'], 'concatInfo',  ...
    'subFolders','trialLength',...
    'ecgselect','ecgSampleRate','ecgTrial','ecgRunLength','ecgInterpMethod',...
    'ecg','ecgPulseRate','interpPulseRate',...
    'respselect','resp',...
    'rwdPulseTrials','rwdRvTrials','meanPulse','meanRV',...
    'designMatPulse','designMatResp','designMatRespPulse');

physioRunRwd
toc

% %% plot polar histogram of task related phase
% figure
% rows=2;
% cols=numSubs;
% nbins=20;
% nVox=1000;
% thresh=0.015;
% for iSub=1:numSubs
%     %     bestVox = sort(allVoxTaskCo{iSub,1}+allVoxTaskCo{iSub,2},'descend');
%     for rwd=1:2
%         subplot(rows,cols,iSub+(rwd-1)*cols)
%         %         polarhistogram(allVoxTaskPhase{iSub,rwd}(bestVox(1:nVox)),nbins);
%         polarhistogram(allVoxTaskPhase{iSub,rwd}(allVoxTaskCo{iSub,rwd}>thresh),nbins);
%         set(gca,'RTickLabel',[]);
%         set(gca,'ThetaTickLabel',[]);
%     end
% end
%%
function newConcat = zscoreConcat(concatData, concatInfo)%voxels,TRs
newConcat = [];
for iRun=1:size(concatInfo.runTransition,1)
    thisRun = concatData(:,concatInfo.runTransition(iRun,1):concatInfo.runTransition(iRun,2));
    newConcat = [newConcat zscore(thisRun,0,2)];
end
% newConcat = concatData;
end
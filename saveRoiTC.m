close all
mrQuit
clear all
tic

% dataFolder = '/Volumes/MH02086153MACDT-Drobo/rwdFmri/';
dataFolder = '/Volumes/TaskDrive/rwdFmri/';
saveFolder = '~/rwdFmri/';
subFolders = {'010620190823','001920190515','008020190509','007420190529',...
    '007520190715','008320190523','009020190515','009720190603','010420190805',...
    '010720190911',...
    '005520190411','008620190410','008720190503','009920190627'};

% subFolders = {'001920190515','005520190411','007420190529','007520190715','008020190509',...
%     '008320190523', '008620190410','008720190503','009020190515',...
%     '009720190603','009920190627','010420190805', '010620190823',...
%     '010720190911'};
% subFolders = {'001920190515','008020190509','007420190529',...
%     '007520190715','008320190523','009020190515','009720190603','010420190805',...
%     '010720190911',...
%     '005520190411','008620190410','008720190503','009920190627'};
% roiNames = {'benson','leftBenson', 'rightBenson'};
% roiNames = {'V1','V2', 'V3','Benson','lh_17Networks_16','rh_17Networks_16','leftCerebellarCortex','rightCerebellarCortex'};
roiNames = {'Benson','V2', 'V3','V1'};

eccMin = 0.2;
eccMax = 70;
nbins = 12;
binBorders = logspace(log10(eccMin),log10(eccMax),nbins+1);

for i=2:length(binBorders)
    binCenters(i-1) = (binBorders(i)+binBorders(i-1))/2;
end
curFolder = pwd;
cd(dataFolder);

trialsPerRun=16;
trialLength = 10;
nVolumes = trialsPerRun*trialLength;
numContrasts = 2;
numFreqs = 5;
frames=170;
junkedFrames = 10;
TR=1.5;
numSubs = length(subFolders);
rwdString = {'H','L'};
plotColors = {[1 0 0], [0 0 1]};
linewidth = 1;
markersize=10;

for toZscore=0:1
    for concatProj= 0:1
        
        for iSub=1:numSubs
            cd(subFolders{iSub});
            v=newView;
            if concatProj
                concatGroupNum = viewGet(v,'groupNum','concatProj'); %concatenation of multiple sessions
            else
                concatGroupNum = viewGet(v,'groupNum','concatSessions'); %concatenation of multiple sessions
                if isempty(concatGroupNum)%single session
                    concatGroupNum = viewGet(v,'groupNum','Concatenation');
                end
            end
            
            %get Benson maps
            v = viewSet(v, 'curGroup', 'templates');
            templateGroup = viewGet(v,'curGroup');
            v = loadAnalysis(v, 'mrDispOverlayAnal/templateRet.mat');
            for iRoi = 1:length(roiNames)
                bensonData = loadROIbensonMatching(v,roiNames{iRoi},1,templateGroup,1,concatGroupNum);
                if size(bensonData,1)>0%no benson data for other ROIs
                    eccen{iSub,iRoi} = bensonData{1}.eccen;
                    ang{iSub,iRoi} = bensonData{1}.ang;
                    areas{iSub,iRoi} = bensonData{1}.areas;
                end
            end
            
            
            % switch to the concatenation group
            v = viewSet(v, 'curGroup', concatGroupNum);
            
            nScans = viewGet(v, 'nscans');
            clear params
            
            for iscan = 1:2%nScans%2 concatenations, 1 for each reward type
                s = viewGet(v, 'stimfile', iscan);
                rwdType = s{1}.myscreen.stimulus.rewardType;
                if strcmp(rwdType, 'H')
                    rwd = 1;
                elseif strcmp(rwdType, 'L')
                    rwd = 2;
                else
                    keyboard
                end
                numRuns(iSub,rwd) = length(s);
                concatInfo{iSub,rwd} = viewGet(v, 'concatInfo', iscan);
                numTRs(iSub,rwd) = length(concatInfo{iSub,rwd}.whichScan);
                
                %         v = viewSet(v, 'curGroup', 'Concatenation');
                %initialize concatenated trial events
                nullTrials{iSub,rwd}=[]; contrastTrials{iSub,rwd}=[]; freqTrials{iSub,rwd}=[];
                for r=1:numRuns(iSub,rwd)
                    nullTrialsRun{iSub,rwd}(r,:) = s{r}.task{2}{1}.randVars.nullTrial(2:trialsPerRun+1);%skipping the first because of junked frames.0=stim, 1=null
                    nullTrials{iSub,rwd} = [nullTrials{iSub,rwd} nullTrialsRun{iSub,rwd}(r,:)];
                    
                    contrastTrialsRun{iSub,rwd}(r,:) = s{r}.task{2}{1}.randVars.contrast(2:trialsPerRun+1).*(1-nullTrialsRun{iSub,rwd}(r,:));
                    freqTrialsRun{iSub,rwd}(r,:) = s{r}.task{2}{1}.randVars.spatFreq(2:trialsPerRun+1).*(1-nullTrialsRun{iSub,rwd}(r,:));
                    
                    contrastTrials{iSub,rwd} = [contrastTrials{iSub,rwd} contrastTrialsRun{iSub,rwd}(r,:)];
                    freqTrials{iSub,rwd} = [freqTrials{iSub,rwd} freqTrialsRun{iSub,rwd}(r,:)];
                    
                end
                temp = repmat(nullTrials{iSub,rwd},trialLength,1);
                nullTrialsTRs{iSub,rwd} = temp(:);
                
                for iRoi = 1:length(roiNames)
                    roiTC{iSub,iRoi,rwd} = loadROITSeries(v, roiNames{iRoi}, iscan, [], 'keepNAN',true);
                    roiTC{iSub,iRoi,rwd}.tSeries = zscoreConcat(roiTC{iSub,iRoi,rwd}.tSeries, concatInfo{iSub,rwd});
                    numVox(iSub,iRoi) = size(roiTC{iSub,iRoi,rwd}.tSeries,1);
                    nullTseries{iSub,iRoi,rwd} = roiTC{iSub,iRoi,rwd}.tSeries(:,nullTrialsTRs{iSub,rwd}==1);
                    nullTrialTseries{iSub,iRoi,rwd} = reshape(nullTseries{iSub,iRoi,rwd},numVox(iSub,iRoi),trialLength,[]);%iSub,T,trial
                    stimTseries{iSub,iRoi,rwd} = roiTC{iSub,iRoi,rwd}.tSeries(:,nullTrialsTRs{iSub,rwd}==0);
                    stimTrialTseries{iSub,iRoi,rwd} = reshape(stimTseries{iSub,iRoi,rwd},numVox(iSub,iRoi),trialLength,[]);%iSub,T,trial
                end
            end
            
            deleteView(v);
            
            
            %% return to home directory
            cd('..');
        end
        %%
        zScoreString = '';
        if toZscore
            zScoreString = '_zscored';
        end
        concatProjStr = '';
        if concatProj
            concatProjStr = 'proj';
        end
        save([saveFolder 'roiTC_' zScoreString concatProjStr '.mat'], 'dataFolder', 'subFolders', 'roiNames', ...
            'numRuns','numTRs','concatInfo',...
            'frames', 'junkedFrames', 'TR', 'trialsPerRun', 'trialLength', 'nVolumes',...
            'nbins','binBorders','binCenters',...
            'eccen','ang','areas',...
            'numVox','roiTC','nullTrialsTRs',...
            'nullTrialsRun','nullTrials','contrastTrialsRun','freqTrialsRun','contrastTrials','freqTrials',...
            'nullTseries','nullTrialTseries','stimTseries','stimTrialTseries','-v7.3');
        toc
    end
end
%% mean null trials

figure(1); clf
rows=length(roiNames);
cols = numSubs;
for iSub=1:numSubs
    for iRoi=1:length(roiNames)
        subplot(rows,cols,iSub + (iRoi-1)*cols)
        for rwd=1:2
            temp = squeeze(nanmean(nullTrialTseries{iSub,iRoi,rwd}));
            plot(mean(temp,2),'color',plotColors{rwd});
            hold on
        end
    end
end


%%
function newConcat = zscoreConcat(concatData, concatInfo)%voxels,TRs
newConcat = [];
for iRun=1:size(concatInfo.runTransition,1)
    thisRun = concatData(:,concatInfo.runTransition(iRun,1):concatInfo.runTransition(iRun,2));
    
    newConcat = [newConcat zscore(thisRun,0,2)];
end

end




%extract all null trials, concatenate and perform correlation analysis
dataFolder = '~/data/rwdFmri/';
dataFolder = '/Volumes/MH02086153MACDT-Drobo/rwdFmri/';
% subFolders = {'009020190515', '009020190325/', '009720190603', '009720190531', '007520190715','007520190318/', '001920190515','001920190403/', '007420190529', '007420190221/', '008320190523', '008320190320/', '008720190503', '008620190410', '005520190411', '008020190509'};
% subFolders = {'010720190911','010720190829'};
subFolders = {'010720190911'};
subFolders = {'010620190823'};
subFolders = {'010620190823','001920190515','008020190509','007420190529',...
    '007520190715','008320190523','009020190515','009720190603','010420190805',...
    '010720190911',...
    '005520190411','008620190410','008720190503','009920190627'};
subFolders = {'008020190509','007420190529',...
    '007520190715','008320190523','009020190515','009720190603','010420190805',...
    '005520190411','008620190410','008720190503','009920190627'};


curFolder = pwd;
cd(dataFolder);
mrQuit;
if exist('v')
    deleteView(v);
end
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
for iSub=1:numSubs
    cd(subFolders{iSub});
    v=newView;
    
    
    for concatProj=0:1 %create concatenations for null and stim trials, both without and with global signal regress
        concatProjString = '';
        if concatProj
            concatProjString = '_proj';
        end
        
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
        clear params
        
        for iscan = 1:2%nScans%2 concatenations, 1 for each reward type
            s = viewGet(v, 'stimfile', iscan);
            if isstruct(s{1})
                rwdType = s{1}.myscreen.stimulus.rewardType;
            else
                rwdType = s{1}{1}.myscreen.stimulus.rewardType;
            end
            
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
            nullTrials{iSub,rwd}=[];
            for r=1:numRuns(iSub,rwd)
                nullTrialsRun = s{r}.task{2}{1}.randVars.nullTrial(2:trialsPerRun+1);%skipping the first because of junked frames.0=stim, 1=null
                nullTrials{iSub,rwd} = [nullTrials{iSub,rwd} nullTrialsRun];
            end
            temp = repmat(nullTrials{iSub,rwd},trialLength,1);
            nullTrialsTRs = temp(:);
            %load all voxels timecourse
            origTSeries = loadTSeries(v,iscan);
            origSize{iSub,rwd} = size(origTSeries);
            allTseries{iSub,rwd} = reshape(origTSeries,[],numTRs(iSub,rwd));%voxels,TRs
            nullTseries{iSub,rwd} = allTseries{iSub,rwd}(:,nullTrialsTRs==1);
            stimTseries{iSub,rwd} = allTseries{iSub,rwd}(:,nullTrialsTRs==0);
            
        end
        allNullTseries{iSub} = [nullTseries{iSub,1} nullTseries{iSub,2}];
        allStimTseries{iSub} = [stimTseries{iSub,1} stimTseries{iSub,2}];
        
        if concatProj==0
            v=viewSet(v,'deletegroup','nullConcat');
            v=viewSet(v,'newGroup','nullConcat');
        end
        v=viewSet(v,'curGroup','nullConcat');
        hdrScan = 1;
        hdrGroup = 'Raw';
        fullhdr = viewGet(v,'niftihdr',hdrScan,hdrGroup);
        hdr=cbiCreateNiftiHeader('qform44',fullhdr.qform44, 'sform44',fullhdr.sform44,'sform_code',fullhdr.sform_code,'vox_offset',fullhdr.vox_offset);

        
        nullTC = allNullTseries{iSub};
        nullTC = reshape(nullTC,origSize{iSub,rwd}(1),origSize{iSub,rwd}(2),origSize{iSub,rwd}(3),[]);
        
        scanNum = concatProj*2+1;
        scanCycles(scanNum) = size(nullTC,4)/trialLength;
        scanParams.fileName = ['null' concatProjString '.nii'];
        v=saveNewTSeries(v,nullTC,scanParams,hdr);
        
        stimTC = allStimTseries{iSub};
        stimTC = reshape(stimTC,origSize{iSub,rwd}(1),origSize{iSub,rwd}(2),origSize{iSub,rwd}(3),[]);
        scanNum = concatProj*2+2;
        scanCycles(scanNum) = size(stimTC,4)/trialLength;
        scanParams.fileName = ['stim' concatProjString '.nii'];
        v=saveNewTSeries(v,stimTC,scanParams,hdr);
    end
    [v, corParams] = corAnal(v,[],'justGetParams=1','defaultParams=1');
    corParams.ncycles = [scanCycles];%2X scans,  both before and fter global mean projection
    corParams.recompute = ones(size(corParams.recompute));
    v = corAnal(v,corParams);
    
    deleteView(v);
    
    
    %% return to home directory
    cd('..');
end

nullTrialsPredictor=0;%either only null trials (1) or all trials (0)
regressType = 0;%1=pca, 2=global mean
regressScaling = 1;

dataFolder = '~/data/rwdFmri/';
samplerate=500;
subFolders = {'009720190603', '009720190531', '008020190509','008020190411', '001920190515','001920190403/', '007420190529', '007420190221/', '008320190523', '008320190320/', '008720190503', '008620190410', '005520190411', '007520190318/', '009020190325/'};
subFolders = {'009020190515', '009020190325/', '009720190603', '009720190531', '008020190509','008020190411', '001920190515','001920190403/', '007420190529', '007420190221/', '008320190523', '008320190320/', '008720190503', '008620190410', '005520190411', '007520190318/'};
subFolders = {'009020190515', '009020190325/', '009720190603', '009720190531', '001920190515','001920190403/', '007420190529', '007420190221/', '008320190523', '008320190320/', '008720190503', '008620190410', '005520190411', '007520190318/', '008020190509'};
% subFolders = {'009020190325/', '009720190531',  '008320190320/', '008620190410', '005520190411', '007520190318/', '008020190509'};

% subFolders = {'009020190325', '009720190531'};
% subFolders = {'005520190411', '007520190318/', '008020190509','008020190411'};

roiNames = {'leftBenson','rightBenson'};
roiNames = {'benson'};

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
trialsPerRun=16;
trialLength = 10;
nVolumes = trialsPerRun*trialLength;
numContrasts = 2;
numFreqs = 5;

plotColors = {[0 0 1], [1 0 0], [0 1 0], [0.5 1 0.2]};
plotStyles = {'-','--',':','-.','-','--',':','-.'};

clear binMeanTseries binTrialTC numBinVoxels binHrf
clear contrastBetas sfBetas sfContrastBetas nullBetas nullDeconv contrastDeconv freqDeconv contrastFreqDeconv

deconvLength = 10;

%bins by log eccentricity
% eccMin = 0.5;
% eccMax = 20;
% nbins = 5;
% binBorders = logspace(log10(eccMin),log10(eccMax),nbins+1);
% binBorders = [0 1.3 3 7 20 100];
% binBorders = [0 1 3 12 100];
% binBorders = [0 1.3 3 7 20 100];
% binBorders = [0 2.5 6 15 100];
binBorders = [0 0.4 1.3 2 3 5 7 12 20 40 100];
% binBorders = [0 1 3 7 20 100];
nbins = length(binBorders)-1;

for iSub = 1:numSubs
    cd(subFolders{iSub});
    v=newView;
    
    %load benson eccentricity maps
    v = viewSet(v, 'curGroup', 'Concatenation');
    concatGroupNum = viewGet(v,'curGroup');
    v = viewSet(v, 'curGroup', 'templates');
    templateGroup = viewGet(v,'curGroup');
    v = loadAnalysis(v, 'mrDispOverlayAnal/templateRet.mat');
    for iRoi = 1:length(roiNames)
        bensonData = loadROIbensonMatching(v,roiNames{iRoi},1,templateGroup,1,concatGroupNum);
        eccen{iSub,iRoi} = bensonData{1}.eccen;
        ang{iSub,iRoi} = bensonData{1}.ang;
        areas{iSub,iRoi} = bensonData{1}.areas;
    end
    
    % switch to the concatenation group
    v = viewSet(v, 'curGroup', 'Concatenation');
    nScans = viewGet(v, 'nscans');
    
    for iscan = 1:nScans%2 concatenations, 1 for each reward type
%         d{iscan} = viewGet(v, 'd', iscan);
        
        %             concatInfo{iSub,iScan} = viewGet(v, 'concatInfo', iScan);
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
        
        %load all voxels timecourse, to get global PCA, and global task response
        temp = reshape(loadTSeries(v,iscan),[],numTRs(iSub,rwd));%voxels,TRs
        allTseries{rwd} = zscoreConcat(temp, concatInfo{iSub,rwd});
        
%         p=pca(allTseries{rwd});
%         taskPca{iSub,rwd} = zscore(p(:,1));
        globalMean{iSub,rwd} = nanmean(allTseries{rwd})';
        
        for iRoi=1:numRois
            roiTC{iSub,iRoi,rwd} = loadROITSeries(v, roiNames{iRoi}, iscan, [], 'keepNAN',true);    
            roiTC{iSub,iRoi,rwd}.tSeries = zscoreConcat(roiTC{iSub,iRoi,rwd}.tSeries, concatInfo{iSub,rwd});
            
            nvox = roiTC{iSub,iRoi,rwd}.n;
            for ibin=1:nbins
                binVoxels = eccen{iSub,iRoi}>binBorders(ibin) & eccen{iSub,iRoi}<=binBorders(ibin+1);
                numBinVoxels(iSub,iRoi,ibin) = sum(binVoxels);
                %                 binVoxels = binVoxels & areas{iSub,iRoi}==1;%ONLY V1
                binVoxels = binVoxels & areas{iSub,iRoi}>0;%ONLY V1,V2,V3
                binMeanTseries{iSub,iRoi,ibin,rwd} = nanmean(roiTC{iSub,iRoi,rwd}.tSeries(binVoxels,:));%mean timecourse across voxels
                
                %perform regression of PCA on the mean bin timecourse
                regressBetasPca(iSub,rwd,iRoi,ibin) = taskPca{iSub,rwd}\binMeanTseries{iSub,iRoi,ibin,rwd}';
                pcaCorr(iSub,rwd,iRoi,ibin) = corr(taskPca{iSub,rwd},binMeanTseries{iSub,iRoi,ibin,rwd}');

                %perform regression of Global Signal on the mean bin timecourse
                regressBetasGlobal(iSub,rwd,iRoi,ibin) = globalMean{iSub,rwd}\binMeanTseries{iSub,iRoi,ibin,rwd}';
                globalCorr(iSub,rwd,iRoi,ibin) = corr(globalMean{iSub,rwd},binMeanTseries{iSub,iRoi,ibin,rwd}');

                if regressType==1
                    binMeanTseries{iSub,iRoi,ibin,rwd} = binMeanTseries{iSub,iRoi,ibin,rwd} - regressScaling*regressBetasPca(iSub,rwd,iRoi,ibin)'*taskPca{iSub,rwd}';
                elseif regressType==2
                    binMeanTseries{iSub,iRoi,ibin,rwd} = binMeanTseries{iSub,iRoi,ibin,rwd} - regressScaling*regressBetasGlobal(iSub,rwd,iRoi,ibin)'*globalMean{iSub,rwd}';
                end

                temp = reshape(binMeanTseries{iSub,iRoi,ibin,rwd},trialLength,[]);
                binTrialTC(iSub,iRoi,ibin,rwd,:) = mean(temp,2);
            end
        end
        
        
        %extract conditions for all scans within concatenation
        deconvMatNull{rwd} = [];
        deconvMatContrast{rwd} = [];
        deconvMatFreq{rwd} = [];
        deconvMatContrastFreq{rwd} = [];
        
        %create  vector with a 1 for each trial in a single run, to be convolved with the
        %task response (null trials). This is identical for all runs!
        taskSmallMat = ones(1,trialsPerRun);
        taskMat = upsample(taskSmallMat,trialLength);%add zeros in between
%         clear contrastMat sfMat contrastSfMat nullTrialMat
        clear designMatContrastFreqBinRun designMatFreqBinRun designMatContrastBinRun designMatNullBinRun 
        for r=1:numRuns(iSub,rwd)%create design matrices per run
            nullTrials = s{r}.task{2}{1}.randVars.nullTrial(2:trialsPerRun+1);%skipping the first because of junked frames.0=stim, 1=null
            runContrastEvents = s{r}.task{2}{1}.randVars.contrast(2:trialsPerRun+1);
            runFreqEvents = s{r}.task{2}{1}.randVars.spatFreq(2:trialsPerRun+1);
            runContrastEvents = runContrastEvents.*(1-nullTrials);
            runFreqEvents = runFreqEvents.*(1-nullTrials);
            
            
            %change the task predictor to correspond only to null trials
            if nullTrialsPredictor
                taskSmallMat(nullTrials==0) = 0;
            end
            
            %null vs stimulated trials
            smallNullMat = zeros(2,trialsPerRun);%number per trial
            for c=0%:1
                smallNullMat(c+1,nullTrials==c) = ones; %0=stim, 1=null
            end
            smallNullMat(end,:) = taskSmallMat;%all trial onsets
            nullTrialMat{r,rwd} = upsample(smallNullMat,trialLength);%number per volume

            %contrast
            smallContrastMat = zeros(numContrasts,trialsPerRun);
            for c=1:numContrasts
                smallContrastMat(c,runContrastEvents==c) = ones;%1=low, 2=high
            end
            smallContrastMat(end+1,:) = taskSmallMat;%all trial onsets
            contrastMat{r,rwd} = upsample(smallContrastMat,trialLength);%number per volume
            
            %spatial frequency
            smallFreqMat = zeros(numFreqs,trialsPerRun);%number per trial
            for c=1:numFreqs
                smallFreqMat(c,runFreqEvents==c) = ones;
            end
            smallFreqMat(end+1,:) = taskSmallMat;%all trial onsets
            freqMat{r,rwd} = upsample(smallFreqMat,trialLength);%number per volume
            
            
            %spatial frequency & contrast!
            smallContrastSfMat = zeros(numFreqs*numContrasts,trialsPerRun);%number per trial
            for c=1:numFreqs*numContrasts
                curCond = (runContrastEvents-1)*numFreqs+runFreqEvents;
                smallContrastSfMat(c,curCond==c) = ones;
            end
            smallContrastSfMat(end+1,:) = taskSmallMat;%all trial onsets
            contrastSfMat{r,rwd} = upsample(smallContrastSfMat,trialLength);%number per volume
            
            
            
            
            %nullTrial deconvolution matrix
            deconvMatNullRun = zeros((size(nullTrialMat{r,rwd},1)-1)*deconvLength,size(nullTrialMat{r,rwd},2));
            for c=1:size(nullTrialMat{r,rwd},1)
                for j=1:deconvLength
                    deconvMatNullRun((c-1)*deconvLength+j,:) = [zeros(1,j-1) nullTrialMat{r,rwd}(c,1:end-j+1)];%shift and pad with zeros
                end
            end
%             deconvMatNullRun(end+1,:) = ones;
            %high-pass filter
            deconvMatNullRun = highpassConcat(deconvMatNullRun, concatInfo{iSub,rwd}, r);
            deconvMatNullRun(end+1,:) = ones;
            deconvMatNull{rwd} = [deconvMatNull{rwd} deconvMatNullRun];
            
            %contrast deconvolution matrix
            deconvMatContrastRun = zeros((size(contrastMat{r,rwd},1)-1)*deconvLength,size(contrastMat{r,rwd},2));
            for c=1:size(contrastMat{r,rwd},1)
                for j=1:deconvLength
                    deconvMatContrastRun((c-1)*deconvLength+j,:) = [zeros(1,j-1) contrastMat{r,rwd}(c,1:end-j+1)];%shift and pad with zeros
                end
            end
%             deconvMatContrastRun(end+1,:) = ones;
            %high-pass filter
            deconvMatContrastRun = highpassConcat(deconvMatContrastRun, concatInfo{iSub,rwd}, r);
            deconvMatContrastRun(end+1,:) = ones;
            deconvMatContrast{rwd} = [deconvMatContrast{rwd} deconvMatContrastRun];
            
            %sf deconvolution matrix
            deconvMatFreqRun = zeros((size(freqMat{r,rwd},1)-1)*deconvLength,size(freqMat{r,rwd},2));
            for c=1:size(freqMat{r,rwd},1)
                for j=1:deconvLength
                    deconvMatFreqRun((c-1)*deconvLength+j,:) = [zeros(1,j-1) freqMat{r,rwd}(c,1:end-j+1)];%shift and pad with zeros
                end
            end
%             deconvMatFreqRun(end+1,:) = ones;
            %high-pass filter
            deconvMatFreqRun = highpassConcat(deconvMatFreqRun, concatInfo{iSub,rwd}, r);
            deconvMatFreqRun(end+1,:) = ones;
            deconvMatFreq{rwd} = [deconvMatFreq{rwd} deconvMatFreqRun];
            
            %sf & contrast deconvolution matrix
            deconvMatContrastFreqRun = zeros((size(contrastSfMat{r,rwd},1)-1)*deconvLength,size(contrastSfMat{r,rwd},2));
            for c=1:size(contrastSfMat{r,rwd},1)
                for j=1:deconvLength
                    deconvMatContrastFreqRun((c-1)*deconvLength+j,:) = [zeros(1,j-1) contrastSfMat{r,rwd}(c,1:end-j+1)];%shift and pad with zeros
                end
            end
%             deconvMatContrastFreqRun(end+1,:) = ones;
            %high-pass filter
            deconvMatContrastFreqRun = highpassConcat(deconvMatContrastFreqRun, concatInfo{iSub,rwd}, r);
            deconvMatContrastFreqRun(end+1,:) = ones;
            deconvMatContrastFreq{rwd} = [deconvMatContrastFreq{rwd} deconvMatContrastFreqRun];
            
            
        end
    end
    
    
    %get mean response per ROI/bin - to be used instead of canonical HRF
    for iRoi=1:numRois
        for ibin=1:nbins
            binHrf(iSub,iRoi,ibin,:) = (numTRs(iSub,1)*binTrialTC(iSub,iRoi,ibin,1,:) + numTRs(iSub,2)*binTrialTC(iSub,iRoi,ibin,2,:))./(sum(numTRs(iSub,:)));
            %         voxMeanResp{iSub,iRoi} = (voxMeanRespRwd{iSub,iRoi,1}+voxMeanRespRwd{iSub,iRoi,2})/2;
        end
    end
    
    for rwd=1:2
        for iRoi=1:numRois
            %create design matrices for each eccentricity bin, using mean response
            %in place of canonical hrf
            
            for ibin=1:nbins
                tc = squeeze(binMeanTseries{iSub,iRoi,ibin,rwd})';
                hrf = squeeze(zscore(binHrf(iSub,iRoi,ibin,:)));
                
                designMatContrastBin=[];
                designMatFreqBin=[];
                designMatContrastFreqBin=[];
                designMatNullBin=[];
                
                
                for r=1:numRuns(iSub,rwd)%create design matrices per run. at the beginning of each run the predictors must be zeroed
                    
                    %contrast
                    clear designMatContrastBinRun
                    for c=1:size(contrastMat{r,rwd},1)
                        temp = conv(contrastMat{r,rwd}(c,:),hrf);
                        designMatContrastBinRun(c,:) =  temp(1:trialsPerRun*trialLength);
                    end
%                     designMatContrastBinRun(numContrasts+1,:) = ones;
                    %high-pass filter
                    designMatContrastBinRun = highpassConcat(designMatContrastBinRun, concatInfo{iSub,rwd}, r);
                    designMatContrastBinRun(end+1,:) = ones;
                    designMatContrastBin = [designMatContrastBin designMatContrastBinRun];
                    
                    
                    %sf
                    clear designMatFreqBinRun
                    for c=1:size(freqMat{r,rwd},1)
                        temp = conv(freqMat{r,rwd}(c,:),hrf);
                        designMatFreqBinRun(c,:) =  temp(1:trialsPerRun*trialLength);
                    end
%                     designMatFreqBinRun(numFreqs+1,:) = ones;
                    %high-pass filter
                    designMatFreqBinRun = highpassConcat(designMatFreqBinRun, concatInfo{iSub,rwd}, r);
                    designMatFreqBinRun(end+1,:) = ones;
                    designMatFreqBin = [designMatFreqBin designMatFreqBinRun];
                    
                    
                    %contrast & sf
                    clear designMatContrastFreqBinRun
                    for c=1:size(contrastSfMat{r,rwd},1)
                        temp = conv(contrastSfMat{r,rwd}(c,:),hrf);
                        designMatContrastFreqBinRun(c,:) =  temp(1:trialsPerRun*trialLength);
                    end
%                     designMatContrastFreqBinRun(numFreqs*numContrasts+1,:) = ones;
                    %high-pass filter
                    designMatContrastFreqBinRun = highpassConcat(designMatContrastFreqBinRun, concatInfo{iSub,rwd}, r);
                    designMatContrastFreqBinRun(end+1,:) = ones;
                    designMatContrastFreqBin = [designMatContrastFreqBin designMatContrastFreqBinRun];
                    
                    
                    %null vs stim
                    clear designMatNullBinRun
                    for c=1:size(nullTrialMat{r,rwd},1)
                        temp = conv(nullTrialMat{r,rwd}(c,:),hrf);
                        designMatNullBinRun(c,:) =  temp(1:trialsPerRun*trialLength);
                    end
%                     designMatNullBinRun(2+1,:) = ones;
                    %high-pass filter
                    designMatNullBinRun = highpassConcat(designMatNullBinRun, concatInfo{iSub,rwd}, r);
                    designMatNullBinRun(end+1,:) = ones;
                    designMatNullBin = [designMatNullBin designMatNullBinRun];
                    
                end
                
                %perform GLM to get betas
                temp= designMatContrastBin'\tc;
                contrastBetas(iSub,iRoi,ibin,rwd,:) =temp(1:end-1);
                temp=designMatFreqBin'\tc;
                sfBetas(iSub,iRoi,ibin,rwd,:) = temp(1:end-1);
                temp=designMatContrastFreqBin'\tc;
                sfContrastBetas(iSub,iRoi,ibin,rwd,:)= temp(1:end-1);
                temp = designMatNullBin'\tc;
                nullBetas(iSub,iRoi,ibin,rwd,:) = temp(1:end-1); 
                
                %deconvolution with stim/null predictors
                temp = deconvMatNull{rwd}'\tc;
                temp = temp(1:end-1);
                temp=reshape(temp,deconvLength,[])';
                nullDeconv(iSub,iRoi,ibin,rwd,:,:) = temp;
%                 nullDeconv(iSub,iRoi,ibin,rwd,:) = deconvMatNull{rwd}'\tc;
                
                %deconvolution with contrast predictors
                temp = deconvMatContrast{rwd}'\tc;
                temp = temp(1:end-1);
                temp=reshape(temp,deconvLength,[])';
                contrastDeconv(iSub,iRoi,ibin,rwd,:,:) = temp;
%                 contrastDeconv(iSub,iRoi,ibin,rwd,:) = deconvMatContrast{rwd}'\tc;
                
                %deconvolution with frequency predictors
                temp = deconvMatFreq{rwd}'\tc;
                temp = temp(1:end-1);
                temp=reshape(temp,deconvLength,[])';
                freqDeconv(iSub,iRoi,ibin,rwd,:,:) = temp;
%                 freqDeconv(iSub,iRoi,ibin,rwd,:) = deconvMatFreq{rwd}'\tc;
                
                %deconvolution with frequency & contrast predictors
                temp = deconvMatContrastFreq{rwd}'\tc;
                temp = temp(1:end-1);
                temp=reshape(temp,deconvLength,[])';
                contrastFreqDeconv(iSub,iRoi,ibin,rwd,:,:) = temp;
%                 contrastFreqDeconv(iSub,iRoi,ibin,rwd,:) = deconvMatContrastFreq{rwd}'\tc;
                
            end
        end
     end
    
    deleteView(v);
    
    
    %% return to home directory
    cd('..');
end

%Get amplitudes of deconvolved responses, by taking std
nullDeconvStd = squeeze(std(nullDeconv,0,6));
contrastDeconvStd = squeeze(std(contrastDeconv,0,6));
freqDeconvStd = squeeze(std(freqDeconv,0,6));
contrastFreqDeconvStd = squeeze(std(contrastFreqDeconv,0,6));

%get amplitudes of deconvolved responses, by taking FFT amplitude
temp = fft(nullDeconv,deconvLength,6);
nullDeconvFftAmp = squeeze(abs(temp(:,:,:,:,:,2)));
nullDeconvFftPh = squeeze(angle(temp(:,:,:,:,:,2)));

temp = fft(contrastDeconv,deconvLength,6);
contrastDeconvFftAmp = squeeze(abs(temp(:,:,:,:,:,2)));
contrastDeconvFftPh = squeeze(angle(temp(:,:,:,:,:,2)));

temp = fft(freqDeconv,deconvLength,6);
freqDeconvFftAmp = squeeze(abs(temp(:,:,:,:,:,2)));
freqDeconvFftPh = squeeze(angle(temp(:,:,:,:,:,2)));

temp = fft(contrastFreqDeconv,deconvLength,6);
contrastFreqDeconvFftAmp = squeeze(abs(temp(:,:,:,:,:,2)));
contrastFreqDeconvFftPh = squeeze(angle(temp(:,:,:,:,:,2)));


%Amplitude of mean response, with std
trialTCstd = squeeze(std(binTrialTC,0,5));

%%

if nullTrialsPredictor
    nullTrialsString = 'onlyNulls';
else
    nullTrialsString = 'allTrials';
end
save([dataFolder 'rwdRapidData_regress_' nullTrialsString num2str(regressType) '.mat'], 'dataFolder', 'subFolders', 'samplerate', 'roiNames', ...
    'frames', 'junkedFrames', 'TR', 'trialsPerRun', 'trialLength', 'nVolumes',...
    'numContrasts', 'numFreqs', 'numSubs', 'numRois', 'onlyCorrect', ...
    'concatInfo', 'nullBetas', 'contrastBetas','sfBetas','sfContrastBetas', ...
    'deconvLength',...
    'contrastBetas', 'sfBetas', 'sfContrastBetas', 'nullBetas', ...
    'nullDeconv', 'contrastDeconv', 'freqDeconv', 'contrastFreqDeconv', ...
    'binMeanTseries','binTrialTC','numBinVoxels','binHrf',...
    'nullDeconvStd', 'contrastDeconvStd', 'freqDeconvStd', 'contrastFreqDeconvStd',...
    'nbins','binBorders','trialTCstd',...
    'nullDeconvFftAmp','nullDeconvFftPh','contrastDeconvFftAmp','contrastDeconvFftPh',...
    'freqDeconvFftAmp','freqDeconvFftPh','contrastFreqDeconvFftAmp','contrastFreqDeconvFftPh',...
    'regressBetasPca', 'pcaCorr', 'regressBetasGlobal', 'globalCorr','regressScaling');


nullTrialsPredictor%either only null trials (1) or all trials (0)
regressType%1=pca, 2=global mean
regressScaling
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
function filteredData = highpassConcat(concatData, concatInfo, irun)
%remove mean
% concatData = concatData - mean(concatData,2);
% filteredData = concatData;
%high-pass filter
filteredData = real(ifft(fft(concatData') .* repmat(concatInfo.hipassfilter{irun}', 1, size(concatData,1)) ))';

end

%%
function newConcat = zscoreConcat(concatData, concatInfo)%voxels,TRs
    newConcat = [];
    for iRun=1:size(concatInfo.runTransition,1)
                thisRun = concatData(:,concatInfo.runTransition(iRun,1):concatInfo.runTransition(iRun,2));
%         nanmean(std(thisRun'))
%         nanmean(std(zscore(thisRun,0,2)'))
    newConcat = [newConcat zscore(thisRun,0,2)];
%         newConcat = [newConcat zscore(concatData(:,concatInfo.runTransition(iRun,1):concatInfo.runTransition(iRun,2))')'];
    end
%     figure(1); 
%     clf
%     subplot(1,2,1); plot(nanmean(newConcat)); 
%     subplot(1,2,2); plot(nanmean(concatData));
%     nanmean(std(newConcat'))
%     keyboard
end
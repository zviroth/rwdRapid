%%
close all;
clear all;
dataFolder = '/Volumes/MH02086153MACDT-Drobo/rwdFmri/';

load([dataFolder 'rwdRapid_microsaccade.mat'], 'dataFolder', 'subFolders', 'samplerate',  ...
    'trialsPerRun', 'trialLength', ...
    'numSubs', 'onlyCorrect', ...
    'sacRwd','binnedSac', 'smoothSac', ...
    'rwdPupil','meanPupil','taskTimes',...
    'sacRun','pupilRun', 'numRuns','startTimes','endTimes','eyetrackerTime',...
    'nullTrials','meanNullPupil','meanStimPupil',...
    'sacRunSmooth', 'sacRunSmoothTrialZeroFilled','sacRunSmoothTrial','L',...
    'smoothSac', 'sacPupilFuncRun', 'sacPupilFunc','sumBinnedSac','sumSmoothSac','allMeanPupil', ...
    'allResidMeanPupil', 'allMeanStimPupil','allResidStimPupil','allMeanNullPupil','allResidNullPupil','allRegressMeanPupil',...
    'concatSac','residRwdPupil','regressRwdPupil','residRun','sacPupilRun','residMeanPupil','meanResidNullPupil','meanResidNullPupil');
keyboard

%%
rwdTypes = {'high','low'};
figure(1); clf%pupil size

rows=4;
cols=numSubs;
for iSub=1:numSubs
    minPup = min(min(meanPupil{iSub,1}),min(meanPupil{iSub,2})) - 100;
    maxPup = max(max(meanPupil{iSub,1}),max(meanPupil{iSub,2})) + 100;
    for rwd=1:2
        subplot(rows,cols,iSub + (rwd-1)*cols)
        plot(rwdPupil{iSub,rwd}')
        hold all
        plot(meanPupil{iSub,rwd}, 'color', plotColors{rwd}, 'linewidth', 2)
        
        titleString = [rwdTypes{rwd} ' reward: ' num2str(size(rwdPupil{iSub,rwd},1)) ' trials'];
        title(titleString);
%             ylim([minPup maxPup]);
        plotRwdStimSegments(meanPupil{iSub,rwd}, taskTimes(iSub,:))
        
        xlabel('time (ms)')
        ylabel('pupil size (pixels)');
        xticks(0:1000:8000);
        xticklabels(2*[0:1000:8000]);
    end

    subplot(rows,cols,iSub+2*cols)
    for rwd=1:2
        plot(meanPupil{iSub,rwd}, 'Color', plotColors{rwd}, 'linewidth', 1)
        hold on
    end
    legend('high', 'low')
    xlabel('time (ms)')
    ylabel('pupil size (pixels)');
    ylim([minPup maxPup]);
    xticks(0:1000:8000);
    xticklabels(2*[0:1000:8000]);
    subplot(rows,cols,iSub+3*cols)

    l=subMeanPupilLength(iSub)-100;
    diffPupil = meanPupil{iSub,1}(1:l) - meanPupil{iSub,2}(1:l);
    plot(zeros(l,1),'--');
    hold all
    plot(diffPupil, 'color', 'k', 'linewidth', 1)
    plotRwdStimSegments(diffPupil, taskTimes(iSub,:))
    title('high minus low');
    xlabel('time (ms)')
    xticks(0:1000:8000);
    xticklabels(2*[0:1000:8000]);
    ylabel('difference in pupil size (pixels)');
    set(gcf, 'position', [130 130 600 1000]);
end

figure(2); clf%microsaccades
rows=2;
cols=numSubs;
% L = 100;
% filter = ones(L,1);
% % filter = [1:L/2 L/2:-1:1];
% filter = filter./sum(filter);
for iSub=1:numSubs
    subplot(rows,cols,iSub)
    for rwd=1:2
        plot(binnedSac{iSub,rwd})
        hold on
    end
    ylim([0 20]);
    xlabel('time (ms)')
    ylabel('# of microsaccade onsets');
    title('microsaccade onsets: original data');
    xticks(0:4000:8000);
    xticklabels(2*[0:4000:8000]);
    
    subplot(rows,cols,iSub + cols)
    for rwd=1:2
%         smoothSac{iSub,rwdTypeNum} = conv(binnedSac{iSub,rwdTypeNum}, filter, 'same');
        plot(smoothSac{iSub,rwd},'linewidth', 1)
        hold on
    end
    legend('high', 'low')
    xlabel('time (ms)')
    ylabel('# of microsaccade onsets');
    title('microsaccade onsets: smoothed data');
    xticks(0:4000:8000);
    xticklabels(2*[0:4000:8000]);
    set(gcf, 'position', [100 100 600 500]);
    
end


%%
figure(3); clf
rows=2;
cols=numSubs;
for iSub=1:numSubs
    subplot(rows,cols,iSub)
    for rwd=1:2
        for r=1:numRuns(iSub)
            plot(squeeze(sacPupilFuncRun{iSub}(rwd,r,1:end-1)) , 'color',plotColors{rwd},'linewidth',0.5);
            hold all
        end
        plot(squeeze(mean(sacPupilFuncRun{iSub}(rwd,:,1:end-1))) , 'color',plotColors{rwd},'linewidth',2);
    end
    temp=squeeze(mean(sacPupilFuncRun{iSub}(:,:,1:end-1)));
    plot(squeeze(mean(temp)) , 'color','k','linewidth',3);
    yl = ylim;
    
    subplot(rows,cols,iSub+cols)
    plot(sacPupilFunc(iSub,1:end-1) , 'linewidth',2);% ignore constant
    ylim(yl);
end

%%
figure(4); clf
rows=2;
cols=2;
subplot(rows,cols,1)
plot(sumBinnedSac');
title('mean microsaccade onsets');
subplot(rows,cols,2)
plot(sumSmoothSac');
title('smoothed mean ms rate');

subplot(rows,cols,3)
plot(squeeze(nanmean(sacPupilFunc(:,1:end-1))),'linewidth', 3);%ignore constant
title('pupil-microsaccade function');

clear meanConv
for rwd=1:2
    temp = conv(sumSmoothSac(rwd,:)',squeeze(nanmean(sacPupilFunc(:,1:end-1))));
    meanConv(rwd,:) = temp(1:length(sumSmoothSac));
end

subplot(rows,cols,4)
plot(meanConv')
title('pupil-ms function convolved with mean ms rate');

%%
rwdTypes = {'high','low'};
figure(5); clf%pupil size

rows=4;
cols=numSubs;
for iSub=1:numSubs
    minPup = min(min(residMeanPupil{iSub,1}),min(residMeanPupil{iSub,2})) - 100;
    maxPup = max(max(residMeanPupil{iSub,1}),max(residMeanPupil{iSub,2})) + 100;
    for rwd=1:2
        subplot(rows,cols,iSub + (rwd-1)*cols)
        plot(residRwdPupil{iSub,rwd}')
        hold all
        plot(residMeanPupil{iSub,rwd}, 'color', plotColors{rwd}, 'linewidth', 2)
        
        titleString = [rwdTypes{rwd} ' reward: ' num2str(size(residRwdPupil{iSub,rwd},1)) ' trials'];
        title(titleString);
%             ylim([minPup maxPup]);
        plotRwdStimSegments(residMeanPupil{iSub,rwd}, taskTimes(iSub,:))
        
        xlabel('time (ms)')
        ylabel('pupil size (pixels)');
        xticks(0:1000:8000);
        xticklabels(2*[0:1000:8000]);
    end

    subplot(rows,cols,iSub+2*cols)
    for rwd=1:2
        plot(residMeanPupil{iSub,rwd}, 'Color', plotColors{rwd}, 'linewidth', 1)
        hold on
    end
    legend('high', 'low')
    xlabel('time (ms)')
    ylabel('pupil size (pixels)');
    ylim([minPup maxPup]);
    xticks(0:1000:8000);
    xticklabels(2*[0:1000:8000]);
    subplot(rows,cols,iSub+3*cols)
    l = min(length(residMeanPupil{iSub,1}), length(residMeanPupil{iSub,2}));
    l=l-100;
    diffPupil = residMeanPupil{iSub,1}(1:l) - residMeanPupil{iSub,2}(1:l);
    plot(zeros(l,1),'--');
    hold all
    plot(diffPupil, 'color', 'k', 'linewidth', 1)
    plotRwdStimSegments(diffPupil, taskTimes(iSub,:))
    title('high minus low');
    xlabel('time (ms)')
    xticks(0:1000:8000);
    xticklabels(2*[0:1000:8000]);
    ylabel('difference in pupil size (pixels)');
    set(gcf, 'position', [130 130 600 1000]);
end

%%

figure(6)
clf
rows=2;
cols=5;

subplot(rows,cols,1)
plot(allMeanPupil') %original mean pupil size
title('original mean pupil')
miny = min(allMeanPupil(:)) -400;
maxy = max(allMeanPupil(:)) + 400;
y = [miny maxy];
ylim(y);

subplot(rows,cols,cols+1)
plot(allResidMeanPupil') %mean pupil size after removing microsaccades convolved with pupil function
title('pupil after regressing out ms')
ylim(y);
% ylim([miny maxy]);


subplot(rows,cols, 2)
plot(allMeanPupil(1,:)'-allMeanPupil(2,:)');%difference between high and low
title('high minus low')
y2 = ylim;

subplot(rows,cols, cols+2)
plot(allResidMeanPupil(1,:)'-allResidMeanPupil(2,:)');%difference between high and low after removing ms
title('high minus low after regression')
ylim(y2);

subplot(rows,cols,3)
plot(allMeanStimPupil') %original mean pupil size for stim trials
title('original mean stim')
ylim(y);

subplot(rows,cols, cols+3)
plot(allResidStimPupil') %mean pupil size for stim trials after removing microsaccades
title('mean stim after regression')
ylim(y);

subplot(rows,cols,4)
plot(allMeanNullPupil') %original mean pupil size for null trials
title('original mean null')
ylim(y);

subplot(rows,cols, cols+4);
plot(allResidNullPupil') %mean pupil size for null trials after removing microsaccades
title('mean null after regression')
ylim(y);

subplot(rows,cols,rows*cols)
plot(allRegressMeanPupil') %what we regressed out
title('mean regressed out')




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
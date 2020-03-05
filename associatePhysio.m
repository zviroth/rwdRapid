%stimfiles must be already associated
clear d d2 d3 
v=newView;

v = viewSet(v,'groupname','Raw');
nScans = viewGet(v, 'nscans');

% %find logfile for mepreproc
% temp = dir([pwd '/afnipreproc/output.proc.*']);
% preprocfilename = fullfile(temp.folder, temp.name);
% [fileID msg] = fopen(preprocfilename,'r');
% filetext = fileread(preprocfilename);
% expr = 'chan_001.nii[4.';





% curfolder = getLastDir(pwd);
% % f = regexp(pwd,'/'); p=pwd; curfolder = p(f(end)+1:end);
% expr = ['3dTcat -prefix ' curfolder '.results/pb00.' curfolder '.r*.e01.tcat mr_00*_'];
% expr = ['3dTcat -prefix ' curfolder '.results/pb00.' curfolder '.r'];
% [startInd endInd] = regexp(filetext,expr);%3dTcat -prefix 005520190411.results/pb00.005520190411.r01.e01.tcat mr_0004_task_chan_001.nii[4..$]
% for i=1:length(startInd)
%    preprocnewrun(i) = str2num(filetext(endInd(i)+1:endInd(i)+2));
%    preprocoldrun(i) = str2num(filetext(endInd(i)+18:endInd(i)+19));
% end

if isfile('origFilenames.mat')
    % load file with original filenames for all Raw timeseries (excluding blip)
    load('origFilenames.mat');
    
    %find realtime directory
    d=dir('*.tar');
    d2 = dir([d.name(1:12) '*']);%beginning of the name of the zipped data file
    d3 = dir(d2(1).name);%should include only one folder
    realtimePath = fullfile(pwd, d2(1).name, d3(end).name, 'realtime/');
    
    for iscan = 1:nScans
        s = viewGet(v, 'stimfile', iscan);
        filename = viewGet(v,'tseriesfile', iscan);
        if strcmp(filename(1),'_')%filename is uncorrected, doesn't include the original
            numrun = str2num(filename(9:10));
            origName = listOfOriginalFilenames{numrun};
            nameStart = strfind(origName, 'mr_0');
            origNum = str2num(origName(nameStart+5:nameStart+6));
            
        elseif strcmp(filename(1:5),'mr_00')%filename has been correcly named
            %this can be changed - we can get the run from the filename, no
            %need for listOfOriginalFilenames
            k = strfind(filename,'_aproc');
            numrun = str2num(filename(k+8:k+9));
            origName = listOfOriginalFilenames{numrun};
            nameStart = strfind(origName, 'mr_0');
            origNum = str2num(origName(nameStart+5:nameStart+6));
            
        else% runs imported from another session ???
            keyboard
            %need to run this on the original session, then copy the
            %stimfiles
        end
        origNum = origNum+1;
        origName(6:7) = num2str(origNum,'%02u');
        ecgfilename = ['ECG_epiRTnih_scan_' origName(4:7)];
        d = dir([realtimePath '*scan_' origName(4:7) '*']);%d(1) is ecg, d(2) is resp
        
        s{1}.myscreen.originalmlrdir = pwd;
        s{1}.myscreen.physiodir = fullfile(d2(1).name, d3(end).name, 'realtime/');
        %     s{1}.myscreen.ecgfilename = fullfile(d2(1).name, d3(end).name, 'realtime', d(1).name);
        %     s{1}.myscreen.respirationfilename = fullfile(d2(1).name, d3(end).name, 'realtime', d(2).name);
        s{1}.myscreen.ecgfilename = d(1).name;
        s{1}.myscreen.respirationfilename = d(2).name;
        stimfilename = viewGet(v,'stimfilename',iscan);
        [FILEPATH,NAME,EXT] = fileparts(stimfilename{1});
        %     save(stimfilename{1},s{1}.myscreen, s{1}.task, s{1}.
        temp = s{1};
%         save(['~/temp/' [NAME EXT]],'-struct','temp');

        save(stimfilename{1},'-struct','temp');
    end
else
    %origFilenames.mat doesn't exist
    keyboard
end
deleteView(v);
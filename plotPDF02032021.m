addpath('../fileSet');
addpath('../iniconfig');
addpath('../cbrewer');

%% get list of files

try
    configFileName;
catch
    % config file path
    configFileName = '20210203_DNA_FISH_plot_config.ini';
end

% create file list object
fs = fileSet;
fs.buildFileSetFromConfig(configFileName);

% reformat channel, replicate and FOV as numbers
convertNonNumeralStrings = 1; 
fs.toNumeric({'Channel','TechnicalReplicate','FOV'},convertNonNumeralStrings);

%read parameters from config file relative to analysis into a parameters object
params = readParamsFromIniFile(configFileName);

%% Plot PDF and CDF for dis between airlocalize center

%find list of conditions
cond = find_cond(fs.fList);
colorID = [];

if params.Settings.testing
    cond = 1;
end

if params.channelDescription.ConditionDescription
    %group conditions and generate color scheme for ploting
    des = parse_description(params.channelDescription.ConditionDescription);
    cond = group_conditions(cond,des);
    size0 = size(cond);
    color1 = color_generator(size0);
    idx = 1:numel(cond);
    idx = reshape(idx,size0);
    [cond,idx1] = flat_array(cond,'column',string(idx));
    [row,col]=ind2sub(size0,str2double(idx1));
    colorID = zeros(length(idx1),3);
    for i = 1:length(idx1)
        colorID(i,:) = color1(row(i),col(i),:);
    end
    cond = mat2cell(cond,ones(length(cond),1));
end

%read in dis
dist = {};
for i = 1:numel(cond(:,1))
    if params.Settings.testing
        curFileList = fs.fList{:,4};
    else
        curFileList = fs.getFileName({'Condition'},cond(i),'match','all');
    end
    
    dis = [];
    for j = 1:length(curFileList)
        %load match file
        match = load(curFileList{j},'-ascii');
        if isempty(match)
            continue
        end
        %save dist
        dis = [dis;match(:,6)];    
    end
    dist{i} = dis;%each column is dis of a different condition
end

plot_fig(dist,cond,'distance/nm',colorID,params.Settings.hist_binSize,params.Output.outFolder,'dis_Airlocalize')

%% Plot PDF and CDF for dx, dy, dz of airlocalize center
voxsize = [params.Settings.voxSize_dxy,...
            params.Settings.voxSize_dxy,...
            params.Settings.voxSize_dz];
        
try
    cond;
    colorID;
catch
    %find list of conditions
    cond = find_cond(fs.fList);
    colorID = [];
    if params.channelDescription.ConditionDescription
        des = parse_description(params.channelDescription.ConditionDescription);
        cond = group_conditions(cond,des);
        size0 = size(cond);
        color1 = color_generator(size0);
        idx = 1:numel(cond);
        idx = reshape(idx,size0);
        [cond,idx1] = flat_array(cond,'column',string(idx));
        [row,col]=ind2sub(size0,str2double(idx1));
        colorID = zeros(length(idx1),3);
        for i = 1:length(idx1)
            colorID(i,:) = color1(row(i),col(i),:);
        end
        cond = mat2cell(cond,ones(length(cond),1));
    end
end

%read in dis
dx = {};
dy = {};
dz = {};
for i = 1:numel(cond(:,1))
    if params.Settings.testing
        curFileList = fs.fList{:,4};
    else
        curFileList = fs.getFileName({'Condition'},cond(i),'match','all');
    end
    dis = [];
    for j = 1:length(curFileList)
        %load match file
        match = load(curFileList{j},'-ascii');
        if isempty(match)
            continue
        end
        %save dist
        dis = [dis;match(:,8:10)];    
    end
    if isempty(dis)
        dx{i} = [];
        dy{i} = [];
        dz{i} = [];
        continue
    end
    dis = int64(convert_loc_pix_to_nm(dis,1./voxsize));
    dx{i} = double(dis(:,1));%each column is dis of a different condition
    dy{i} = double(dis(:,2));
    dz{i} = double(dis(:,3));
end

plot_fig(dx,cond,'pixel',colorID,1,params.Output.outFolder,'dx')
plot_fig(dy,cond,'pixel',colorID,1,params.Output.outFolder,'dy')
plot_fig(dz,cond,'pixel',colorID,1,params.Output.outFolder,'dz')

%% plot PDF and CDF for overlap

try
    cond;
    colorID;
catch
    %find list of conditions
    cond = find_cond(fs.fList);
    colorID = [];
    if params.channelDescription.ConditionDescription
        des = parse_description(params.channelDescription.ConditionDescription);
        cond = group_conditions(cond,des);
        size0 = size(cond);
        color1 = color_generator(size0);
        idx = 1:numel(cond);
        idx = reshape(idx,size0);
        [cond,idx1] = flat_array(cond,'column',string(idx));
        [row,col]=ind2sub(size0,str2double(idx1));
        colorID = zeros(length(idx1),3);
        for i = 1:length(idx1)
            colorID(i,:) = color1(row(i),col(i),:);
        end
        cond = mat2cell(cond,ones(length(cond),1));
    end
end

%read in overlap
dist = {};
for i = 1:numel(cond(:,1))
    if params.Settings.testing
        curFileList = fs.fList{:,4};
    else
        curFileList = fs.getFileName({'Condition'},cond(i),'match','all');
    end
    dis = [];
    for j = 1:length(curFileList)
        %load overlap file
        [~,f,~] = fileparts(curFileList{j});
        try
            overlap = load(fullfile(params.InputFolders.OverlapFolder,[f,'_overlap.txt']),'-ascii');
        catch
            continue
        end
        %save dist
        dis = [dis;round(overlap,10)];    
    end
    if isempty(dis)
        dist{i} = [];
        continue
    end
    dis0 = 1./dis;
    dis0(isinf(dis0)) = 100;
    dis0(dis0>100) = 100;
    dist{i} = dis0;%each column is dis of a different condition
end

plot_fig(dist,cond,'1/correlation',colorID,0.05,params.Output.outFolder,'corr')

%% Plot PDF and CDF for dis between loci center
voxsize = [params.Settings.voxSize_dxy,...
            params.Settings.voxSize_dxy,...
            params.Settings.voxSize_dz];

try
    cond;
    colorID;
catch
    %find list of conditions
    cond = find_cond(fs.fList);
    colorID = [];
    if params.channelDescription.ConditionDescription
        des = parse_description(params.channelDescription.ConditionDescription);
        cond = group_conditions(cond,des);
        size0 = size(cond);
        color1 = color_generator(size0);
        idx = 1:numel(cond);
        idx = reshape(idx,size0);
        [cond,idx1] = flat_array(cond,'column',string(idx));
        [row,col]=ind2sub(size0,str2double(idx1));
        colorID = zeros(length(idx1),3);
        for i = 1:length(idx1)
            colorID(i,:) = color1(row(i),col(i),:);
        end
        cond = mat2cell(cond,ones(length(cond),1));
    end
end

%calculate dis between center
dist = {};
distx = {};
disty = {};
distz = {};
for i = 1:numel(cond(:,1))
    if params.Settings.testing
        curFileList = fs.fList{:,4};
    else
        curFileList = fs.getFileName({'Condition'},cond(i),'match','all');
    end
    dis1 = [];
    dis2 = [];
    for j = 1:length(curFileList)
        %load match file
        match = load(curFileList{j},'-ascii');        
        [d,f,~] = fileparts(curFileList{j});
        disp('processing file:')
        disp(f)
        if isempty(match)
            continue
        end
        
        % find the conditions in the current file
        fl = fs.getFileName({},{},'match','all');
        idx = find( ismember(fl,curFileList{j}) );
        curCond = {};
        for k = 1:numel(fs.conditions)
            curCond = [curCond,fs.fList.(fs.conditions{k})(idx)];
        end
        %load loc3 of centroid
        C1 = curCond{1};
        C2 = curCond{2};
        CropFolder1 = fullfile(params.InputFolders.CropFolder,...
            strcat('C',C1,f(5:length(f)),'_crop'));
        CropFolder2 = fullfile(params.InputFolders.CropFolder,...
            strcat('C',C2,f(5:length(f)),'_crop'));
        cen1 = load(fullfile(CropFolder1,'centroid0.loc3'),'-ascii');
        cen2 = load(fullfile(CropFolder2,'centroid0.loc3'),'-ascii');
        deltanm = match(:,3:5);
        delta = convert_loc_pix_to_nm(deltanm,1./voxsize);
        
        dis = [];
        
        for k = 1:length(match(:,1))
            dis3 = cen2(match(k,2),:)-cen1(match(k,1),:);
            dis2 = [dis2;dis3];
            dis3 = convert_loc_pix_to_nm(dis3,voxsize);
            dis0 = sqrt(sum(dis3.^2));
            %save dist
            dis = [dis;dis0];
        end
        save(fullfile(d,strcat(f,'_cen_dist.txt')),'dis','-ascii')
        dis1 = [dis1;dis];
    end
    if isempty(dis1)
        dist{i} = [];
        distx{i} = [];
        disty{i} = [];
        distz{i} = [];
        continue
    end
    dist{i} = dis1;%each column is dis of a different condition
    distx{i} = dis2(:,1);
    disty{i} = dis2(:,2);
    distz{i} = dis2(:,3);
end

%plot CDF and PDF
plot_fig(dist,cond,'distance/nm',colorID,params.Settings.hist_binSize,params.Output.outFolder,'dis_cen')
plot_fig(distx,cond,'pixel',colorID,1,params.Output.outFolder,'dx_cen')
plot_fig(disty,cond,'pixel',colorID,1,params.Output.outFolder,'dy_cen')
plot_fig(distz,cond,'pixel',colorID,1,params.Output.outFolder,'dz_cen')

%% Plot PDF and CDF for Rg and volumn
voxsize = [params.Settings.voxSize_dxy,...
            params.Settings.voxSize_dxy,...
            params.Settings.voxSize_dz];
        
try
    cond;
    colorID;
catch
    %find list of conditions
    cond = find_cond(fs.fList);
    colorID = [];
    if params.channelDescription.ConditionDescription
        des = parse_description(params.channelDescription.ConditionDescription);
        cond = group_conditions(cond,des);
        size0 = size(cond);
        color1 = color_generator(size0);
        idx = 1:numel(cond);
        idx = reshape(idx,size0);
        [cond,idx1] = flat_array(cond,'column',string(idx));
        [row,col]=ind2sub(size0,str2double(idx1));
        colorID = zeros(length(idx1),3);
        for i = 1:length(idx1)
            colorID(i,:) = color1(row(i),col(i),:);
        end
        cond = mat2cell(cond,ones(length(cond),1));
    end
end

R1 = {};
R2 = {};
R1x = {};
R1y = {};
R1z = {};
R2x = {};
R2y = {};
R2z = {};
V1 = {};
V2 = {};
I1 = {};
I2 = {};
Vi1 ={};
Vi2 ={};
for i = 1:numel(cond(:,1))
    if params.Settings.testing
        curFileList = fs.fList{:,4};
    else
        curFileList = fs.getFileName({'Condition'},cond(i),'match','all');
    end
    r1 = [];
    r2 = [];
    v1 = [];
    v2 = [];
    i1 = [];
    i2 = [];
    for j = 1:length(curFileList)
        [~,f,~] = fileparts(curFileList{j});
        % find the conditions in the current file
        fl = fs.getFileName({},{},'match','all');
        idx = find( ismember(fl,curFileList{j}) );
        curCond = {};
        for k = 1:numel(fs.conditions)
            curCond = [curCond,fs.fList.(fs.conditions{k})(idx)];
        end
        %load Rg and volume
        C1 = curCond{1};
        C2 = curCond{2};
        CropFolder1 = fullfile(params.InputFolders.CropFolder,...
            strcat('C',C1,f(5:length(f)),'_crop'));
        CropFolder2 = fullfile(params.InputFolders.CropFolder,...
            strcat('C',C2,f(5:length(f)),'_crop'));            
        try
            r1 = [r1;load(fullfile(CropFolder1,'Rg.loc3'),'-ascii')];
            r2 = [r2;load(fullfile(CropFolder2,'Rg.loc3'),'-ascii')];
            v1 = [v1;load(fullfile(CropFolder1,'volumn.txt'),'-ascii')];
            v2 = [v2;load(fullfile(CropFolder2,'volumn.txt'),'-ascii')];
            loc = load(fullfile(CropFolder1,'location_mask_sorted.loc3'),'-ascii');
            i1 = [i1;loc(:,5)];
            loc = load(fullfile(CropFolder2,'location_mask_sorted.loc3'),'-ascii');
            i2 = [i2;loc(:,5)];
        catch
            continue
        end
    end
    if isempty(r1)
        R1x{i} = [];
        R1y{i} = [];
        R1z{i} = [];
        R1{i} = [];
    else
        R1x{i} = r1(:,1);
        R1y{i} = r1(:,2);
        R1z{i} = r1(:,3);
        R1{i} = sqrt(sum(convert_loc_pix_to_nm(r1,voxsize).^2,2));
    end
    if isempty(r2)
        R2x{i} = [];
        R2y{i} = [];
        R2z{i} = [];
        R2{i} = [];
    else
        R2x{i} = r2(:,1);
        R2y{i} = r2(:,2);
        R2z{i} = r2(:,3);
        R2{i} = sqrt(sum(convert_loc_pix_to_nm(r2,voxsize).^2,2));
    end
    V1{i} = v1;
    V2{i} = v2;
    I1{i} = i1;
    I2{i} = i2;
    Vi1{i} = v1./i1;
    Vi2{i} = v2./i2;
end
plot_fig(R1,cond,'Rg/nm',colorID,0.5,params.Output.outFolder,'Rg_127')
plot_fig(R2,cond,'Rg/nm',colorID,0.5,params.Output.outFolder,'Rg_130')
plot_fig(R1x,cond,'Rg/pixel',colorID,0.02,params.Output.outFolder,'Rg_x_127')
plot_fig(R2x,cond,'Rg/pixel',colorID,0.02,params.Output.outFolder,'Rg_x_130')
plot_fig(R1y,cond,'Rg/pixel',colorID,0.02,params.Output.outFolder,'Rg_y_127')
plot_fig(R2y,cond,'Rg/pixel',colorID,0.02,params.Output.outFolder,'Rg_y_130')
plot_fig(R1z,cond,'Rg/pixel',colorID,0.02,params.Output.outFolder,'Rg_z_127')
plot_fig(R2z,cond,'Rg/pixel',colorID,0.02,params.Output.outFolder,'Rg_z_130')
plot_fig(V1,cond,'Volume/pixel',colorID,params.Settings.hist_binSize,params.Output.outFolder,'Volume_127')
plot_fig(V2,cond,'Volume/pixel',colorID,params.Settings.hist_binSize,params.Output.outFolder,'Volume_130')

% comparing Rg and V
figure;
for i = 1:length(cond)
    h = scatter(cell2mat(R1(:,i)),cell2mat(V1(:,i)));
    try
        set(h,'color',colorID(i,:))             
    catch
    end
    set(h,'LineWidth',1.0)
    hold on
end
legend(cond)
xlabel('Rg/nm')
ylabel('Volume/pixel')
hold off
savefig(fullfile(params.Output.outFolder,'Rg_Volume_127.fig'))

figure;
for i = 1:length(cond)
    h = scatter(cell2mat(R2(:,i)),cell2mat(V2(:,i)));
    try
        set(h,'color',colorID(i,:))             
    catch
    end
    set(h,'LineWidth',1.0)
    hold on
end
legend(cond)
xlabel('Rg/nm')
ylabel('Volume/pixel')
hold off
savefig(fullfile(params.Output.outFolder,'Rg_Volume_130.fig'))

% plot Rg against intensity
figure;
for i = 1:length(cond)
    h = scatter(cell2mat(I1(:,i)),cell2mat(R1(:,i)));
    try
        set(h,'color',colorID(i,:))             
    catch
    end
    set(h,'LineWidth',1.0)
    hold on
end
legend(cond)
ylabel('Rg/nm')
xlabel('Intensity')
hold off
savefig(fullfile(params.Output.outFolder,'Rg_Intensity_127.fig'))

figure;
for i = 1:length(cond)
    h = scatter(cell2mat(I2(:,i)),cell2mat(R2(:,i)));
    try
        set(h,'color',colorID(i,:))             
    catch
    end
    set(h,'LineWidth',1.0)
    hold on
end
legend(cond)
ylabel('Rg/nm')
xlabel('Intensity')
hold off
savefig(fullfile(params.Output.outFolder,'Rg_Intensity_130.fig'))

% plot V against intensity
figure;
for i = 1:length(cond)
    h = scatter(cell2mat(I1(:,i)),cell2mat(V1(:,i)));
    try
        set(h,'color',colorID(i,:))             
    catch
    end
    set(h,'LineWidth',1.0)
    hold on
end
legend(cond)
ylabel('Volume/pixel')
xlabel('Intensity')
hold off
savefig(fullfile(params.Output.outFolder,'Volume_Intensity_127.fig'))

figure;
for i = 1:length(cond)
    h = scatter(cell2mat(I2(:,i)),cell2mat(V2(:,i)));
    try
        set(h,'color',colorID(i,:))             
    catch
    end
    set(h,'LineWidth',1.0)
    hold on
end
legend(cond)
ylabel('Volume/pixel')
xlabel('Intensity')
hold off
savefig(fullfile(params.Output.outFolder,'Volume_Intensity_130.fig'))

% plot Vi against Rg
figure;
for i = 1:length(cond)
    h = scatter(cell2mat(R1(:,i)),cell2mat(Vi1(:,i)));
    try
        set(h,'color',colorID(i,:))             
    catch
    end
    set(h,'LineWidth',1.0)
    hold on
end
legend(cond)
ylabel('Volume/Intensity')
xlabel('Rg/nm')
hold off
savefig(fullfile(params.Output.outFolder,'Rg_Vi_127.fig'))

figure;
for i = 1:length(cond)
    h = scatter(cell2mat(R2(:,i)),cell2mat(Vi2(:,i)));
    try
        set(h,'color',colorID(i,:))             
    catch
    end
    set(h,'LineWidth',1.0)
    hold on
end
legend(cond)
ylabel('Volume/Intensity')
xlabel('Rg/nm')
hold off
savefig(fullfile(params.Output.outFolder,'Rg_Vi_130.fig'))

%% plot cenroid dis v.s. airocalize dis
%read in air localize dis
%if params.Settings.testing
%    curFileList = fs.fList{:,4};
%else
%    curFileList = fs.getFileName({},{},'match','all');
%end
    
%dis = [];
%for j = 1:length(curFileList)
    %load match file
%    match = load(curFileList{j},'-ascii');
%    if isempty(match)
%        continue
%    end
    %save dist
%    dis = [dis;match(:,6)];    
%end

%read in centroid dis
%disc = [];
%for j = 1:length(curFileList)
    %load distance
%    [d,f,~] = fileparts(curFileList{j});
%    dist = load(fullfile(d,strcat(f,'_cen_dist.txt')),'-ascii');
    %save dist
%    disc = [disc;dist];    
%end

%figure;
%scatter(dis,disc)
%xlabel('Airlocalize dis/nm')
%ylabel('Centroid dis/nm')

%% for testing
%curFileList = fs.getFileName({},{},'match','all');
%cen = [];
%dis = [];
%for j = 1:length(curFileList)
%    [~,f,~] = fileparts(curFileList{j});
    %load file
%    cen = [cen;load(fullfile('R:\lionnt01lab\lionnt01labspace\Yi_Fu\test\100_noStd_noNoise\res\loci_match',strcat(f,'_cen_dist.txt')))];
    %dis = [dis;load(fullfile('R:\lionnt01lab\lionnt01labspace\Yi_Fu\test\100_noStd_noNoise\res\',strcat('C2-',num2str(j),'_crop'),'location.loc3'))];
%end
%figure;
%hold on
%scatter(dis(:,1)-floor(dis(:,1)),dis(:,1)-cen(:,1))
%scatter(dis(:,2)-floor(dis(:,2)),dis(:,2)-cen(:,2))
%scatter(dis(:,3)-floor(dis(:,3)),dis(:,3)-cen(:,3))

%% functions for plot
function cond = find_cond(fList)
    %find list of conditions
    cond = unique(fList(:,3));
    cond = cond.(1);%convert table to array
end
function cond2=group_conditions(cond,key)
    cond1 = strings(length(cond),length(key));
    k = zeros(1,length(key));
    for i = 1:length(cond)
        for j = 1:length(key)
            if contains(cond(i),key(j))
                cond1(k(j)+1,j) = string(cond(i));
                k(j) = k(j) + 1;
            end
        end
    end

    %remove empty strings
    cond1(all(strcmp(cond1,""),2),:) = [];
    
    %sort by drug concentration
    cond2 = strings(size(cond1));
    for i = 1:length(key)
        con = regexp(cond1(:,i),'\d+(?:\.\d+)?uM','match');
        con = double(erase([con{:}],'uM'));
        [~,idx] = sort(con);
        for j = 1:k(i)
            if contains(cond1(j,i),'ctrl')%only 1 ctrl for each cond
                cond2(1,i)=cond1(j,i);
            else
                cond2(find(idx==j,1)+1,i)=cond1(j,i);
            end
        end
        cond = setdiff(cond1(:,i),cond2);
        cond2(length(idx)+2:length(idx)+length(cond)+1,i)=cond;%add all cond don't cotain number or ctrl at the end
    end
end
function colorID = color_generator(size)
    %generate an array of colors, each column is a different series of
    %colors(max 6), each row has a different shade from darker to brighter
    if size(2) > 10
        disp('Only support 10 color series!')
        colorID = [];
        return
    end
    colorID = zeros(size(1),size(2),3);
    scheme = ["Blues","Oranges","Greens","Purples","Greys","Reds","RdPu","PuBuGn","PuRd","BuPu"];
    for i = 1:size(2)
        colors = cbrewer('seq',scheme(i),size(1)+2);
        %discard 2 most light colors and invert, so that darker colors come first
        colorID(:,i,:) = flip(colors(3:length(colors),:),1);
    end
end
function [array, array2] = flat_array(array,direction,array2)
    %flat string array aong specified direction and remove ""
    %Also manipulate array2 as did array
    %Array2 must be same size as array
    
    try
        if size(array2) ~= size(array)
            disp('Sizes of arrays don''t match!')
            return
        end
    catch
        disp('Dimensions of arrays don''t match!')
        return
    end
    if strcmp(direction,'column')
        array = reshape(array,[],1);
        array2 = reshape(array2,[],1);
     
    else
        if strcmp(direction,'row')
            array = reshape(array,1,[]);
            array2 = reshape(array2,1,[]);
        else
            disp('Direction must be by ''column'' or by ''row''.')
            return
        end        
    end
    array2 = array2((~strcmp(array,""))&(~ismissing(array)));
    array = array((~strcmp(array,""))&(~ismissing(array)));%remove "" and <missing>
end
function plot_fig(varible,cond,label,colorID,binsize,outFolder,figname)
    %This function group varible by conditions and plot a PDF and a CDF for
    %each condition in two plots (one PDF, one CDF)
    %color must be a matrix of RGB value with each row as a color
    %varible is an array with each column content results from a different
    %condition
    %figname is name of figs WITHOUT extention
    
    try
        close(figure(1))
        close(figure(2))
    catch
    end
    
    idx=[];
    for i = 1:numel(cond(:,1))
        if length(varible{:,i}) < 10
            idx=[idx,i];
            continue
        end
        %plot PDF
        if sum(cell2mat(varible(:,1)) < 0)
            x = cell2mat(varible(:,i));
            binmax = ceil(max(max(x),-min(x))/binsize)*binsize;
            [n,x] = hist(x,-binmax:binsize:binmax);
        else
            [n,x] = hist(cell2mat(varible(:,i)),0:binsize:max(cell2mat(varible(:,i))));
        end
        
        n = n/length(cell2mat(varible(:,i)));
        figure(1);
        try 
            h = plot(x,n,'color',colorID(i,:));            
        catch
            h = plot(x,n);
        end
        set(h,'LineWidth',1.0)
        hold on
        
        %plot CDF
        figure(2);
        h = cdfplot(cell2mat(varible(:,i)));
        try 
            set(h,'color',colorID(i,:))             
        catch
        end
        set(h,'LineWidth',1.0)
        hold on
    end
    
    if ~isempty(idx)
        disp('Not enough data for:')
        disp(cond(idx))
        cond = cond(setdiff(1:length(cond),idx));
    end
    %save PDF
    figure(1)
    try
        legend(cond)
        xlabel(label)
    catch
    end
    ylabel('Probability')
    hold off
    savefig(figure(1),fullfile(outFolder,['PDF_',figname,'.fig']))
        
    %save CDF
    figure(2)
    try
        legend(cond,'Location','southeast')
        xlabel(label)
    catch
    end
    ylabel('Cumulative Probability')
    hold off
    savefig(figure(2),fullfile(outFolder,['CDF_',figname,'.fig']))
end

%%
function s = add_digits(n,nDigitsMax)
    nDigits = ceil(log(n+1)/log(10));
    nZerosToAdd = nDigitsMax - nDigits;
    s = [repmat('0',1,nZerosToAdd),num2str(n)];
end
function des1 = parse_description(des)
    idx = strfind(des,', ');
    des1 = string([]);
    for i = 1:length(idx)
        if i == 1
            des1 = [des1;des(1:idx(1)-1)];
        else
            des1 = [des1;des(idx(i-1)+2:idx(i)-1)];
            if i == length(idx)
                des1 = [des1;des(idx(i)+2:length(des))];
            end            
        end
    end
end
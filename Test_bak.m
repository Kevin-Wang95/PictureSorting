%% sort
%需将工作路径设置为当前project的路径
%输入的图片尺寸为640*480较为合适
%Preparation
curentPath  = pwd;
fileFolder=fullfile(pwd,'\unsorted_images');
dirOutput=dir(fullfile(fileFolder,'*.jpg'));
filenames={dirOutput.name}';
scale_test = imread(filenames{1});
scale = 0.125 * 480/size(scale_test,1);
Length = length(filenames);
src=cell(1,Length);
sorted=cell(1,Length);
% sample Preparation
samplepath = fullfile(pwd, 'Trainsample.jpg');
sample = imread(samplepath);
sample = imresize(sample,scale);
sample = sample(50:1:55, 1:1:50, 1:3);
%Data Preparation
waithand1=waitbar(0,'Loading Image');
for ii = 1:Length
    str=['Loading Image  ',num2str(ii),'/',num2str(Length)];
    waitbar(ii/Length,waithand1,str);
    src{ii} = imread(filenames{ii});
    sorted{ii} = imresize(src{ii},scale);
    sorted{ii} = sorted{ii}(50:1:55, 1:1:50, 1:3);
end
close(waithand1);
A1_cnt=0;A2_cnt=0;
test_rgb=zeros(3,1);
sample_rgb=zeros(3,1);
A1_list=cell(1,Length);
A2_list=cell(1,Length);
for i=1:3
    [r c]=size(sample(:,:,i));
    num=r*c;
    sample_rgb(i,1)=sum(sum(sample(:,:,i)))/num;
end
% sample = imresize(sample,scale);
fid1 = fopen('A111.txt','wt');
fid2 = fopen('A121.txt','wt');
PATH1 = 'E:\0ClassWork\信号与系统\大作业\1\';
PATH2 = 'E:\0ClassWork\信号与系统\大作业\2\';


waithand2=waitbar(0,'Processing Sort');
for i = 1:Length
    str=['Processing Sort  ',num2str(i),'/',num2str(Length)];
    waitbar(i/Length, waithand2, str);
%     err = image_sort(sample,sorted{i});
    for ii=1:3
    [r c]=size(sorted{i}(:,:,ii));
    num=r*c;
    test_rgb(ii,1)=sum(sum(sorted{i}(:,:,ii)))/num;
    end
    err = test_rgb - sample_rgb;
    [pathstr, name, ext] = fileparts(filenames{i});
    dis = sqrt(sum(err.^2));
    % The threshold is tested best when it is 0.2
    if(dis<2.5)
        fprintf(fid1,'%s\n',name);
        A1_cnt = A1_cnt+1;
        A1_list{A1_cnt}=filenames{i};

        OBJECT=fullfile(fileFolder,filenames{i});
        copyfile(OBJECT,PATH1);
    else
        fprintf(fid2,'%s\n',name);
        A2_cnt = A2_cnt+1;
        A2_list{A2_cnt}=filenames{i};
        OBJECT=fullfile(fileFolder,filenames{i});
        copyfile(OBJECT,PATH2);
    end
end
close(waithand2);
fclose(fid1);
fclose(fid2);

%% sort in time of speaker
% Paperation
% quene_start = 0;
% quene_end = 0;
% queue = cell(1,A1_cnt);
% surf = cell(1,A1_cnt);
% for i=1:A1_cnt
%     surf
% end




%% sort in time of earth
% The earth is at the area of 55~119,397~437


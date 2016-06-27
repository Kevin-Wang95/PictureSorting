clear all; 
%% sort
%需将工作路径设置为当前project的路径
%输入的图片尺寸为640*480较为合适
%Preparation
curentPath  = pwd;
fileFolder=fullfile(pwd,'\unsorted_images');
dirOutput=dir(fullfile(fileFolder,'*.jpg'));
filenames={dirOutput.name}';
Length = length(filenames);
src=cell(1,Length);
originimage=cell(1,Length);
sorted=cell(1,Length);
% sample Preparation
samplepath = fullfile(pwd, 'Trainsample.jpg');
sample = imread(samplepath);
sample = imresize(sample,0.125);
sample = sample(50:1:55, 1:1:50, 1:3);
%Data Preparation
waithand1=waitbar(0,'Loading Image');
for ii = 1:Length
    str=['Loading Image  ',num2str(ii),'/',num2str(Length)];
    waitbar(ii/Length,waithand1,str);
    src{ii} = imread(filenames{ii});
    originimage{ii}=src{ii};
    src{ii}=imresize(src{ii},0.5);
    sorted{ii} = imresize(src{ii},0.25);
    sorted{ii} = sorted{ii}(50:1:55, 1:1:50, 1:3);
end
close(waithand1);
A1_cnt=0;A2_cnt=0;
test_rgb=zeros(3,1);
sample_rgb=zeros(3,1);
for i=1:3
    [r c]=size(sample(:,:,i));
    num=r*c;
    sample_rgb(i,1)=sum(sum(sample(:,:,i)))/num;
end

% sample = imresize(sample,scale);
fid1 = fopen('A11.txt','wt');
fid2 = fopen('A12.txt','wt');
% PATH1 = 'E:\0ClassWork\信号与系统\大作业\1\';
% PATH2 = 'E:\0ClassWork\信号与系统\大作业\2\';

% Processing Sort
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
        A1_list{A1_cnt}=i;

%         OBJECT=fullfile(fileFolder,filenames{i});
%         copyfile(OBJECT,PATH1);
    else
        fprintf(fid2,'%s\n',name);
        A2_cnt = A2_cnt+1;
        A2_list{A2_cnt}=i;
%         OBJECT=fullfile(fileFolder,filenames{i});
%         copyfile(OBJECT,PATH2);
    end
end
close(waithand2);
fclose(fid1);
fclose(fid2);
%% sort in time of speakers
% Paperation
 queue_end = 0;
 queue=cell(1,0); 
 waithand3=waitbar(0,'Preparing Data');
 for i=1:A1_cnt
    str=['Preparing Data   ',num2str(i),'/',num2str(A1_cnt)];
    waitbar(i/A1_cnt, waithand3, str);
    tmp = im2double(rgb2gray(src{A1_list{i}}));
%     tmp = [tmp(50:109,1:110) tmp(50:109,160:320);tmp(130:193,1:110) tmp(130:193,160:320)];
%     tmp = [tmp(50:193,1:110) tmp(50:193,160:320)];
%     tmp = im2double(tmp);
    speakers(:,i)=reshape(tmp,size(tmp,1)*size(tmp,2),1);
 end
 close(waithand3);
 
 % Get the relationship
%  coef = corrcoef(speakers);
%  for i=1:size(coef)
%      coef(i,i)=-inf;
%  end
 D=pdist(speakers');
 l=1;
 for i=1:size(speakers,2)-1
     for j=i+1:size(speakers,2)
       coef(i,j)=D(l);
       coef(j,i)=coef(i,j);
       l=l+1;
     end
 end
%  coef = corrcoef(speakers);
 for i=1:size(coef)
     coef(i,i)=inf;
 end
% Start restoring two speaker
 waithand4=waitbar(0,'Restoring Speakers');
 rightpoint=45;
 leftpoint=45;
 % Run for the first point
 test_list=A1_list;
 queue_end=queue_end+1;
 queue{queue_end}=test_list{rightpoint};
 testcoef=coef;
 % Run for the second point
 [B1, I1] = sort(testcoef(:,rightpoint));
 queue_end=queue_end+1;
 queue{queue_end}=test_list{I1(1)};
 rightpoint=find(cell2mat(test_list)==queue{queue_end});
 while(size(test_list,2)~=2)
     str=['Restoring Speakers  ',num2str(A1_cnt-length(test_list)+2),'/',num2str(A1_cnt)];
     waitbar((A1_cnt-length(test_list))/A1_cnt, waithand4, str);
     [B1, I1] = sort(testcoef(:,rightpoint));
     [B2, I2] = sort(testcoef(:,leftpoint));
     if(I1(1)==leftpoint)
         B_1=B1(2);
         I_1=I1(2);
     else
         B_1=B1(1);
         I_1=I1(1);
     end
     if(I2(1)==rightpoint)
         B_2=B2(2);
         I_2=I2(2);
     else
         B_2=B2(1);
         I_2=I2(1);
     end
     if(B_1<B_2)
         queue_end=queue_end+1;
         queue{queue_end}=test_list{I_1};
         if(rightpoint<leftpoint)
             leftpoint=leftpoint-1;
         end
         test_list(rightpoint)=[];
         testcoef(rightpoint,:)=[];
         testcoef(:,rightpoint)=[];
         rightpoint=find(cell2mat(test_list)==queue{queue_end});
     else
         queue_end=queue_end+1;
         queue=[test_list{I_2} queue];
         if(leftpoint<rightpoint)
             rightpoint=rightpoint-1;
         end
         test_list(leftpoint)=[];
         testcoef(leftpoint,:)=[];
         testcoef(:,leftpoint)=[];
         leftpoint=find(cell2mat(test_list)==queue{1});
     end
 end
 close(waithand4);
 
 fid1 = fopen('A2.txt','wt'); 
 for i=length(queue):-1:1
    [pathstr, name, ext] = fileparts(filenames{queue{i}});
    imshow(originimage{queue{i}});
    fprintf(fid1,'%s\n',name);
 end
 close(figure(gcf));

%% sort the lens
% The earth is at the area of 55~119,397~437

% Data Preparetion
queue_speaker=queue;
clear queue;
waithand5=waitbar(0,'Data Preparetion For Sorting Lens');
for i=1:A2_cnt
    str=['Data Preparetion For Sorting Lens  ',num2str(i),'/',num2str(A2_cnt)];
    waitbar(i/A2_cnt, waithand5, str);
    tmp=imresize(src{A2_list{i}},0.4);
    tmp1=src{A2_list{i}}(197:224,22:64,:);
    scenes(:,i)=reshape(im2double(rgb2gray(tmp)),1,size(tmp,1)*size(tmp,2));
    globalscenes(:,i)=reshape(im2double(rgb2gray(tmp1)),1,size(tmp1,1)*size(tmp1,2));
end
close(waithand5);

% waithand6=waitbar(100,'Pca is Processing Data');
% [pcacoeff,score,latent]=pca(scenes');
% close(waithand6);

% for i=1:size(latent)
%     if(latent(i,1)>0.01)
%         cluster(:,i)=score(:,i);
%         clusterlatent(i,1)=latent(i,1);
%     end
% end
waithand6=waitbar(100,'KMeans is Processing Data');
[IDX, C] = kmeans(scenes',35);
close(waithand6);

PATH1 = 'E:\0ClassWork\信号与系统\大作业\ext\';
PATH2 = 'E:\0ClassWork\信号与系统\大作业\ext.bak1\';

% First Time cluster data preparation
waithand7=waitbar(0,'Cluster Data Orginazing');
cluster_length=zeros(1,35);
for i=1:A2_cnt
    str=['Cluster Data Orginazing  ',num2str(i),'/',num2str(A2_cnt)];
    waitbar(i/A2_cnt, waithand7, str);
    cluster_length(1,IDX(i,1))=cluster_length(1,IDX(i,1))+1;
%     OBJECT=fullfile(fileFolder,filenames{A2_list{i}});
%     str1=[PATH1,num2str(IDX(i,1)),'\'];
%     copyfile(OBJECT,str1);
    kmeans_sorted(cluster_length(1,IDX(i,1)),IDX(i,1))=A2_list{i};
    kmeans_sorted_image{IDX(i,1)}(:,cluster_length(1,IDX(i,1)))=scenes(:,i);
    kmeans_sorted_global{IDX(i,1)}(:,cluster_length(1,IDX(i,1)))=globalscenes(:,i);
end
close(waithand7);

% Second Time data preparation
cluster_num=35;clear queue_end queue;
queue_end=cell(1,0);
queue=cell(1,0);
i=1;
while(i<=cluster_num)
    if(cluster_length(i)>2)
        clear coef test_list test_image;
        test_list=zeros(cluster_length(1,i),1);
        test_list(:,1)=kmeans_sorted(1:cluster_length(1,i),i);
        tested_list=test_list;
        test_image=kmeans_sorted_image{i};
        tested_image=test_image;
        test_global=kmeans_sorted_global{i};
        tested_global=test_global;
        D=pdist(test_image');
        l=1;
        for ii=1:size(test_image,2)-1
            for j=ii+1:size(test_image,2)
                coef(ii,j)=D(l);
                coef(j,ii)=coef(ii,j);
                l=l+1;
            end
        end
        for ii=1:size(coef)
            coef(ii,ii)=inf;
        end
        coef2=corrcoef(test_global);
        for ii=1:size(coef2)
            coef2(ii,ii)=-inf;
        end
        queue_end{i}=0;
        queue{i}=cell(1,0);
        rightpoint=1;
        leftpoint=1;
        % Run for the first point
        queue_end{i}=queue_end{i}+1;
        queue{i}{queue_end{i}}=test_list(rightpoint,1);
        testcoef=coef;
        testglobal=coef2;
        % Run for the second point
        [B1, I1] = sort(testcoef(:,rightpoint));
        queue_end{i}=queue_end{i}+1;
        queue{i}{queue_end{i}}=test_list(I1(1),1);
        rightpoint=find(test_list==queue{i}{queue_end{i}});
        flag = true;
        while(size(test_list,1)~=2 && flag)
    %         str=['Restoring Speakers  ',num2str(A1_cnt-length(test_list)+2),'/',num2str(A1_cnt)];
    %         waitbar((A1_cnt-length(test_list))/A1_cnt, waithand7, str);
            [B1, I1] = sort(testcoef(:,rightpoint));
            [B2, I2] = sort(testcoef(:,leftpoint));
            if(I1(1)==leftpoint)
                B_1=B1(2);
                I_1=I1(2);
            else
                B_1=B1(1);
                I_1=I1(1);
            end
            if(I2(1)==rightpoint)
                B_2=B2(2);
                I_2=I2(2);
            else
                B_2=B2(1);
                I_2=I2(1);
            end
            if(B_1<B_2)
                if(testglobal(I_1,rightpoint)>0.920)
                    queue_end{i}=queue_end{i}+1;
                    queue{i}{queue_end{i}}=test_list(I_1,1);
                    if(rightpoint<leftpoint)
                        leftpoint=leftpoint-1;
                    end
                    test_list(rightpoint)=[];
                    testcoef(rightpoint,:)=[];
                    testcoef(:,rightpoint)=[];
                    testglobal(rightpoint,:)=[];
                    testglobal(:,rightpoint)=[];
                    rightpoint=find(test_list==queue{i}{queue_end{i}});
                else
                    flag = false;
                end
            else
                if(testglobal(I_2,leftpoint)>0.920)
                    queue_end{i}=queue_end{i}+1;
                    queue{i}=[test_list(I_2,1) queue{i}];
                    if(leftpoint<rightpoint)
                        rightpoint=rightpoint-1;
                    end
                    test_list(leftpoint)=[];
                    testcoef(leftpoint,:)=[];
                    testcoef(:,leftpoint)=[];
                    testglobal(leftpoint,:)=[];
                    testglobal(:,leftpoint)=[];
                    leftpoint=find(test_list==queue{i}{1});
                else
                    flag = false;
                end
            end
        end
        if(flag==false && size(test_list,1)~=2 )
            if(rightpoint<leftpoint)
                test_list(leftpoint)=[];
                test_list(rightpoint)=[];
            else
                test_list(rightpoint)=[];
                test_list(leftpoint)=[];
            end
            cluster_num=cluster_num+1;
            cluster_length(1,cluster_num)=size(test_list,1);
            cluster_length(1,i)=cluster_length(1,i)-size(test_list,1);
            kmeans_sorted(:,cluster_num)=zeros(size(kmeans_sorted,1),1);
            kmeans_sorted(1:size(test_list,1),cluster_num)=test_list;
            kmeans_sorted_image{cluster_num}(:,:)=zeros(size(test_image,1),cluster_length(1,cluster_num));
            kmeans_sorted_global{cluster_num}(:,:)=zeros(size(test_global,1),cluster_length(1,cluster_num));
            for k=1:size(test_list,1)
                pos=find(tested_list==test_list(k,1));
                tested_list(pos,:)=[];
                kmeans_sorted_image{cluster_num}(:,k)=tested_image(:,pos);
                kmeans_sorted_global{cluster_num}(:,k)=tested_global(:,pos);
                tested_image(:,pos)=[];
                tested_global(:,pos)=[];
            end
            kmeans_sorted_image{i}=tested_image;
            kmeans_sorted_global{i}=tested_global;
        end
    else
        queue{i}=kmeans_sorted(1:cluster_length(1,i),i);
    end
    i=i+1;
end

for iii=1:cluster_num
    for j=1:cluster_length(1,iii)
%       OBJECT=fullfile(fileFolder,filenames{queue{iii}{j}});
%       str2=[PATH2,num2str(iii),'\'];
%       copyfile(OBJECT,str2);
        if(cluster_length(1,iii)==1)
          imshow(originimage{queue{iii}});
        else
          imshow(originimage{queue{iii}{j}});
        end
    end
    close(figure(gcf));
end
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
 queue_start = 0;
 queue_end = 0;
 queue_left=cell(1,0);
 leftFlag=true;
 rightFlag=true;
 speakers = zeros(size(src{A1_list{1}},1)*size(src{A1_list{1}},2),A1_cnt);
 waithand3=waitbar(0,'Preparing Data');
 for i=1:A1_cnt
    str=['Preparing Data   ',num2str(i),'/',num2str(A1_cnt)];
    waitbar(i/A1_cnt, waithand3, str);
    tmp = im2double(rgb2gray(src{A1_list{i}}));
    speakers(:,i)=reshape(tmp,size(tmp,1)*size(tmp,2),1);
 end
 close(waithand3);
 
 % Get the relationship
 coef = corrcoef(speakers);
 for i=1:size(coef)
     coef(i,i)=0;
 end
 
% Start restoring two speaker
 waithand4=waitbar(0,'Restoring Speakers');
 startpoint=1;
 % Processing Right
 test_list=A1_list;
 while(rightFlag && size(test_list,2)~=0)
     testpoint=startpoint;
     queue_end=queue_end+1;
     queue_right{queue_end}=A1_list{startpoint};
     test_list=A1_list;
     testcoef=coef;
     startcoef=testcoef(:,startpoint);
     while(rightFlag && size(test_list,2)~=0)
         if(queue_end>19)
             str=['Restoring Speakers  ',num2str(A1_cnt-length(test_list)),'/',num2str(A1_cnt)];
             waitbar((A1_cnt-length(test_list))/A1_cnt, waithand4, str);
         end
         [B, I] = sort(testcoef(:,testpoint),'descend');
         if((1-B(1))<=0.007)
             queue_end=queue_end+1;
             queue_right{queue_end}=test_list{I(1)};
             test_list(testpoint)=[];
             testcoef(testpoint,:)=[];
             testcoef(:,testpoint)=[];
             startcoef(testpoint,:)=[];
             testpoint=find(cell2mat(test_list)==queue_right{queue_end});
         else
            if(queue_end<19)  
                queue_end=0;
                queue_right(:)=[];
                startpoint=startpoint+1;
                break;
            else
                test_list(testpoint)=[];
                testcoef(testpoint,:)=[];
                testcoef(:,testpoint)=[];
                startcoef(testpoint,:)=[];
                rightFlag=false;
            end
         end
     end
 end
 % Before Processing Left 
 if(size(test_list,2)~=0)
    [B, I] = sort(startcoef,'descend');
    if((1-B(1))<=0.007)
        queue_start=queue_start+1;
        queue_left{queue_start}=test_list{I(1)};
        testpoint=find(cell2mat(test_list)==queue_left{queue_start});
    else
        leftFlag=false;
    end
 else
    leftFlag =false; 
 end
 % Processing Left
 while(size(test_list,2)~=0 && leftFlag)
       str=['Restoring Speakers  ',num2str(A1_cnt-length(test_list)),'/',num2str(A1_cnt)];
       waitbar((A1_cnt-length(test_list))/A1_cnt, waithand4, str);
       
       [B, I] = sort(testcoef(:,testpoint),'descend');
       if((1-B(1))<=0.007)
           queue_start=queue_start+1;
           queue_left{queue_start}=test_list{I(1)};
           test_list(testpoint)=[];
           testcoef(testpoint,:)=[];
           testcoef(:,testpoint)=[];
           testpoint=find(cell2mat(test_list)==queue_left{queue_start});
       else
           leftFlag=false;
           test_list(testpoint)=[];
           testcoef(testpoint,:)=[];
           testcoef(:,testpoint)=[];
       end
 end
 
 % Process remained left
 remain_right = true;
 l=1;
 while(size(test_list,2)~=0 && remain_right)
      str=['Restoring Speakers  ',num2str(A1_cnt-length(test_list)),'/',num2str(A1_cnt)];
      waitbar((A1_cnt-length(test_list))/A1_cnt, waithand4, str);
      
      max=0;
      for i=1:length(queue_right)-1
          tmp0 = im2double(rgb2gray(src{test_list{l}}));
          tmp1 = im2double(rgb2gray(src{queue_right{i}}));
          tmp2 = im2double(rgb2gray(src{queue_right{i+1}}));
          tmp0=reshape(tmp0,size(tmp0,1)*size(tmp0,2),1);
          tmp1=reshape(tmp1,size(tmp1,1)*size(tmp1,2),1);
          tmp2=reshape(tmp2,size(tmp2,1)*size(tmp2,2),1);
          tmpcoef1=corrcoef([tmp1 tmp0]);
          tmpcoef2=corrcoef([tmp2 tmp0]);
          tmpcoef3=corrcoef([tmp1 tmp2]);
          temp=tmpcoef1(2,1)+tmpcoef2(2,1)-tmpcoef3(2,1);
          if(max<temp)
              max=temp;
              pos=i;
          end
      end
      tmp1 = im2double(rgb2gray(src{queue_right{pos}}));
      tmp2 = im2double(rgb2gray(src{queue_right{pos+1}}));
      tmp1=reshape(tmp1,size(tmp1,1)*size(tmp1,2),1);
      tmp2=reshape(tmp2,size(tmp2,1)*size(tmp2,2),1);
      tmpcoef1=corrcoef([tmp1 tmp0]);
      tmpcoef2=corrcoef([tmp2 tmp0]);
      if(tmpcoef1(2,1)<0.99||tmpcoef2(2,1)<0.99)
          l=l+1;
          if(l>size(test_list,2))
              remain_right=false;
          end
      else
          queue_right=[queue_right(1:pos) test_list{l} queue_right(pos+1:queue_end)];
          queue_end=queue_end+1;
          test_list(l)=[];
      end
 end
 
 % Process remain left
 while(size(test_list,2)~=0)
      str=['Restoring Speakers  ',num2str(A1_cnt-length(test_list)),'/',num2str(A1_cnt)];
      waitbar((A1_cnt-length(test_list))/A1_cnt, waithand4, str);
      
      max=0;
      for i=1:length(queue_left)-1
          tmp0 = im2double(rgb2gray(src{test_list{1}}));
          tmp1 = im2double(rgb2gray(src{queue_left{i}}));
          tmp2 = im2double(rgb2gray(src{queue_left{i+1}}));
          tmp0=reshape(tmp0,size(tmp0,1)*size(tmp0,2),1);
          tmp1=reshape(tmp1,size(tmp1,1)*size(tmp1,2),1);
          tmp2=reshape(tmp2,size(tmp2,1)*size(tmp2,2),1);
          tmpcoef1=corrcoef([tmp1 tmp0]);
          tmpcoef2=corrcoef([tmp2 tmp0]);
          tmpcoef3=corrcoef([tmp1 tmp2]);
          temp=tmpcoef1(2,1)+tmpcoef2(2,1)-tmpcoef3(2,1);
          if(max<temp)
              max=temp;
              pos=i;
          end
      end
      
       queue_left=[queue_left(1:pos) test_list{1} queue_left(pos+1:queue_start)];
       queue_start=queue_start+1;
       test_list(1)=[];
 end
  close(waithand4);
  
 if(size(queue_left,2)~=0)
    queue = [queue_left(length(queue_left):-1:1) queue_right(1:length(queue_right))];  
 else
    queue = queue_right;
 end
 fid1 = fopen('A2.txt','wt'); 
 for i=1:length(queue)
    [pathstr, name, ext] = fileparts(filenames{queue{i}});
     fprintf(fid1,'%s\n',name);
 end

%% sort the lens
% The earth is at the area of 55~119,397~437

% Data Preparetion
waithand5=waitbar(0,'Data Preparetion');
for i=1:A2_cnt
    str=['Data Preparetion  ',num2str(i),'/',num2str(A2_cnt)];
    waitbar(i/A2_cnt, waithand5, str);
    tmp=imresize(src{A2_list{i}},0.4);
    scenes(:,i)=reshape(im2double(rgb2gray(tmp)),1,size(tmp,1)*size(tmp,2));
end
close(waithand5);
[pcacoeff,score,latent]=pca(scenes');

for i=1:size(latent)
    if(latent(i,1)>0.02)
        cluster(:,i)=score(:,i);
        clusterlatent(i,1)=latent(i,1);
    end
end

[IDX, C] = kmeans(cluster,35);
PATH = 'E:\0ClassWork\信号与系统\大作业\ext\';
 
waithand6=waitbar(0,'Cluster Data Orginazing');
cluster_length=zero(1,35);
for i=1:A2_cnt
%       OBJECT=fullfile(fileFolder,filenames{A2_list{i}});
%       str=[PATH,num2str(IDX(i,1)),'\']
%       copyfile(OBJECT,str);
    str=['Cluster Data Orginazing  ',num2str(i),'/',num2str(A2_cnt)];
    waitbar(i/A2_cnt, waithand6, str);
    cluster_length(1,IDX(i,1))=cluster_length(1,IDX(i,1))+1;
    kmeans_sorted(cluster_length(1,IDX(i,1)),IDX(i,1))=A2_list{i};
    kmeans_sorted_image{IDX(i,1)}(:,cluster_length(1,IDX(i,1)))=scenes(:,i);
end
close(waithand6);


for i=1:35
    tested_list(:,:)=[];tested_image(:,:)=[];
    tested_list(:,1)=kmeans_sorted(1:cluster_length(1,IDX(i,1)),i);
    tested_image=kmeans_sorted_image{i};
    
    
end

%% sort
%需将工作路径设置为当前project的路径
%输入的图片尺寸为640*480较为合适
%Preparation
curentPath  = pwd;
fileFolder=fullfile(pwd,'\Inorder2');
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
fid1 = fopen('A111.txt','wt');
fid2 = fopen('A121.txt','wt');
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

%% sort in time of speaker
% Paperation
 queue_start = 0;
 queue_end = 0;
 queue_left=cell(1,0);
 leftFlag=true;
 rightFlag=true;
 
% start restoring two speaker
 waithand5=waitbar(0,'Restoring Speakers');
 coef=zeros(1,A1_cnt);
 startpoint=1;
 while(rightFlag)
     testpoint=startpoint;
     queue_end=queue_end+1;
     queue_right{queue_end}=A1_list{startpoint};
     test_list=A1_list;
     test_vector=src{A1_list{startpoint}};
     start_vector=test_vector;
     while(rightFlag)
         if(queue_end>12)
             str=['Restoring Speakers  ',num2str(A1_cnt-length(test_list)),'/',num2str(A1_cnt)];
             waitbar((A1_cnt-length(test_list))/A1_cnt, waithand5, str);
         end
         coef=zeros(1,A1_cnt);
         test_list(testpoint)=[];
         for i=1:length(test_list)
            coef(1,i)=Feature(test_vector,src{test_list{i}});
         end
         [B, I] = sort(coef,'descend');
         if((1-B(1))<=0.01)
             testpoint=I(1);
             test_vector=src{test_list{testpoint}};
             queue_end=queue_end+1;
             queue_right{queue_end}=test_list{I(1)};
         else
            if(queue_end<12)  
                queue_end=0;
                queue_right(:)=[];
                startpoint=startpoint+1;
                break;
            else
                rightFlag=false;
            end
         end
     end
 end
 
  coef=zeros(1,A1_cnt);
  for i=1:length(test_list)
      coef(1,i)=Feature(start_vector,src{test_list{i}});
  end
       [B, I] = sort(coef,'descend');
       if((1-B(1))<=0.01)
           testpoint=I(1);
           test_vector=src{test_list{testpoint}};
           queue_start=queue_start+1;
           queue_left{queue_start}=test_list{I(1)};
       else
           leftFlag=false;
       end
  
 while(size(test_list,2)~=1 && leftFlag)
       str=['Restoring Speakers  ',num2str(A1_cnt-length(test_list)),'/',num2str(A1_cnt)];
       waitbar((A1_cnt-length(test_list))/A1_cnt, waithand5, str);
       coef=zeros(1,A1_cnt);
       test_list(testpoint)=[];
       for i=1:length(test_list)
          coef(1,i)=Feature(test_vector,src{test_list{i}});
       end
       [B, I] = sort(coef,'descend');
       if((1-B(1))<=0.01)
          testpoint=I(1);
          test_vector=src{test_list{testpoint}};
          queue_start=queue_start+1;
          queue_left{queue_start}=test_list{I(1)};
       else
           leftFlag=false;
       end
 end
 
 while(size(test_list,2)~=0)
      str=['Restoring Speakers  ',num2str(A1_cnt-length(test_list)),'/',num2str(A1_cnt)];
      waitbar((A1_cnt-length(test_list))/A1_cnt, waithand5, str);
      coef=zeros(1,A1_cnt);
      test_vector=src{test_list{1}};
      for i=1:length(queue_right)
        coef(1,i)=Feature(test_vector,src{queue_right{i}});
      end
      max=0;
      for i=1:length(queue_right)-1
          temp=coef(1,i)+coef(1,i+1)-Feature(src{queue_right{i+1}},src{queue_right{i}});
          if(max<temp)
              max=temp;
              pos=i;
          end
      end
      queue_right=[queue_right(1:pos) test_list{1} queue_right(pos+1:queue_end)];
      queue_end=queue_end+1;
      test_list(1)=[];
 end
  close(waithand5);
  
 if(size(queue_left,2)~=0)
    queue = [queue_left(length(queue_left):-1:1) queue_right(1:length(queue_right))];  
 else
    queue = queue_right;
 end
 fid1 = fopen('A21.txt','wt'); 
 for i=1:length(queue)
    [pathstr, name, ext] = fileparts(filenames{queue{i}});
     fprintf(fid1,'%s\n',name);
 end

%% sort the lens
% The earth is at the area of 55~119,397~437

% Data Preparetion
coef = zeros(A2_cnt, A2_cnt);
for i=1:A2_cnt
   for j=1:(i-1)
       coef(i,j)=Feature(src{i},src{j});
   end
end

clc
clear all
for i=1:68
%%%%%%%%i for i-th file
 if i<10
     fileindex=strcat('0',int2str(i));
% elseif i>9 && i<100
%     fileindex=strcat('0',int2str(i));
 end
%fileindex=int2str(i);
filename=strcat('point-pressure',fileindex);
filename=strcat(filename,'.out');
%%%%%%%%%store file number into fileindex
if exist(filename, 'file')==0
    continue
end
%%%%%%%%%check file exist or not
filec=importdata(filename);
data=filec.data;
%%%%%%%%%load i-th file
if i==1
NUM=data(:,2);
NUMcap=[0;NUM];
xlswrite('pressure.xlsx',NUMcap,'sheet1','A');
end
%%%%%%%%%initialize the first column of xls file 
datapressure=data(:,3);
datapressurecap=[i;datapressure];
if i<26
    COL=char('A'+i);
elseif i>25 && i<676
    COL=strcat(char('A'-1+i/26),char('A'+rem(i,26)));
end
%%%%%%%%%check the index of columns
xlswrite('pressure.xlsx',datapressurecap,'sheet1',COL);
%%%%%%%%%output to corresponding columns
end
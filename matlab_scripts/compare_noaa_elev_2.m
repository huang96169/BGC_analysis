%comparison of model output and obs for the 63 coastal stations
clear all
%-----------------------inputs---------------------------
%start_datenum=datenum(2017,8,4);%Harvey
start_datenum=datenum(2018,8,1);%Florence start time of model outputs
%start_datenum=datenum(2011,7,27);
time_1=datenum(2018,8,1)+30; %start time of user's time period
time_2=datenum(2018,8,1)+40; %end time of user's time period
ndays=40;%days of model outputs

runo='RUN10k'; 
station_file_nameo='Pacific_2sta_o';%'upbay_sta53.moved2';%'upbay_sta53';%coast63.moved6 

run='RUN10k'; 
station_file_name='Pacific_2sta_mv1';%'upbay_sta53.moved2';%'upbay_sta53';%coast63.moved5 

%dirs
%prj_dir='/sciclone/home20/whuang07/schism10/NWM/Case1/NOAATIDES/Forecast/results/';  %one level up RUN*
prj_dir='/sciclone/home20/whuang07/schism20/Pacific/Case1/';

%noaa_obs_dir='/sciclone/home20/whuang07/schism10/NWM/Case1/NOAATIDES/Harvey2017/';
%noaa_obs_dir='~feiye/schism10/work/NOAA_TIDE/';
%noaa_obs_dir='/sciclone/home20/whuang07/schism10/NWM/Case1/NOAATIDES/Florence/';
noaa_obs_dir='/sciclone/home20/whuang07/schism20/Pacific/NOAAobs/';
%noaa_obs_dir='/sciclone/home20/whuang07/schism10/NWM/Case1/NOAATIDES/Forecast/';
%noaa_obs_dir='/sciclone/home20/whuang07/schism10/NWM/Case1/NOAATIDES/NOAA_TIDE_Irene/';
station_dir='/sciclone/home20/whuang07/git/NWM_scripts/matlab_scripts/Elev/BPfiles/';
saved_mat_dir='/sciclone/schism10/whuang07/NWM/Case1/SAVED_MAT/'; %
saved_mat_dir='/sciclone/home20/whuang07/schism20/Pacific/Case1/SAVED_MAT/';


%---------------------end inputs-------------------------


outname=['elev.' station_file_name '.' run];
outnameo=['elev.' station_file_nameo '.' runo];

% station id and name  
f1=fopen([station_dir '/stations.txt']);
[tmp]=textscan(f1,'%s%s','delimiter',',');
stIds=tmp{1,1};
stNames=tmp{1,2};
fclose(f1);

% read lat/lon for stations 
fid=fopen([station_dir '/' station_file_name '.bp']);
[tmp]=textscan(fid,'%d',1,'headerlines',1); nf = double(tmp{1});
[tmp]=textscan(fid,'%d%f%f%f%d');
sa_lon=tmp{1,2};
sa_lat=tmp{1,3};
sa_id=tmp{1,5};

%load SCHISM output
output=load([prj_dir '/'  run '/PostP/' outname]);
dt = 1/24; %round to integer hours, then /24
%output=output(1:ndays*(dt*24*24)-1,:);
time_model=start_datenum+1/24:dt:start_datenum+dt*(length(output(:,1))-1); 
idx=find(time_model>=time_1&time_model<=time_2);
time_model=time_model(idx);
elem2=output(idx,2:end);
tt=1/24:dt:size(elem2,1)/24;

%load SCHISM output to be compared
output=load([prj_dir '/' runo '/PostP/' outnameo]);
elem1=output(idx,2:end);

eleo=cell(1,nf);
tt1=cell(1,nf);

ff=fopen([prj_dir '/NOAAtidecompare.' outname '.' outnameo '.log'],'w');
fprintf(ff,'%s\n',[outname outnameo]);
VD=cell(1,nf);
staname=cell(1,nf);
for i=1:nf
    
    id2=find(str2double(stIds)==sa_id(i));
    fname1=[noaa_obs_dir '/NAVD/' stNames{id2}];
    fname2=[noaa_obs_dir '/MSL/' stNames{id2}];
    if (exist(fname1,'dir')~=0)
       fprintf(ff,'%d %d %s\n',[i sa_id(i) 'datum: NAVD']);
       VD{1,i}=cellstr('NAVD');
       fname=[fname1 '/' stNames{id2} '.csv'];
       staname{1,i}=stNames{id2};
       fileID=fopen(fname);
       C=textscan(fileID,'%s %f %d %d %d %d %d %s','Delimiter',','); 
       tmp=[cell2mat(C{1,1}) repmat(':00',size(cell2mat(C{1,1}),1),1)];
       tmp2=cellstr(string(tmp));
       tmp3=DateStr2Num(cellstr(tmp2),31);
       ids=find(abs(tmp3-time_model(1))<=6/60/24);
       %ids=ids(1);
       ide=find(abs(tmp3-time_model(end))<=10^-3);
       
       if(isempty(ids)==0&&isempty(ide)==0)
         eleo{1,i}=interp1(tmp3(ids:ide),C{1,2}(ids:ide),time_model);
         tt1{1,i}=C{1,1}(ids:ide);
       else
         eleo{1,i}=[];
         tt1{1,i}=[];
       end
       fclose(fileID);
    %end
    elseif (exist(fname2,'dir')~=0)
       fprintf(ff,'%d %d %s\n',[i sa_id(i) 'datum: MSL']);
       VD{1,i}=cellstr('MSL');
       fname=[fname2 '/' stNames{id2} '.csv'];
       fileID=fopen(fname);
       C=textscan(fileID,'%s %f %d %d %d %d %d %s','Delimiter',','); 
       tmp=[cell2mat(C{1,1}) repmat(':00',size(cell2mat(C{1,1}),1),1)];
       tmp2=cellstr(string(tmp));
       tmp3=DateStr2Num(cellstr(tmp2),31);
       ids=find(abs(tmp3-time_model(1))<=10^-3);
       ide=find(abs(tmp3-time_model(end))<=10^-3);

       if(isempty(ids)==0&&isempty(ide)==0)
         eleo{1,i}=interp1(tmp3(ids:ide),C{1,2}(ids:ide),time_model);
         tt1{1,i}=C{1,1}(ids:ide);
       else
         eleo{1,i}=[];
         tt1{1,i}=[];
       end
       fclose(fileID);
    end

end

n_nemp=0;%number of stations has observed values
for i=1:nf
    if (isempty(eleo{1,i})==0)      
        n_nemp=n_nemp+1;
    end
end

tmp=intersect(station_file_name,station_file_nameo,'stable');
fName_savedObs=[saved_mat_dir '/' tmp '_' datestr(time_model(1),'yyyymmdd') '-' datestr(time_model(end),'yyyymmdd') '.mat']; %###not a good name, use bpfilename
if ~exist (fName_savedObs,'file')
    save(fName_savedObs,'eleo');
else
    load(fName_savedObs);
end

figure  
dn=1;
cc2=zeros(1,nf);
rmse2=zeros(1,nf);
mae2=zeros(1,nf);
SS2=zeros(1,nf);
cc1=zeros(1,nf);
rmse1=zeros(1,nf);
mae1=zeros(1,nf);
SS1=zeros(1,nf);
ncol=ceil(nf/10);

n=0;
for i=1:nf
     
       
       i
       h(i)=subplot(10,ncol,i);    
       %plot(tt,(detrend(eleos{1,i})),'r');hold on;
       %plot(tt,(detrend(elem1(:,i))),'g'); 
       %plot(tt,(detrend(elem2(:,i))),'b'); 
       tmpo=eleo{1,i}; 
     if(isempty(eleo{1,i})==0&&isempty(tt1{1,i})==0)
       id=find(isnan(tmpo)==0);
       aa=corrcoef(tmpo(id),elem2(id,i));
       cc=round(aa(1,2),2);
       mae=sum(abs(detrend(tmpo(id))-detrend(elem2(id,i)')))/length(time_model(id));
       %mae=sum(abs(tmpo(id)-elem2(id,i)'))/length(time_model(id));
       mae2(i)=round(mae,2);
       plot(tt,detrend(tmpo),'r');hold on;
       plot(tt,detrend((elem1(:,i))),'g');hold on
       plot(tt,detrend((elem2(:,i))),'b');
       n=n+1;
       if (n==1)
          legend({'obs',outnameo,outname});
       end
       
     else
       mae2(i)=NaN;
       cc=NaN;
       plot(tt,((elem1(:,i))),'g');hold on
       plot(tt,((elem2(:,i))),'b');
     end
       title({['CC=',num2str(cc),'MAE=',num2str(mae2(i))]
             [' ID=',num2str(sa_id(i)),staname{1,i},' Datum:',cell2mat(VD{1,i})]})
       xlim([1 10]);
       %ylim([-1 1]);
     
end

fclose(ff);

%set(h(1),'position',[0.01,0.9,0.2,0.1])

name_fig=strcat(prj_dir,'/Figures','/elev_',outname,'_',outnameo,'.fig');  
savefig(name_fig)
% save the new selected stations into new *.bp file and xyz file.

%figure
%plot(4:7,cc2(4:7),'bo');hold on;plot(4:7,cc1(4:7),'go');
%figure
%plot(4:7,mae2(4:7),'bo');hold on;plot(4:7,mae1(4:7),'go');
%figure
%plot(4:7,rmse2(4:7),'bo');hold on;plot(4:7,rmse1(4:7),'go');
%figure
%plot(4:7,SS2(4:7),'bo');hold on;plot(4:7,SS1(4:7),'go')


%figure
%geoplot(ss_lat,ss_lon,'ro')
%hold on
%for i=1:n_nemp
%    text(ss_lat(i),ss_lon(i),num2str(i),'FontSize',12);
%end

%ff=fopen(strcat('./stations/coast_s',num2str(n_nemp),'_origin.xyz','w'));
%for i=1:n_nemp
%    fprintf(ff,'%d %f %f %d\n',[i ss_lon(i) ss_lat(i) ss_id(i)]);
%end
%fclose(ff);

clear all;

%--------------inputs-----------------
%Irene
%start_datenum=datenum(2011,7,27);
%ndays=30;

%Harvey
%start_datenum=datenum(2017,8,4);
%ndays=70;

%Florence
start_datenum1=datenum(2018,6,17);
start_datenum2=datenum(2021,4,1);
ndays=20;
time_obs1=start_datenum1:1/24:start_datenum1+20;
time_obs2=start_datenum2:1/24:start_datenum2+20;

prj_dir='/sciclone/home20/whuang07/schism10/NWM/Case1/';  %one level up RUN*

%noaa_obs_dir='/sciclone/home20/whuang07/schism10/NWM/Case1/NOAATIDES/Harvey2017/';
%noaa_obs_dir='/sciclone/home20/whuang07/schism10/NWM/Case1/NOAATIDES/NOAA_TIDE_Irene/';
%noaa_obs_dir='/sciclone/schism10/feiye/work/NOAA_TIDE/20180823-20181001/';
noaa_obs_dir1='/sciclone/home20/whuang07/schism10/NWM/Case1/NOAATIDES/Florence/';
noaa_obs_dir2='/sciclone/home20/whuang07/schism10/NWM/Case1/NOAATIDES/Forecast/';
station_dir='/sciclone/home20/whuang07/git/NWM_scripts/matlab_scripts/Elev/BPfiles/';

station_file_name='8454000';
%-------------end inputs--------------

t_tide_outdir=[prj_dir '/T_Tide_out/'];

addpath ./T_Tides;
if (~exist(t_tide_outdir,'dir')) 
  system(['mkdir ' t_tide_outdir]);
end

% station id and name  
f1=fopen([station_dir '/stations.txt']);
[tmp]=textscan(f1,'%s%s','delimiter',',');
stIds=tmp{1,1};
stNames=tmp{1,2};
fclose(f1);


fid=fopen([station_dir '/' station_file_name '.bp']);
[tmp]=textscan(fid,'%d',1,'headerlines',1); nf = double(tmp{1});
[tmp]=textscan(fid,'%d%f%f%f%d');
sa_lon=tmp{1,2};
sa_lat=tmp{1,3};
sa_id=tmp{1,5};

%load SCHISM output

time_ha=1/24:1/24:ndays; %for T_TIDE


start_y1=str2num(datestr(time_obs1(1),'yyyy'));
start_m1=str2num(datestr(time_obs1(1),'mm'));
start_d1=str2num(datestr(time_obs1(1),'dd'));

start_y2=str2num(datestr(time_obs2(1),'yyyy'));
start_m2=str2num(datestr(time_obs2(1),'mm'));
start_d2=str2num(datestr(time_obs2(1),'dd'));

eleo=cell(1,nf);
tt1=cell(1,nf);
for i=1:nf
    i
    id2=find(str2double(stIds)==sa_id(i));
    fname1=[noaa_obs_dir1 'NAVD/' stNames{id2}];
    fname2=[noaa_obs_dir1 'MSL/' stNames{id2}];
    if (exist(fname1,'dir')~=0)
       fname=[fname1 '/' stNames{id2} '.csv'];
       fileID=fopen(fname);
       C=textscan(fileID,'%s %f %d %d %d %d %d %s','Delimiter',','); 
       tmp=[cell2mat(C{1,1}) repmat(':00',size(cell2mat(C{1,1}),1),1)];
       tmp2=cellstr(string(tmp));
       tmp3=DateStr2Num(cellstr(tmp2),31);
       tin=find(tmp3>=time_obs1(1) & tmp3<=time_obs1(end));
       if(isempty(tin)==0)
         eleo{1,i}=interp1(tmp3(tin),C{1,2}(tin),time_obs1);
         tt1{1,i}=C{1,1}(tin);
       else
         eleo{1,i}=[];
         tt1{1,i}=[];
       end
       fclose(fileID);
    
    elseif (exist(fname2,'dir')~=0)
       fname=[fname2 '/' stNames{id2} '.csv'];
       fileID=fopen(fname);
       C=textscan(fileID,'%s %f %d %d %d %d %d %s','Delimiter',','); 
       tmp=[cell2mat(C{1,1}) repmat(':00',size(cell2mat(C{1,1}),1),1)];
       tmp2=cellstr(string(tmp));
       tmp3=DateStr2Num(cellstr(tmp2),31);
       tin=find(tmp3>=time_obs1(1) & tmp3<=time_obs1(end));
       if(isempty(tin)==0&&isempty(find(isnan(C{1,2}(tin))))==1)
         eleo{1,i}=interp1(tmp3(tin),C{1,2}(tin),time_obs1);
         tt1{1,i}=C{1,1}(tin);
       else
         eleo{1,i}=[];
         tt1{1,i}=[];
       end
       fclose(fileID);
    end
end
n_nemp=nf;
eleos1=eleo;
for i=1:nf
    i
    id2=find(str2double(stIds)==sa_id(i));
    fname1=[noaa_obs_dir2 'NAVD/' stNames{id2}];
    fname2=[noaa_obs_dir2 'MSL/' stNames{id2}];
    if (exist(fname1,'dir')~=0)
       fname=[fname1 '/' stNames{id2} '.csv'];
       fileID=fopen(fname);
       C=textscan(fileID,'%s %f %d %d %d %d %d %s','Delimiter',',');
       tmp=[cell2mat(C{1,1}) repmat(':00',size(cell2mat(C{1,1}),1),1)];
       tmp2=cellstr(string(tmp));
       tmp3=DateStr2Num(cellstr(tmp2),31);
       tin=find(tmp3>=time_obs2(1) & tmp3<=time_obs2(end));
       if(isempty(tin)==0)
         eleo{1,i}=interp1(tmp3(tin),C{1,2}(tin),time_obs2);
         tt1{1,i}=C{1,1}(tin);
       else
         eleo{1,i}=[];
         tt1{1,i}=[];
       end
       fclose(fileID);

    elseif (exist(fname2,'dir')~=0)
       fname=[fname2 '/' stNames{id2} '.csv'];
       fileID=fopen(fname);
       C=textscan(fileID,'%s %f %d %d %d %d %d %s','Delimiter',',');
       tmp=[cell2mat(C{1,1}) repmat(':00',size(cell2mat(C{1,1}),1),1)];
       tmp2=cellstr(string(tmp));
       tmp3=DateStr2Num(cellstr(tmp2),31);
       tin=find(tmp3>=time_obs2(1) & tmp3<=time_obs2(end));
       if(isempty(tin)==0&&isempty(find(isnan(C{1,2}(tin))))==1)
         eleo{1,i}=interp1(tmp3(tin),C{1,2}(tin),time_obs2);
         tt1{1,i}=C{1,1}(tin);
       else
         eleo{1,i}=[];
         tt1{1,i}=[];
       end
       fclose(fileID);
    end
end
eleos2=eleo;

%list_con={'O1','K1','Q1','P1','K2','N2','M2','S2'}; %list of consti. to be extracted by HA
%list_con={'O1','K1','Q1','K2','N2','M2','S2'}; %list of consti. to be extracted by HA
list_con={'O1','K1','Q1','N2','M2','S2'};

%system(['rm ' t_tide_outdir '/*']); 
for i=1:n_nemp
  obs1=eleos1{1,i};
  obs2=eleos2{1,i};
%Harmonic analysis
  if (~isempty(obs1))
    [cons,xout]=t_tide(obs1,'interval',1.,'start time',[start_y1,start_m1,start_d1,0,0,0],'output',strcat(t_tide_outdir,'/obs.H.',num2str(sa_id(i))));
    [cons,xout]=t_tide(obs2,'interval',1.,'start time',[start_y2,start_m2,start_d2,0,0,0],'output',strcat(t_tide_outdir,'/obs.F.',num2str(sa_id(i))));
  end
%low-pass filter
    %obs_lo=simple_filter_data([time_ha;obs]',1./3,'low'); %3 day cut-off
    %mod_lo=simple_filter_data([time_ha;m2]',1./3,'low');
end

fname=[t_tide_outdir '/obs.H.' num2str(sa_id(i))];
fileID=fopen(fname);
obs1=textscan(fileID,'%s %f %f %f %f %f %f','headerlines',16);
    %com=obs{1,1};
amp_o1{1,i}=obs1{1,3};
pha_o1{1,i}=obs1{1,5};
fclose(fileID);
cm='M2';
tc=length(amp_o1{1,i});
for nc=1:tc
    %mod{1,1}(nc)
    if(contains(obs1{1,1}(nc),cm))
      am1=nc
    end
end
fname=[t_tide_outdir '/obs.F.' num2str(sa_id(i))];
fileID=fopen(fname);
obs2=textscan(fileID,'%s %f %f %f %f %f %f','headerlines',16);
    %com=obs{1,1};
amp_o2{1,i}=obs2{1,3};
pha_o2{1,i}=obs2{1,5};
fclose(fileID);

tc=length(amp_o2{1,i});
for nc=1:tc
    %mod{1,1}(nc)
    if(contains(obs2{1,1}(nc),cm))
      am2=nc
    end
end

amp_o1=amp_o1{1,i}(am1);
amp_o2=amp_o2{1,i}(am2);
pha_o1=pha_o1{1,i}(am1);
pha_o2=pha_o2{1,i}(am2);
figure
subplot(2,1,1)
plot(time_obs1,eleos1{1,1},'b')
datetick('x')
subplot(2,1,2)
plot(time_obs2,eleos2{1,1},'b')
datetick('x')
saveas(gcf,strcat(prj_dir,'/Figures/elev.',station_file_name,cm,'.png'))
figure
subplot(2,1,1)
scatter(1,amp_o1);hold on;scatter(1,amp_o2);
title([cm,'amp'])
subplot(2,1,2)
scatter(1,pha_o1);hold on;scatter(1,pha_o2);
title([cm,'pha'])
saveas(gcf,strcat(prj_dir,'/Figures/',station_file_name,cm,'.png'))

clear all;

%--------------inputs-----------------
%Irene
start_datenum=datenum(2011,7,27);
ndays=50;

%Harvey
%start_datenum=datenum(2017,8,4);
%ndays=70;

%Florence
%start_datenum=datenum(2018,8,24);
%ndays=36;

prj_dir='/sciclone/home20/whuang07/schism10/NWM/Case1/';  %one level up RUN*

%noaa_obs_dir='/sciclone/home20/whuang07/schism10/NWM/Case1/NOAATIDES/Harvey2017/';
noaa_obs_dir='/sciclone/home20/whuang07/schism10/NWM/Case1/NOAATIDES/NOAA_TIDE_Irene/';
%noaa_obs_dir='/sciclone/schism10/feiye/work/NOAA_TIDE/20180823-20181001/';
station_dir='/sciclone/home20/whuang07/git/NWM_scripts/matlab_scripts/Elev/BPfiles/';

run='RUN16a'; 
station_file_name='7station_harvey'; 
%-------------end inputs--------------

outname=['elev.' station_file_name '.' run];
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

start_y=str2num(datestr(start_datenum,'yyyy'));
start_m=str2num(datestr(start_datenum,'mm'));
start_d=str2num(datestr(start_datenum,'dd'));

fid=fopen([station_dir '/' station_file_name '.bp']);
[tmp]=textscan(fid,'%d',1,'headerlines',1); nf = double(tmp{1});
[tmp]=textscan(fid,'%d%f%f%f%d');
sa_lon=tmp{1,2};
sa_lat=tmp{1,3};
sa_id=tmp{1,5};

%load SCHISM output
output=load([prj_dir '/'  run '/PostP/' outname]);
time=start_datenum:1/24:start_datenum+1/24*(length(output(:,1))-1);
elem=output(:,2:end);
tt=output(:,1);

time_ha=1/24:1/24:ndays; %for T_TIDE

eleo=cell(1,nf);
tt1=cell(1,nf);

for i=1:nf
    id2=find(str2double(stIds)==sa_id(i));
    fname1=[noaa_obs_dir 'NAVD/' stNames{id2}];
    fname2=[noaa_obs_dir 'MSL/' stNames{id2}];
    if (exist(fname1,'dir')~=0)
       fname=[fname1 '/' stNames{id2} '.csv'];
       fileID=fopen(fname);
       C=textscan(fileID,'%s %f %d %d %d %d %d %s','Delimiter',','); 
       tt1{1,i}=C{1,1};
       %ts=find(datenum(tt1)==time(1)); te=find(datenum(tt1)==time(end));
       %eleo{1,i}=interp1(datenum(tt1(ts:te)),C{1,2}(ts:te),time);
       eleo{1,i}=C{1,2};
       fclose(fileID);
    
    elseif (exist(fname2,'dir')~=0)
       fname=[fname2 '/' stNames{id2} '.csv'];
       fileID=fopen(fname);
       C=textscan(fileID,'%s %f %d %d %d %d %d %s','Delimiter',','); 
       tt1{1,i}=C{1,1};
       %ts=find(datenum(tt1)==time(1)); te=find(datenum(tt1)==time(end));
       %eleo{1,i}=interp1(datenum(tt1(ts:te)),C{1,2}(ts:te),time);
       eleo{1,i}=C{1,2};
       fclose(fileID);
    end
    
end


for i=1:nf  
    if (isempty(eleo{1,i})==1)
       id3=find(str2double(stIds)==sa_id(i));
       fname3=[noaa_obs_dir 'PREDICTION/MSL/' stNames{id3}];
       if (exist(fname3,'dir')~=0)
           fname=[fname3 '/' stNames{id3} '.csv'];
           fileID=fopen(fname);
           C=textscan(fileID,'%s %f %d %d %d %d %d %s','Delimiter',','); 
           tt1{1,i}=C{1,1};
           eleo{1,i}=C{1,2};
           fclose(fileID);
       end
    end
end

n_nemp=nf;
eleos=eleo;
tt2=tt1;
ss_id=sa_id;ss_lat=sa_lat;ss_lon=sa_lon;

%{
n_nemp=0;%number of stations has observed values
for i=1:nf
    if (isempty(eleo{1,i})==0)      
        n_nemp=n_nemp+1;
    end
end

eleos=cell(1,n_nemp);tt2=cell(1,n_nemp);
ss_id=zeros(1,n_nemp);
ss_lat=zeros(1,n_nemp);
ss_lon=zeros(1,n_nemp);
j=0;
for i=1:nf
    if (isempty(eleo{1,i})==0)
        j=j+1;
        eleos{1,j}=eleo{1,i};
        tt2{1,j}=tt1{1,i};
        ss_id(j)=sa_id(i);
        ss_lat(j)=sa_lat(i);
        ss_lon(j)=sa_lon(i);
    end
    
end
%}

for i=1:n_nemp

    i
    
    ts=find(datenum(tt2{1,i})==time(1)); te=find(datenum(tt2{1,i})==time(end));
    if (isempty(ts)==0&&isempty(te)==0)
        tt3=tt2{1,i}(ts:te);
        elo=eleos{1,i}(ts:te);
        eleos{1,i}=interp1(datenum(tt3),elo,time);
    end
    if isempty(ts)==0&&isempty(te)==1
       id3=find(str2double(stIds)==ss_id(i));
       fname3=[noaa_obs_dir 'PREDICTION/MSL/' stNames{id3}];
       if (exist(fname3,'dir')~=0)
           fname=[fname3 '/' stNames{id3} '.csv'];
           fileID=fopen(fname);
           C=textscan(fileID,'%s %f %d %d %d %d %d %s','Delimiter',','); 
           tt1=C{1,1};
           tp1=find(datenum(tt1)==datenum(tt2{1,i}(end))); tp2=find(datenum(tt1)==time(end));
           elep=interp1([datenum(tt2{1,i}); datenum(tt1(tp1+1:tp2))],[eleos{1,i}; C{1,2}(tp1+1:tp2)],time);
           eleos{1,i}=elep;
           fclose(fileID);
       end 
    end   
    if isempty(ts)==1&&isempty(te)==0
       id3=find(str2double(stIds)==ss_id(i));
       fname3=[noaa_obs_dir 'PREDICTION/MSL/' stNames{id3}];
       if (exist(fname3,'dir')~=0)
           fname=[fname3 '/' stNames{id3} '.csv'];
           fileID=fopen(fname);
           C=textscan(fileID,'%s %f %d %d %d %d %d %s','Delimiter',','); 
           tt1=C{1,1};
           tp1=find(datenum(tt1)==time(1)); tp2=find(datenum(tt1)==datenum(tt2{1,i}(1)));
           elep=interp1([datenum(tt1(tp1:tp2-1));datenum(tt2{1,i})],[C{1,2}(tp1:tp2-1);eleos{1,i}],time);
           eleos{1,i}=elep;
           fclose(fileID);
       end 
    end   
end

%filling the missing values with predicted values
for i=1:n_nemp
    id_nan=find(isnan(eleos{1,i}));
    if (isempty(id_nan)==0)
       id3=find(str2double(stIds)==ss_id(i));
       fname3=[noaa_obs_dir 'PREDICTION/MSL/' stNames{id3}];
       if (exist(fname3,'dir')~=0)
           fname=[fname3 '/' stNames{id3} '.csv'];
           fileID=fopen(fname);
           C=textscan(fileID,'%s %f %d %d %d %d %d %s','Delimiter',','); 
           tt1=C{1,1};
           ts=find(datenum(tt1)==time(1)); te=find(datenum(tt1)==time(end));
           elep=interp1(datenum(tt1(ts:te)),C{1,2}(ts:te),time);
           eleos{1,i}(id_nan)=elep(id_nan);
           fclose(fileID);
       end
    end
end

list_con={'O1','K1','Q1','P1','K2','N2','M2','S2'}; %list of consti. to be extracted by HA
%list_con={'O1','K1','Q1','K2','N2','M2','S2'}; %list of consti. to be extracted by HA


 
for i=1:n_nemp
%Harmonic analysis
    obs=eleos{1,i};
    m2=elem(:,i)';
    [cons,xout]=t_tide(obs,'interval',1.,'start time',[start_y,start_m,start_d,0,0,0],'output',strcat(t_tide_outdir,'/obs.H.',num2str(sa_id(i))));
    [cons2,xout2]=t_tide(m2,'interval',1.,'start time',[start_y,start_m,start_d,0,0,0],'output',strcat(t_tide_outdir,'/mod.',outname,'.',num2str(sa_id(i))));
    
%low-pass filter
    obs_lo=simple_filter_data([time_ha;obs]',1./3,'low'); %3 day cut-off
    mod_lo=simple_filter_data([time_ha;m2]',1./3,'low');
end


COMP_Harmonic_plot(station_file_name,station_file_name,station_dir,prj_dir,run,run);

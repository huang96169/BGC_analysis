clear all;

%--------------inputs-----------------
%Irene
%start_datenum=datenum(2011,7,27);
%ndays=30;

%Harvey
%start_datenum=datenum(2017,8,4);
%ndays=30;

%Florence
start_datenum=datenum(2018,6,17);
ndays=30;

prj_dir='/sciclone/home20/whuang07/schism10/NWM/Case1/';  %one level up RUN*

%noaa_obs_dir='/sciclone/home20/whuang07/schism10/NWM/Case1/NOAATIDES/Harvey2017/';
%noaa_obs_dir='/sciclone/home20/whuang07/schism10/NWM/Case1/NOAATIDES/NOAA_TIDE_Irene/';
%noaa_obs_dir='/sciclone/schism10/feiye/work/NOAA_TIDE/20180823-20181001/';
noaa_obs_dir='/sciclone/home20/whuang07/schism10/NWM/Case1/NOAATIDES/Florence/';
station_dir='/sciclone/home20/whuang07/git/NWM_scripts/matlab_scripts/Elev/BPfiles/';

run='RUN1k_ZG'; 
station_file_name='coast80_6b';%'coast65_moved2';%'GOME7';%'coast80_6b';
nh=24;%hours per day from the output
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


fid=fopen([station_dir '/' station_file_name '.bp']);
[tmp]=textscan(fid,'%d',1,'headerlines',1); nf = double(tmp{1});
[tmp]=textscan(fid,'%d%f%f%f%d');
sa_lon=tmp{1,2};
sa_lat=tmp{1,3};
sa_id=tmp{1,5};

%load SCHISM output
output=load([prj_dir '/'  run '/PostP/' outname]);
dt=(output(2,1)-output(1,1));%/86400;
time=start_datenum+double(output(:,1));%/86400;
time=time(nh*1+1:nh*(ndays+1));
%time=time(1:2:end);
elem=output(nh*1+1:nh*(ndays+1),2:end);
%elem=elem(1:2:end,:);

start_y=str2num(datestr(time(1),'yyyy'));
start_m=str2num(datestr(time(1),'mm'));
start_d=str2num(datestr(time(1),'dd'));


eleo=cell(1,nf);
tt1=cell(1,nf);

for i=1:nf
    i
    id2=find(str2double(stIds)==sa_id(i));
    fname1=[noaa_obs_dir 'NAVD/' stNames{id2}];
    fname2=[noaa_obs_dir 'MSL/' stNames{id2}];
    if (exist(fname1,'dir')~=0)
       fname=[fname1 '/' stNames{id2} '.csv'];
       fileID=fopen(fname);
       C=textscan(fileID,'%s %f %d %d %d %d %d %s','Delimiter',','); 
       tmp=[cell2mat(C{1,1}) repmat(':00',size(cell2mat(C{1,1}),1),1)];
       tmp2=cellstr(string(tmp));
       tmp3=DateStr2Num(cellstr(tmp2),31);
       tin=find(tmp3>=time(1) & tmp3<=time(end));
       if(isempty(tin)==0)
         eleo{1,i}=interp1(tmp3(tin),C{1,2}(tin),time);
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
       tin=find(tmp3>=time(1) & tmp3<=time(end));
       if(isempty(tin)==0&&isempty(find(isnan(C{1,2}(tin))))==1)
         eleo{1,i}=interp1(tmp3(tin),C{1,2}(tin),time);
         tt1{1,i}=C{1,1}(tin);
       else
         eleo{1,i}=[];
         tt1{1,i}=[];
       end
       fclose(fileID);
    end
    
end



n_nemp=nf;
eleos=eleo;
tt2=tt1;
ss_id=sa_id;ss_lat=sa_lat;ss_lon=sa_lon;



list_con={'O1','K1','Q1','P1','K2','N2','M2','S2'}; %list of consti. to be extracted by HA
%list_con={'O1','K1','Q1','K2','N2','M2','S2'}; %list of consti. to be extracted by HA
%list_con{'O1','K1','Q1','N2','M2','S2'};

%system(['rm ' t_tide_outdir '/*']); 
for i=1:n_nemp
%Harmonic analysis
  
  obs=eleos{1,i};
  if (~isempty(obs))
    m2=elem(:,i)';
    [cons,xout]=t_tide(obs,'interval',1.0,'start time',[start_y,start_m,start_d,0,0,0],'output',strcat(t_tide_outdir,'/obs.H.',num2str(sa_id(i))));
    [cons2,xout2]=t_tide(m2,'interval',1.0,'start time',[start_y,start_m,start_d,0,0,0],'output',strcat(t_tide_outdir,'/mod.',outname,'.',num2str(sa_id(i))));
  end
%low-pass filter
    %obs_lo=simple_filter_data([time_ha;obs]',1./3,'low'); %3 day cut-off
    %mod_lo=simple_filter_data([time_ha;m2]',1./3,'low');
end
COMP_Harmonic_plot_1(station_file_name,station_file_name,station_dir,prj_dir,run,run);

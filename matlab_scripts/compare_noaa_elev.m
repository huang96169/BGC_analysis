%comparison of model output and obs for the 63 coastal stations
clear all
%-----------------------inputs---------------------------
%start_datenum=datenum(2017,8,4);%Harvey
start_datenum=datenum(2018,6,17);%Irene
ndays=50;

runo='CB01a'; 
station_file_nameo='CBDB';%'upbay_sta53.moved2';%'upbay_sta53';%coast63.moved6 

run='CB01a'; 
station_file_name='CBDB';%'upbay_sta53.moved2';%'upbay_sta53';%coast63.moved5 

%dirs
prj_dir='/sciclone/home20/whuang07/schism10/NWM/Case1/';  %one level up RUN*

%noaa_obs_dir='/sciclone/home20/whuang07/schism10/NWM/Case1/NOAATIDES/Harvey2017/';
%noaa_obs_dir='~feiye/schism10/work/NOAA_TIDE/';
noaa_obs_dir='/sciclone/home20/whuang07/schism10/NWM/Case1/NOAATIDES/Florence/';
%noaa_obs_dir='/sciclone/home20/whuang07/schism10/NWM/Case1/NOAATIDES/NOAA_TIDE_Irene/';
station_dir='/sciclone/home20/whuang07/git/NWM_scripts/matlab_scripts/Elev/BPfiles/';
saved_mat_dir='/sciclone/schism10/whuang07/NWM/Case1/SAVED_MAT/'; %


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
dt = round((output(2,1)-output(1,1))*86400/3600)/24; %round to integer hours, then /24
time_model=start_datenum:dt:start_datenum+dt*(length(output(1:ndays*24,1))-1); 

elem2=output(1:ndays*24,2:end);
tt=output(1:ndays*24,1);

%load SCHISM output to be compared
output=load([prj_dir '/' runo '/PostP/' outnameo]);

elem1=output(1:length(time_model),2:end);

eleo=cell(1,nf);
tt1=cell(1,nf);

ff=fopen([prj_dir '/NOAAtidecompare.' outname '.' outnameo '.log'],'w');
fprintf(ff,'%s\n',[outname outnameo]);

for i=1:nf
    
    id2=find(str2double(stIds)==sa_id(i));
    fname1=[noaa_obs_dir '/NAVD/' stNames{id2}];
    fname2=[noaa_obs_dir '/MSL/' stNames{id2}];
    if (exist(fname1,'dir')~=0)
       fprintf(ff,'%d %d %s\n',[i sa_id(i) 'datum: NAVD']);
       fname=[fname1 '/' stNames{id2} '.csv'];
       fileID=fopen(fname);
       C=textscan(fileID,'%s %f %d %d %d %d %d %s','Delimiter',','); 
       tt1{1,i}=C{1,1};
       %ts=find(datenum(tt1)==time_model(1)); te=find(datenum(tt1)==time_model(end));
       %eleo{1,i}=interp1(datenum(tt1(ts:te)),C{1,2}(ts:te),time_model);
       eleo{1,i}=C{1,2};
       fclose(fileID);
    %end
    elseif (exist(fname2,'dir')~=0)
       fprintf(ff,'%d %d %s\n',[i sa_id(i) 'datum: MSL']);
       fname=[fname2 '/' stNames{id2} '.csv'];
       fileID=fopen(fname);
       C=textscan(fileID,'%s %f %d %d %d %d %d %s','Delimiter',','); 
       tt1{1,i}=C{1,1};
       %ts=find(datenum(tt1)==time_model(1)); te=find(datenum(tt1)==time_model(end));
       %eleo{1,i}=interp1(datenum(tt1(ts:te)),C{1,2}(ts:te),time_model);
       eleo{1,i}=C{1,2};
       fclose(fileID);
    end

end

%
% fill the missing eleo{1,i} with the predicted values
for i=1:nf  
    if (isempty(eleo{1,i})==1)
       id3=find(str2double(stIds)==sa_id(i));
       fname3=[noaa_obs_dir '/PREDICTION/MSL/' stNames{id3}];
       if (exist(fname3,'dir')~=0)
           fprintf(ff,'%d %d %s\n',[i sa_id(i) 'Filled with PREDICTION/MSL']);
           fname=[fname3 '/' stNames{id3} '.csv'];
           fileID=fopen(fname);
           C=textscan(fileID,'%s %f %d %d %d %d %d %s','Delimiter',','); 
           tt1{1,i}=C{1,1};
           eleo{1,i}=C{1,2};
           fclose(fileID);
       end
    end
end
%}




%if eleo{1,i} is empty, means sa_id has no observation only has prediction.
%folloing is to remove the sa_id(i) that has no observation

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
n_nemp=nf;
eleos=eleo; 
tt2=tt1;
ss_id=sa_id;ss_lat=sa_lat;ss_lon=sa_lon;

tmp=intersect(station_file_name,station_file_nameo,'stable');
fName_savedObs=[saved_mat_dir '/' tmp '_' datestr(time_model(1),'yyyymmdd') '-' datestr(time_model(end),'yyyymmdd') '.mat']; %###not a good name, use bpfilename
if ~exist (fName_savedObs,'file')
    for i=1:n_nemp
        i
        ts=find(datenum(tt2{1,i})==time_model(1)); te=find(datenum(tt2{1,i})==time_model(end));
        if (isempty(ts)==0&&isempty(te)==0)
            tt3=tt2{1,i}(ts:te);
            elo=eleos{1,i}(ts:te);
            eleos{1,i}=interp1(datenum(tt3),elo,time_model);
        end
        %###print warning msg
        if isempty(ts)==0&&isempty(te)==1
           id3=find(str2double(stIds)==ss_id(i));
           fname3=[noaa_obs_dir '/PREDICTION/MSL/' stNames{id3}];
           if (exist(fname3,'dir')~=0)
               fname=[fname3 '/' stNames{id3} '.csv'];
               fileID=fopen(fname);
               C=textscan(fileID,'%s %f %d %d %d %d %d %s','Delimiter',','); 
               tt1=C{1,1};
               tp1=find(datenum(tt1)==datenum(tt2{1,i}(end))); tp2=find(datenum(tt1)==time_model(end));
               elep=interp1([datenum(tt2{1,i}); datenum(tt1(tp1+1:tp2))],[eleos{1,i}; C{1,2}(tp1+1:tp2)],time_model);
               eleos{1,i}=elep;
               fclose(fileID);
           end 
        end   
        if isempty(ts)==1&&isempty(te)==0
           id3=find(str2double(stIds)==ss_id(i));
           fname3=[noaa_obs_dir '/PREDICTION/MSL/' stNames{id3}];
           if (exist(fname3,'dir')~=0)
               fname=[fname3 '/' stNames{id3} '.csv'];
               fileID=fopen(fname);
               C=textscan(fileID,'%s %f %d %d %d %d %d %s','Delimiter',','); 
               tt1=C{1,1};
               tp1=find(datenum(tt1)==time_model(1)); tp2=find(datenum(tt1)==datenum(tt2{1,i}(1)));
               elep=interp1([datenum(tt1(tp1:tp2-1));datenum(tt2{1,i})],[C{1,2}(tp1:tp2-1);eleos{1,i}],time_model);
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
           fname3=[noaa_obs_dir '/PREDICTION/MSL/' stNames{id3}];
           if (exist(fname3,'dir')~=0)
               fname=[fname3 '/' stNames{id3} '.csv'];
               fileID=fopen(fname);
               C=textscan(fileID,'%s %f %d %d %d %d %d %s','Delimiter',','); 
               tt1=C{1,1};
               ts=find(datenum(tt1)==time_model(1)); te=find(datenum(tt1)==time_model(end));
               elep=interp1(datenum(tt1(ts:te)),C{1,2}(ts:te),time_model);
               eleos{1,i}(id_nan)=elep(id_nan);
               fclose(fileID);
           end
        end
    end
    save(fName_savedObs,'eleos');
else
    load(fName_savedObs);
end
 



%
figure  
dn=1;
cc2=zeros(1,n_nemp);
rmse2=zeros(1,n_nemp);
mae2=zeros(1,n_nemp);
SS2=zeros(1,n_nemp);
cc1=zeros(1,n_nemp);
rmse1=zeros(1,n_nemp);
mae1=zeros(1,n_nemp);
SS1=zeros(1,n_nemp);
ncol=ceil(nf/13);

n=0;
for i=1:n_nemp
    %h(i)=subplot(13,ncol,i);
    %plot(tt,detrend(eleos{1,i}),'r');hold on;plot(tt,detrend(elem1(:,i)),'g');hold on;
    %plot(tt,detrend(elem2(:,i)),'b');

    %if (i==1) 
    %  legend({'obs',outnameo,outname});
    %end

    mean_erro(i)=mean(elem1(dn:end,i)'-eleos{1,i}(dn:end) ); %not detrended, used for adjusting elev ocean boundary
    mean_err(i)=mean(elem2(dn:end,i)'-eleos{1,i}(dn:end) ); %not detrended, used for adjusting elev ocean boundary
    aa1=corrcoef(detrend(elem2(dn:end,i)),detrend(eleos{1,i}(dn:end)));
    cc2(i)=round(aa1(1,2),2);
    rmse2(i)=sqrt(sum((detrend(eleos{1,i}(dn:end))-detrend(elem2(dn:end,i)')).^2)/length(time_model(dn:end)));
    mae=sum(abs(detrend(eleos{1,i}(dn:end))-detrend(elem2(dn:end,i)')))/length(time_model(dn:end));
    mae2(i)=round(mae,2);
    SS2(i)=1-sum((detrend(eleos{1,i}(dn:end))-detrend(elem2(dn:end,i)')).^2)/sum((detrend(eleos{1,i}(dn:end))-mean(eleos{1,i}(dn:end))).^2);
    aa2=corrcoef(detrend(elem1(dn:end,i)),detrend(eleos{1,i}(dn:end)));
    cc1(i)=round(aa2(1,2),2);
    rmse1(i)=sqrt(sum((detrend(eleos{1,i}(dn:end))-detrend(elem1(dn:end,i)')).^2)/length(time_model(dn:end)));
    mae=sum(abs(detrend(eleos{1,i}(dn:end))-detrend(elem1(dn:end,i)')))/length(time_model(dn:end));
    mae1(i)=round(mae,2);
    SS1(i)=1-sum((detrend(eleos{1,i}(dn:end))-detrend(elem1(dn:end,i)')).^2)/sum((detrend(eleos{1,i}(dn:end))-mean(eleos{1,i}(dn:end))).^2);
    if (cc2(i)~=0.0)
       n=n+1;
       h(n)=subplot(13,ncol,n);    
       %plot(tt,(detrend(eleos{1,i})),'r');hold on;
       %plot(tt,(detrend(elem1(:,i))),'g'); 
       %plot(tt,(detrend(elem2(:,i))),'b'); 
       plot(tt,((eleos{1,i})),'r');hold on;
       plot(tt,((elem1(:,i))),'g');
       plot(tt,((elem2(:,i))),'b');
       if (n==1)
          legend({'obs',outnameo,outname});
       end
       title(strcat('CC=',num2str(cc2(i)),'MAE=',num2str(mae2(i)),' ID=',num2str(sa_id(i))))
       xlim([10 ndays-1]);
    end    
    %format short
    %title(strcat('CC=',num2str(cc2(i)),'MAE=',num2str(mae2(i)),' ID=',num2str(sa_id(i))))
    %title(strcat('CC=',num2str(cc),' ID=',num2str(i)))
    %xlim([10 ndays-1]);
        
end

fprintf(ff,'%s\n','station# station_id MAE_oldrun MAE_newrun');
for i=1:n_nemp
    
    fprintf(ff,'%d %d %f %f\n',i, sa_id(i), mae1(i), mae2(i));
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

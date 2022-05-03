%comparison of model output and obs for the 63 coastal stations
clear all
%-----------------------inputs---------------------------
%start_datenum=datenum(2017,8,4);%Harvey
start_datenum=datenum(2021,07,3);%Irene
ndays=8;
%dirs
prj_dir='/sciclone/home20/whuang07/schism10/NWM/Case1/NOAATIDES/Forecast/results3D/';  %one level up RUN*

%noaa_obs_dir='/sciclone/home20/whuang07/schism10/NWM/Case1/NOAATIDES/Harvey2017/';
noaa_obs_dir='/sciclone/home20/whuang07/schism10/NWM/Case1/NOAATIDES/Forecast/';
station_dir='/sciclone/home20/whuang07/git/NWM_scripts/matlab_scripts/Elev/BPfiles/';
saved_mat_dir='/sciclone/home20/whuang07/schism10/NWM/Case1/NOAATIDES/Forecast/results3D/'; %

%outname=['elev.' station_file_name '.' run];
%outnameo=['elev.' station_file_nameo '.' runo];
time_model=start_datenum+1/48:1/48:start_datenum+ndays;
tt=1/48:1/48:ndays;
for i=1:ndays
   outname=[datestr(start_datenum+i-1,'yyyymmdd')]
   if (exist([prj_dir '/fcst/fcst/'  outname],'dir')~=0)
     output=load([prj_dir '/fcst/fcst/'  outname '/staout_1']);
     dt = round((output(2,1)-output(1,1))/3600)/24; %round to integer hours, then /24
     %time_model=start_datenum:dt:start_datenum+dt*(length(output(:,1))-1);
     elem(1+(i-1)*48:i*48,:)=output(49:96,2:end); 
   else
     elem(1+(i-1)*48:i*48,:)=NaN;
   end
end
%runo='RUN20201128'; 
%station_file_nameo='forecast143';%'upbay_sta53.moved2';%'upbay_sta53';%coast63.moved6 
%
%run='RUN20201128'; 
station_file_name='Coast_6b';%'upbay_sta53.moved2';%'upbay_sta53';%coast63.moved5 
%station_s=[8768094 8766072 8764314 8762075 8761724 8761927 8761305 8747437 8735180 8729210]; %for TS Claudette
station_s=[8724580 8723970 8725110 8725520 8726384 8726520 8726724 8727520 8728690];%for TS Elsa

station_s=station_s';
%---------------------end inputs-------------------------



% station id and name  
f1=fopen([station_dir '/stations.txt']);
[tmp]=textscan(f1,'%s%s','delimiter',',');
stIds=tmp{1,1};
stNames=tmp{1,2};
fclose(f1);

% read lat/lon for stations 
fid=fopen([station_dir '/' station_file_name '.bp']);
[tmp]=textscan(fid,'%d',1,'headerlines',1); 
[tmp]=textscan(fid,'%d%f%f%f%d');
sa_lon=tmp{1,2};
sa_lat=tmp{1,3};
sa_id=tmp{1,5};
%only keep selected stations:

[dummy,ia,ib]=intersect(sa_id,station_s);
sa_lon=sa_lon(ia);
sa_lat=sa_lat(ia);
sa_id=sa_id(ia);
nf=length(sa_id);
elem=elem(:,ia);
%-------only keep data at selected stations-------------

%load SCHISM output
%output=load([prj_dir '/'  outname]);
%dt = round((output(2,1)-output(1,1))/3600)/24; %round to integer hours, then /24
%time_model=start_datenum:dt:start_datenum+dt*(length(output(:,1))-1); 

%elem2=output(:,2:end);
%tt=output(:,1);

%load SCHISM output to be compared
%output=load([prj_dir '/' outnameo]);

%elem1=output(1:length(time_model),2:end);

eleo=cell(1,nf);
tt1=cell(1,nf);

ff=fopen([prj_dir '/NOAAtidecompare.' outname '.log'],'w');
fprintf(ff,'%s\n',[outname]);
VD=cell(1,nf);
staname=cell(1,nf);
for i=1:nf

    id2=find(str2double(stIds)==sa_id(i));
    fname1=[noaa_obs_dir '/NAVD/' stNames{id2}];
    fname2=[noaa_obs_dir '/MSL/' stNames{id2}];
    if (exist(fname1,'dir')~=0)
       fprintf(ff,'%d %d %s\n',[i sa_id(i) 'datum: NAVD']);
       VD{1,i}=cellstr('NAVD');
       staname{1,i}=stNames{id2};
       fname=[fname1 '/' stNames{id2} '.csv'];
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
       staname{1,i}=stNames{id2};
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
         idn1=find(isnan(eleo{1,i}));
         if(isempty(idn1)==0)
           idn2=find(~isnan(eleo{1,i}));
           %tt1{1,i}=tt1{1,i}(idn2);
           tmp3=time_model(idn2);
           eleo{1,i}=eleo{1,i}(idn2);
           eleo{1,i}=interp1(tmp3,eleo{1,i},time_model);
         end
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

%tmp=intersect(station_file_name,station_file_nameo,'stable');
ncol=ceil(nf/10); 
color=['b'; 'g';'y';'m';'w';'k' ];

for j=1:1
figure
  for i=1:ndays
    outname=[datestr(start_datenum+i-1,'yyyymmdd')]
    if (exist([prj_dir '/fcst/fcst/'  outname],'dir')~=0)
      output=load([prj_dir '/fcst/fcst/'  outname '/staout_1']);
      dt = ((output(2,1)-output(1,1))/3600)/24; %round to integer hours, then /24
    %time_model=start_datenum:dt:start_datenum+dt*(length(output(:,1))-1);
      time3d=(i-2)+dt:dt:i+1;
      elem3d=output(:,1+ia);
    else
      time3d=(i-2)+1/48:1/48:i+1;
      elem3d(1:48*3)=NaN;
    end
    for nn=(j-1)*16+1:j*16
      if(nn<=nf)
      %if(nn<=nf&&isempty(find(id_s==nn))==0)
        %subplot(5,6,nn-(j-1)*30);
        subaxis(8,2,nn-(j-1)*16,'Spacing', 0.05,'MarginTop',0.05,'MarginBottom',0.01,'MarginLeft',0.03,'MarginRight',0.01);
        tmpo=eleo{1,nn};
        if(isempty(eleo{1,nn})==0&&isempty(tt1{1,nn})==0)
          if(i==1)
            %plot(start_datenum+tt,detrend(tmpo),'r');hold on;
            plot(start_datenum+tt,(tmpo),'r');hold on;
          end
          %plot(start_datenum+time3d,detrend(elem3d(:,nn)),'b');hold on
          plot(start_datenum+time3d,(elem3d(:,nn)),'b');hold on
          %if (nn==1&i==1)
          %  legend({'obs','model'});
          %end
        else
          plot(start_datenum+time3d,((elem3d(:,nn))),'b');
        end
        datetick('x')
        %xticks([0:0.5:ndays+1])
        %xticklabels({num2str(0:0.5:ndays+1)})
        a = get(gca,'XTickLabel');
        set(gca,'XTickLabel',a,'fontsize',5)
        title({[num2str(nn),': ',num2str(sa_id(nn)),' VD ',cell2mat(VD{1,nn})]
              [(staname{1,nn})]},'Interpreter','none','FontSize',6)
        hold on;
        %subaxis(5,6,nn-(j-1)*30,'Margin', 0);
      end
    end%nn 
  end%i
  name_fig=strcat(prj_dir,'/Figures/','selected_elev.3days_',datestr(start_datenum,'yyyymmdd'),'_',datestr(time_model(end),'yyyymmdd'),num2str(j));
  saveas(gcf,[name_fig '.png'])
  savefig([name_fig '.fig'])
end%j
%
dn=1;
cc2=zeros(1,n_nemp);
rmse2=zeros(1,n_nemp);
mae2=zeros(1,n_nemp);
SS2=zeros(1,n_nemp);
cc1=zeros(1,n_nemp);
rmse1=zeros(1,n_nemp);
mae1=zeros(1,n_nemp);
SS1=zeros(1,n_nemp);
ncol=ceil(nf/10);

n=0;
for j=1:1
  figure
  for i=(j-1)*16+1:j*16
    i
    if(i<=nf)
    %if(i<=nf&&isempty(find(id_s==i))==0)
      %subplot(5,6,i-(j-1)*30);
      subaxis(8,2,i-(j-1)*16,'Spacing', 0.05,'MarginTop',0.05,'MarginBottom',0.01,'MarginLeft',0.03,'MarginRight',0.01);
      tmpo=eleo{1,i};
      if(isempty(eleo{1,i})==0&&isempty(tt1{1,i})==0)
        id=find(isnan(tmpo)==0&isnan(elem(:,i)')==0);
        aa=corrcoef(tmpo(id),elem(id,i));
        cc=round(aa(1,2),2);
        mae=sum(abs(detrend(tmpo(id))-detrend(elem(id,i)')))/length(time_model(id));
        mae2(i)=round(mae,2);
        %plot(time_model,detrend(tmpo),'r');hold on;
        %plot(time_model,detrend((elem(:,i))),'b');
        plot(time_model,(tmpo),'r');hold on;
        plot(time_model,((elem(:,i))),'b');
        n=n+1;
        %if (n==1)
        %   legend({'obs','model'});
        %end
 
      else
        mae2(i)=NaN;
        cc=NaN;
        plot(time_model,((elem(:,i))),'b');
      end
        %xticks([0:0.5:ndays+1])
        %xticklabels({num2str(0:0.5:ndays+1)})
        datetick('x')
        a = get(gca,'XTickLabel');
        set(gca,'XTickLabel',a,'fontsize',5)
        title({['CC:',num2str(cc),'MAE:',num2str(mae2(i)),' VD:',cell2mat(VD{1,i})] 
              [num2str(i),': ',num2str(sa_id(i)),staname{1,i}]},'Interpreter','none','FontSize',6)
        xlim([time_model(1) time_model(end)]);
    end
  end
  name_fig=strcat(prj_dir,'/Figures','/selected_elev_',datestr(time_model(1),'yyyymmdd'),'_',datestr(time_model(end),'yyyymmdd'),'_',num2str(j));
  saveas(gcf,[name_fig '.png'])
  savefig([name_fig '.fig'])
end
fprintf(ff,'%s\n','station# station_id MAE_oldrun MAE_newrun');
for i=1:n_nemp
    
    fprintf(ff,'%d %d %f\n',i, sa_id(i), mae2(i));
end
fclose(ff);

%low-pass filter: simple_filter_data([time_ha;obs]',1./3,'low');

fN=1/(2*dt);
[bf,af]=butter(6,0.6/fN);
n=0;
for j=1:1
  figure
  for i=(j-1)*16+1:j*16
    i
    if(i<=nf)
    %if(i<=nf&&isempty(find(id_s==i))==0)
      %subplot(5,6,i-(j-1)*30);
      subaxis(8,2,i-(j-1)*16,'Spacing', 0.05,'MarginTop',0.05,'MarginBottom',0.01,'MarginLeft',0.03,'MarginRight',0.01);
      tmpo=eleo{1,i};
      if(isempty(eleo{1,i})==0&&isempty(tt1{1,i})==0)
        id=find(isnan(tmpo)==0&isnan(elem(:,i)')==0);
        aa=corrcoef(tmpo(id),elem(id,i));
        cc=round(aa(1,2),2);
        mae=sum(abs(detrend(tmpo(id))-detrend(elem(id,i)')))/length(time_model(id));
        mae2(i)=round(mae,2);
        %plot(time_model,detrend(filtfilt(bf,af,double(tmpo))),'r');hold on;
        %plot(time_model,detrend(filtfilt(bf,af,double(elem(:,i)))),'b');
        plot(time_model,(filtfilt(bf,af,double(tmpo))),'r');hold on;
        plot(time_model,(filtfilt(bf,af,double(elem(:,i)))),'b');
        n=n+1;
        %if (n==1)
        %   legend({'obs','model'});
        %end

      else
        mae2(i)=NaN;
        cc=NaN;
        plot(time_model,filtfilt(bf,af,double(elem(:,i))),'b');
      end
        datetick('x')
        a = get(gca,'XTickLabel');
        set(gca,'XTickLabel',a,'fontsize',5)
        title({['CC:',num2str(cc),'MAE:',num2str(mae2(i)),' VD:',cell2mat(VD{1,i})]
              [num2str(i),': ',num2str(sa_id(i)),staname{1,i}]},'Interpreter','none','FontSize',6)
        xlim([time_model(1) time_model(end)]);
    end
  end
  name_fig=strcat(prj_dir,'/Figures','/selected_lpelev_',datestr(time_model(1),'yyyymmdd'),'_',datestr(time_model(end),'yyyymmdd'),'_',num2str(j));
  saveas(gcf,[name_fig '.png'])
  savefig([name_fig '.fig'])
end




%set(h(1),'position',[0.01,0.9,0.2,0.1])

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
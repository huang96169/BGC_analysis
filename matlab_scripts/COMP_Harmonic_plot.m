function COMP_Harmo(varargin)
clc;
%clear all
%-----------------inputs------------------------------
if nargin==0  %change inputs here if running this script independently
  station_file_name='coast80'; 
  station_file_nameo='coast80'; 
  station_dir='/sciclone/home20/whuang07/git/NWM_scripts/matlab_scripts/Elev/BPfiles/';
  prj_dir='/sciclone/home20/whuang07/schism10/NWM/Case1/';  %one level up RUN*
  run='RUN1_ZG'; 
  runo='RUN1_ZG'; 
else %script invoked somewhere else
  iArg=1;
  station_file_name=varargin{1,iArg}; iArg=iArg+1;
  station_file_nameo=varargin{1,iArg}; iArg=iArg+1;
  station_dir=varargin{1,iArg}; iArg=iArg+1;
  prj_dir=varargin{1,iArg}; iArg=iArg+1;
  run=varargin{1,iArg}; iArg=iArg+1;
  runo=varargin{1,iArg}; iArg=iArg+1;
end
%-----------------end inputs-------------------



t_tide_outdir=[prj_dir '/T_Tide_out/'];
t_tide_outdiro=[prj_dir '/T_Tide_out/'];

outname=['elev.' station_file_name '.' run];
outnameo=['elev.' station_file_nameo '.' runo];

%plot phase and amplitude after T-Tide analysis
% read lat/lon for stations 
fid=fopen([station_dir '/' station_file_name '.bp']);
[tmp]=textscan(fid,'%d',1,'headerlines',1); nf = double(tmp{1});
[tmp]=textscan(fid,'%d%f%f%f%d');
sa_lon=tmp{1,2};
sa_lat=tmp{1,3};
sa_id=tmp{1,5};

nst=length(sa_id);
amp_m1=cell(1,nst);pha_m1=cell(1,nst);
amp_m2=cell(1,nst);pha_m2=cell(1,nst);
amp_o=cell(1,nst);pha_o=cell(1,nst);
figure

for i=1:nst 
   
  fname=[t_tide_outdiro '/mod.' outnameo '.' num2str(sa_id(i))];
  if  (exist(fname)~=0)
    fileID=fopen(fname);
    mod=textscan(fileID,'%s %f %f %f %f %f %f','headerlines',16); 
    amp_m1{1,i}=mod{1,3};
    pha_m1{1,i}=mod{1,5};   
    fclose(fileID);
  else
    amp_m1{1,i}=NaN;
    pha_m1{1,i}=NaN;
  end

  fname=[t_tide_outdir '/mod.' outname '.' num2str(sa_id(i))];
  if (exist(fname)~=0)
    fileID=fopen(fname);
    mod=textscan(fileID,'%s %f %f %f %f %f %f','headerlines',16); 
    amp_m2{1,i}=mod{1,3};
    pha_m2{1,i}=mod{1,5};   
    fclose(fileID);
  else
    amp_m2{1,i}=NaN;
    pha_m2{1,i}=NaN;
  end  
  fname=[t_tide_outdir '/obs.H.' num2str(sa_id(i))];
  if (exist(fname)~=0)
    i
    fileID=fopen(fname);
    obs=textscan(fileID,'%s %f %f %f %f %f %f','headerlines',16); 
    %com=obs{1,1};
    amp_o{1,i}=obs{1,3};
    pha_o{1,i}=obs{1,5};   
    fclose(fileID);
  else
    amp_o{1,i}=NaN;
    pha_o{1,i}=NaN;
  end
    %M2 15; K1 8; O1 6; S2 17
    am=17;
    cm='S2';
    plot(i,amp_o{1,i}(am),'ro')
    hold on
    plot(i,amp_m1{1,i}(am),'go')
    hold on
    plot(i,amp_m2{1,i}(am),'bo')
    hold on
    mae(i)=abs(amp_o{1,i}(am)-amp_m1{1,i}(am));
    factor1=(amp_o{1,i}(am)*cosd(pha_o{1,i}(am))-amp_m1{1,i}(am)*cosd(pha_m1{1,i}(am)))^2;
    factor2=(amp_o{1,i}(am)*sind(pha_o{1,i}(am))-amp_m1{1,i}(am)*sind(pha_m1{1,i}(am)))^2;
    complex_e(i)=sqrt(factor1+factor2);
  
end
size(complex_e)
fprintf('%s',['complex error for ',cm,' tide is'])
mean(complex_e)
xticks(1:1:nst)
grid on
%xticklabels(num2str(st58_id))
title([cm 'amplitude'])
%set(gca, 'YScale', 'log')

savefig(strcat(prj_dir,'/Figures/',cm,'tide',num2str(nst),outname,outnameo,'.fig'))


figure;
mae=zeros(1,nst);
ncomp=am;%O1,6; K1,8; M2,15; S2,17;
for i=1:nst
    
    mae(i)=abs(pha_o{1,i}(ncomp)-pha_m1{1,i}(ncomp));
    if mae(i)>180
        i
        pha_m1{1,i}(ncomp)=360-pha_m1{1,i}(ncomp);
    end
    scatter(i,pha_o{1,i}(ncomp),'or'); hold on;
    scatter(i,pha_m1{1,i}(ncomp),'+b'); hold on;
end

for i=1:nst
    plot([i i],[pha_o{1,i}(ncomp) pha_m1{1,i}(ncomp)],'k'); hold on;
end
title([cm 'tidal phase']);
grid on; box on;
xlabel('Station #');
ylabel('Phase');
xlim([0 nst]);
xticks(1:1:nst);
legend('Obs','Model');

savefig(strcat(prj_dir,'/Figures/',cm,'pha',num2str(nst),outname,outnameo,'.fig'))

%figure
%for i=1:nst
%    geoplot(sa_lat(i),sa_lon(i),'ro')
    
%    text(sa_lat(i),sa_lon(i),num2str(i),'FontSize',8)
%    hold on
%end
%}

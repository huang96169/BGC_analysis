function COMP_Harmo(varargin)
clc;
%clear all
%-----------------inputs------------------------------
if nargin==0  %change inputs here if running this script independently
  station_file_name='coast80_6b'; 
  station_file_nameo='coast80_6b'; 
  station_dir='/sciclone/home20/whuang07/git/NWM_scripts/matlab_scripts/Elev/BPfiles/';
  prj_dir='/sciclone/home20/whuang07/schism10/NWM/Case1/';  %one level up RUN*
  run='RUN1k1'; 
  runo='RUN1k1'; 
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
nn=0;
for i=1:nst 
   
  fname=[t_tide_outdiro '/mod.' outnameo '.' num2str(sa_id(i))];
  if  (exist(fname)~=0)
    fileID=fopen(fname);
    mod=textscan(fileID,'%s %f %f %f %f %f %f','headerlines',16); 
    amp_m1{1,i}=mod{1,3};
    pha_m1{1,i}=mod{1,5};   
    fclose(fileID);
  else
    amp_m1{1,i}=[];
    pha_m1{1,i}=[];
  end

  fname=[t_tide_outdir '/mod.' outname '.' num2str(sa_id(i))]
  if (exist(fname)~=0)
    fileID=fopen(fname);
    %fname
    mod=textscan(fileID,'%s %f %f %f %f %f %f','headerlines',16); 
    amp_m2{1,i}=mod{1,3};
    pha_m2{1,i}=mod{1,5};   
    fclose(fileID);
  else
    amp_m2{1,i}=[];
    pha_m2{1,i}=[];
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
    amp_o{1,i}=[];
    pha_o{1,i}=[];
  end
    %M2 15; K1 8; O1 6; S2 17
    %am=15;
    cm='M2';
   
  %for nc=1:29
  %  mod{1,1}(nc)
  %  if(contains(mod{1,1}(nc),cm))
  %    am=nc
  %  end
  %end
  if (isempty(amp_o{1,i})==0&&isempty(amp_m1{1,i})==0)
    nn=nn+1;
    for nc=1:29
       mod{1,1}(nc)
       if(contains(mod{1,1}(nc),cm))
           am=nc
       end
    end
    amp_oo(nn)=amp_o{1,i}(am);
    amp_mm1(nn)=amp_m1{1,i}(am);
    amp_mm2(nn)=amp_m2{1,i}(am);
    plot(i,amp_o{1,i}(am),'ro','MarkerSize',5,'MarkerEdgeColor','r','MarkerFaceColor','r')
    hold on
    plot(i,amp_m1{1,i}(am),'g*','MarkerSize',10)
    hold on
    plot(i,amp_m2{1,i}(am),'b+','MarkerSize',10)
    hold on
    mae(nn)=abs(amp_o{1,i}(am)-amp_m1{1,i}(am));
    mer(nn)=amp_o{1,i}(am)-amp_m1{1,i}(am);
    %factor1=(amp_o{1,i}(am)*cosd(pha_o{1,i}(am))-amp_m1{1,i}(am)*cosd(pha_m1{1,i}(am)))^2;
    %factor2=(amp_o{1,i}(am)*sind(pha_o{1,i}(am))-amp_m1{1,i}(am)*sind(pha_m1{1,i}(am)))^2;
    %complex_e(i)=sqrt(factor1+factor2);
  end
end
fprintf('%s',['MAE for ',cm,' tide is'])
mean(mae)
fprintf('%s',['mean error for ',cm,' tide is'])
mean(mer)
%size(complex_e)
%fprintf('%s',['complex error for ',cm,' tide is'])
%mean(complex_e)
xticks(1:1:nst)
grid on
%xticklabels(num2str(st58_id))
title([cm 'amplitude'])
legend('Obs',runo,run);
%set(gca, 'YScale', 'log')

savefig(strcat(prj_dir,'/Figures/',cm,'tide',num2str(nst),outname,outnameo,'.fig'))
saveas(gcf,strcat(prj_dir,'/Figures/',cm,'tide',num2str(nst),outname,outnameo,'.png'))
%figure;
%for i=1:nst
%   plot(amp_o{1,i}(am),amp_m2{1,i}(am),'bo');hold on;
%end
%savefig(strcat(prj_dir,'/Figures/',cm,'scatter',num2str(nst),outname,'.fig'))

figure;
mae1=zeros(1,nst);
mae2=zeros(1,nst);
ncomp=am;%O1,6; K1,8; M2,15; S2,17;
mm=0;
for i=1:nst

  if (isempty(pha_o{1,i})==0&&isempty(pha_m1{1,i})==0)  
    mae1(i)=abs(pha_o{1,i}(ncomp)-pha_m1{1,i}(ncomp));
    if mae1(i)>180
        pha_m1{1,i}(ncomp)=360-pha_m1{1,i}(ncomp);        
    end
    mae2(i)=abs(pha_o{1,i}(ncomp)-pha_m2{1,i}(ncomp));
    if mae2(i)>180        
        pha_m2{1,i}(ncomp)=360-pha_m2{1,i}(ncomp);        
    end
    mm=mm+1;
    pha_oo(mm)=pha_o{1,i}(am);
    pha_mm1(mm)=pha_m1{1,i}(am);
    pha_mm2(mm)=pha_m2{1,i}(am);
    scatter(i,pha_o{1,i}(ncomp),30,'or','MarkerEdgeColor','r','MarkerFaceColor','r'); hold on;
    scatter(i,pha_m1{1,i}(ncomp),40,'g*'); hold on;
    scatter(i,pha_m2{1,i}(ncomp),40,'+b'); hold on;
    %plot([i i],[pha_o{1,i}(ncomp) pha_m1{1,i}(ncomp)],'k'); hold on;
    factor1=(amp_o{1,i}(am)*cosd(pha_o{1,i}(am))-amp_m1{1,i}(am)*cosd(pha_m1{1,i}(am)))^2;
    factor2=(amp_o{1,i}(am)*sind(pha_o{1,i}(am))-amp_m1{1,i}(am)*sind(pha_m1{1,i}(am)))^2;
    complex_e(mm)=sqrt(factor1/2+factor2/2);
  end
end
%complex_e
size(complex_e)
fprintf('%s',['complex error for ',cm,' tide is'])
mean(complex_e(1:mm))
fprintf('%s',['complex error (GOME) for ',cm,' tide is'])
%mean(complex_e(33:39))/2
%mean(complex_e(end-8:end-1))/2
fprintf('%s',['Mean of Square complex error for ',cm,' tide is'])
%RMSD_amp=sqrt(sum((amp_oo-amp_mm1).^2)/mm)
%RMSD_pha=sqrt(sum((pha_oo-pha_mm1).^2)/mm)
RMSD_complex=mean(complex_e(1:mm).^2)
fprintf('%s',['RMSE (GOME) for ',cm,' tide is'])
%RMSD_complex=sqrt(mean(complex_e(33:39).^2)/2)
%RMSD_complex=sqrt(mean(complex_e(end-8:end-1).^2)/2)
%for i=1:nst
%    plot([i i],[pha_o{1,i}(ncomp) pha_m1{1,i}(ncomp)],'k'); hold on;
%end
title([cm 'tidal phase']);
grid on; box on;
xlabel('Station #');
ylabel('Phase');
xlim([0 nst]);
xticks(1:1:nst);
legend('Obs',runo,run);

savefig(strcat(prj_dir,'/Figures/',cm,'pha',num2str(nst),outname,outnameo,'.fig'))
saveas(gcf,strcat(prj_dir,'/Figures/',cm,'pha',num2str(nst),outname,outnameo,'.png'))
figure
subplot(2,1,1)
plot(pha_oo,pha_mm2,'bo');hold on
plot(pha_oo,pha_mm1,'go');
mean(abs(pha_oo-pha_mm2))
mean(abs(pha_oo-pha_mm1))
subplot(2,1,2)
plot(amp_oo,amp_mm2,'bo');hold on
plot(amp_oo,amp_mm1,'go');
mean(abs(amp_oo-amp_mm2))
mean(abs(amp_oo-amp_mm1))
savefig(strcat(prj_dir,'/Figures/',cm,'scatter',num2str(nst),outname,outnameo,'.fig'))
saveas(gcf,strcat(prj_dir,'/Figures/',cm,'scatter',num2str(nst),outname,outnameo,'.png'))
figure
subplot(2,1,1)
for i=1:nst
  if (isempty(pha_o{1,i})==0&&isempty(pha_m1{1,i})==0)
    plot(i,abs(pha_o{1,i}(am)-pha_m2{1,i}(am)),'bo');hold on;
    plot(i,abs(pha_o{1,i}(am)-pha_m1{1,i}(am)),'go');
  end
end
subplot(2,1,2)
for i=1:nst
  if (isempty(amp_o{1,i})==0&&isempty(amp_m1{1,i})==0)
    plot(i,abs(amp_o{1,i}(am)-amp_m2{1,i}(am)),'bo');hold on;
    plot(i,abs(amp_o{1,i}(am)-amp_m1{1,i}(am)),'go');
  end
end
xticks(1:1:nst)
savefig(strcat(prj_dir,'/Figures/',cm,'error',num2str(nst),outname,outnameo,'.fig'))
saveas(gcf,strcat(prj_dir,'/Figures/',cm,'error',num2str(nst),outname,outnameo,'.png'))
aa_amp=corrcoef((amp_oo),(amp_mm2));
cc_amp2=round(aa_amp(1,2),2)
aa_amp=corrcoef((amp_oo),(amp_mm1));
cc_amp1=round(aa_amp(1,2),2)
aa_amp=corrcoef(pha_oo,pha_mm2);
cc_pha2=round(aa_amp(1,2),2)
aa_amp=corrcoef(pha_oo,pha_mm1);
cc_pha1=round(aa_amp(1,2),2)
%figure
%for i=1:nst
%    geoplot(sa_lat(i),sa_lon(i),'ro')
    
%    text(sa_lat(i),sa_lon(i),num2str(i),'FontSize',8)
%    hold on
%end
%}

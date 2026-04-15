function [tday]=coare35_infile_saildrone20hz_clean(ws20hz,time20hz,hgt20hz,at1min,rh1min,sst1min,bp1min,par1min,time1min,sdnum,tc_name,pth)

hgt20hz(hgt20hz<0.5)=0.5;
hgt20hz(isnan(hgt20hz))=3.4;

%average 20 hz wind to 1 min, interpolate 1 min data to new 1 min wind times
j=1;
for i=1:20*60:length(ws20hz)-20*60
    tday(j)=nanmean(time20hz(i:i+20*60-1));
    ws(j)=nanmean(ws20hz(i:i+20*60-1));
    hgt(j)=nanmean(hgt20hz(i:i+20*60-1));    %wind measurement height
    j=j+1;
end
at=nan*ones(1,length(tday));
indg=find(~isnan(at1min)&~isnan(time1min));
indg2=find(~isnan(tday));
at(indg2)=interp1(time1min(indg),at1min(indg),tday(indg2));

rh=nan*ones(1,length(tday));
indg=find(~isnan(rh1min)&~isnan(time1min));
rh(indg2)=interp1(time1min(indg),rh1min(indg),tday(indg2));

sst=nan*ones(1,length(tday));
indg=find(~isnan(sst1min)&~isnan(time1min));
sst(indg2)=interp1(time1min(indg),sst1min(indg),tday(indg2));

bp=nan*ones(1,length(tday));
indg=find(~isnan(bp1min)&~isnan(time1min));
bp(indg2)=interp1(time1min(indg),bp1min(indg),tday(indg2));

pf=nan*ones(1,length(tday));
indg=find(~isnan(par1min)&~isnan(time1min));
pf(indg2)=interp1(time1min(indg),par1min(indg),tday(indg2));


at(isnan(at))=nanmean(at);
rh(isnan(rh))=nanmean(rh);

%solar radiation
%Appl. Sci. 2020, 10, 8007; doi:10.3390/app10228007
%photon flux: micro-mol/s/m^2 (measured by saildrone)
%SWR: J/s/m^2
%PAR = PAR/(4.57 micro-mol/J) (converts to W/m^2 to match SWR)
%SWR = 2*PAR
swr=2*pf/4.57;
swr(isnan(swr))=0;

%eliminate data with bad time values
j=1;
tref=1;
bad_ind=nan;
for i=1:length(tday)-1
    if tday(i+1)<=tday(i) || tday(i+1)<=tday(tref)
        if j==1
            tref=i;
        end
        bad_ind(j)=i+1;
        j=j+1;
    end
end
if ~isempty((find(~isnan(bad_ind))))
    sst(bad_ind)=nan;
    indg=find(~isnan(sst));
    sst=interp1(tday(indg),sst(indg),tday);
end

lwr=370*ones(1,length(sst));
pr=0*ones(1,length(sst));

bp(isnan(bp))=1015;
sst(isnan(sst))=nanmean(sst);

time=nan*ones(1,length(sst));
for i=1:length(sst)
	t=datevec(tday(i));
	time(i)=t(1)*10^8 + t(2)*10^6 + t(3)*10^4 + t(4)*10^2 + t(5)*10^0;
end


%-------------specific humidity
P=1008;
q = (1.0007 + 3.46e-6*P) * 6.1121*exp(17.502*at./(240.97 + at)) .* rh/100;
q = 1000 * 0.62197 * (q./(P-0.378*q));%g/kg
q(find(at<-9000|rh<-9000))=-9999;
%--------------

%     u = relative wind speed (m/s) at height zu(m)
%     zu 
%     t = bulk air temperature (degC) at height zt(m)
%     zt 
%    rh = relative humidity (%) at height zq(m)
%     zq 
%     P = surface air pressure (mb) (default = 1015)
%    ts = water temperature (degC) see jcool below
%    Rs = downward shortwave radiation (W/m^2) (default = 150) 
%    Rl = downward longwave radiation (W/m^2) (default = 370)
%   lat = latitude (default = +45 N)
%    zi = PBL height (m) (default = 600m)
%  rain = rain rate (mm/hr)
%    cp = phase speed of dominant waves (m/s)  
%  sigH =  significant wave height (m)


data=-9999*ones(length(time),15);
data(:,1)=ws';
data(:,2)=hgt';
data(:,3)=at';
data(:,4)=2.3;
data(:,5)=rh';
data(:,6)=2.3;
data(:,7)=bp';
data(:,8)=sst';
data(:,9)=swr';
data(:,10)=lwr';
data(:,11)=20;
data(:,12)=600;
data(:,13)=pr';
data(:,14)=nan;%15;
data(:,15)=nan;%swh';



fid=fopen([pth 'lhf_input_' int2str(sdnum) tc_name '_qcwind.txt'],'w');
fprintf(fid,'%17.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %6.2f %6.2f %6.2f %6.2f %5.2f \r\n',data');
fclose(fid);


fid=fopen([pth 'lhf_input_rh_' int2str(sdnum) tc_name '_qcwind.txt'],'w');
fprintf(fid,'%8.2f \r\n',rh);
fclose(fid);

fid=fopen([pth 'lhf_input_ws_' int2str(sdnum) tc_name '_v_qcwind.txt'],'w');
fprintf(fid,'%8.2f \r\n',ws);
fclose(fid);

fid=fopen([pth 'lhf_input_airt_' int2str(sdnum) tc_name '_qcwind.txt'],'w');
fprintf(fid,'%8.2f \r\n',at);
fclose(fid);

fid=fopen([pth 'lhf_input_sst_' int2str(sdnum) tc_name '_qcwind.txt'],'w');
fprintf(fid,'%8.2f \r\n',sst);
fclose(fid);


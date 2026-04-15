%choose saildrone number (sdnum) and TC name (tc_name) that match names of
%netcdf files
sdnum=10452;
tc_name='IDALIA';

%whether to save final variables to mat-file
save_file=0;

%---------------------

disp([int2str(sdnum),'  ', tc_name]);

%for loading high-res data
pth='~/Desktop/manuscripts/saildrone_drag/manuscript/final_data/';
fname20hz=[int2str(sdnum) tc_name '20hz.nc'];
fname4hz=[int2str(sdnum) tc_name '4hz.nc'];
fname1min=[int2str(sdnum) tc_name '1min.nc'];

%for writing wind height adjustment data
pth_adj='~/Desktop/manuscripts/saildrone_satellite/';
fname1min_adj=['ws_10m_sd' int2str(sdnum) tc_name '_cdadj_zo25mps_qcwind'];

%IBTrACS file for TC tracks
fname_ibtracs='/Users/gregory.foltz/saildrone_atl21/data/IBTrACS.NA.v04r01a.nc';

%for saving final variables to mat-file
pth_save='/Users/gregory.foltz/Desktop/';
fname_save=['output' int2str(sdnum) tc_name ];


if sdnum==1045
    sdyear=2021;
elseif sdnum==1031||sdnum==1040||sdnum==1059||sdnum==1078
    sdyear=2022;
elseif sdnum==10361||sdnum==10362||sdnum==10411||sdnum==10412||sdnum==10451||sdnum==10452||sdnum==10571||sdnum==10572||sdnum==1068||sdnum==1069||sdnum==1083
    sdyear=2023;
else
    sdyear=2024;
end

y1=sdyear;
y2=sdyear;

%---------------------------------------------------------


ave_int=20; %averaging interval for fluxes (minutes)
ave_overlap=10; %overlap of intervals

stream_wind=1;  %instead of using u and v wind: flux calculated based on fluctuations along and/or across mean wind direction (both components if stream_cross=0 and stream_along=0)
stream_cross=0; %only calculate using cross-stream
stream_along=0; %only calculate using along-stream

%---------------------------------------------------------

%1-min
at=ncread([pth fname1min],'at');
rh=ncread([pth fname1min],'rh');
bp=ncread([pth fname1min],'bp');
swh=ncread([pth fname1min],'swh');
sst=ncread([pth fname1min],'sst');
par=ncread([pth fname1min],'par');
datenum_1min=ncread([pth fname1min],'time');
at(at<-9000)=nan;rh(rh<-9000)=nan;bp(bp<-9000)=nan;swh(swh<-9000)=nan;
sst(sst<-9000)=nan;par(par<-9000)=nan;datenum_1min(datenum_1min<-9000)=nan;

indg=find(~isnan(swh));
swh_filled=interp1(indg,swh(indg),[1:length(swh)]);

%20-Hz
ws=ncread([pth fname20hz],'ws');
wdir=ncread([pth fname20hz],'wdir');
u=ncread([pth fname20hz],'wu');
v=ncread([pth fname20hz],'wv');
w=ncread([pth fname20hz],'w');
sdlat=ncread([pth fname20hz],'lat');
sdlon=ncread([pth fname20hz],'lon');
datenum_20hz=ncread([pth fname20hz],'time');
wind_hgt=ncread([pth fname20hz],'wind_hgt');
ws(ws<-9000)=nan;wdir(wdir<-9000)=nan;u(u<-9000)=nan;v(v<-9000)=nan;
w(w<-9000)=nan;sdlat(sdlat<-9000)=nan;sdlon(sdlon<-9000)=nan;
datenum_20hz(datenum_20hz<-9000)=nan;wind_hgt(wind_hgt<-9000)=nan;

%4-Hz
heave=ncread([pth fname4hz],'heave');
datenum_4hz=ncread([pth fname4hz],'time');
heave(heave<-9000)=nan;datenum_4hz(datenum_4hz<-9000)=nan;


%interpolate to times of 20-hz
indg=find(~isnan(heave)&~isnan(datenum_4hz));
indg2=find(~isnan(datenum_20hz));
heave_interp=interp1(datenum_4hz(indg),heave(indg),datenum_20hz(indg2));

%calculate vehicle velocity, vertical acceleration
%for calculation of wave direction
sd_vel_east=nan*ones(1,length(sdlon));
sd_vel_north=nan*ones(1,length(sdlat));
sd_vel_vert=nan*ones(1,length(heave_interp));
sd_acc_vert=nan*ones(1,length(heave_interp));
sd_vel_east(2:end-1)=(sdlon(3:end)-sdlon(1:end-2))*pi/180.*cos(sdlat(2:end-1)*pi/180)*6.67e6;
sd_vel_north(2:end-1)=(sdlat(3:end)-sdlat(1:end-2))*pi/180*6.67e6;
sd_vel_speed=sqrt(sd_vel_east.^2+sd_vel_north.^2);

sd_vel_dir=180/pi*atan(sd_vel_east./sd_vel_north); %direction toward, clockwise from north
sd_vel_dir(sd_vel_east>0 & sd_vel_north<0)=(180/2 + 180/2 + sd_vel_dir(sd_vel_east>0 & sd_vel_north<0)); 
sd_vel_dir(sd_vel_east<0 & sd_vel_north<0)=(180 + sd_vel_dir(sd_vel_east<0 & sd_vel_north<0));
sd_vel_dir(sd_vel_east<0 & sd_vel_north>0)=(3*180/2 + 180/2 + sd_vel_dir(sd_vel_east<0 & sd_vel_north>0));

sd_vel_vert(2:end-1)=heave_interp(3:end)-heave_interp(1:end-2);
sd_acc_vert(2:end-1)=sd_vel_vert(3:end)-sd_vel_vert(1:end-2);


%identify wave troughs, crests
indg=find(~isnan(heave));
heave_4hz_interp=interp1(indg,heave(indg),[1:length(heave)]);

ihv=find(datenum_4hz>=datenum_20hz(1) & datenum_4hz<datenum_20hz(end));
datenum_4hz=datenum_4hz(ihv);
heave_4hz_interp=heave_4hz_interp(ihv);


i=1;
j=1;
k=1;
while i<length(heave_4hz_interp)-10
    while heave_4hz_interp(i+1)>=heave_4hz_interp(i) & ~isnan(heave_4hz_interp(i)) & ~isnan(heave_4hz_interp(i+1)) & i<length(heave)
        i=i+1;
        if i>=(length(heave_4hz_interp)-1)
            break;
        end
    end
    if ~isnan(heave_4hz_interp(i)) & ~isnan(heave_4hz_interp(i+1)) 
        crest_ind(j)=i;
        j=j+1;
    else
        while isnan(heave_4hz_interp(i)) | isnan(heave_4hz_interp(i+1)) & i<length(heave)
            i=i+1;
            if i>=(length(heave_4hz_interp)-1)
                break;
            end
        end
    end
    while heave_4hz_interp(i+1)<=heave_4hz_interp(i) & ~isnan(heave_4hz_interp(i)) & ~isnan(heave_4hz_interp(i+1)) & i<length(heave)
        i=i+1;
        if i>=(length(heave_4hz_interp)-1)
            break;
        end
    end
    if ~isnan(heave_4hz_interp(i)) & ~isnan(heave_4hz_interp(i+1)) 
        trough_ind(k)=i;
        k=k+1;
    else
         while isnan(heave_4hz_interp(i)) | isnan(heave_4hz_interp(i+1)) & i<length(heave)
            i=i+1;
            if i>=(length(heave_4hz_interp)-1)
                break;
            end
        end
    end
end
%calculate individual wave heights, datenum for each height
for i=2:length(crest_ind)
    sub1=heave_4hz_interp(crest_ind(i));
    sub1t=datenum_4hz(crest_ind(i));
    j=min(find(trough_ind>crest_ind(i)));
    sub2=heave_4hz_interp(trough_ind(j));
    sub2t=datenum_4hz(trough_ind(j));
    if ~isempty(sub1)&~isempty(sub2)
        heave_height(i-1)=sub1-sub2;
        heave_height_datenum(i-1)=nanmean([sub1t,sub2t]);
    else
        heave_height(i-1)=nan;
        heave_height_datenum(i-1)=nan;
    end
end
%create 5-minute time series of sig. wave height using heave_height_datenum
j=1;
for i=datenum_4hz(1):(ave_int-ave_overlap)/(60*24):datenum_4hz(end)
    ind=find(heave_height_datenum>i & heave_height_datenum<=(i+5*4*60/(24*60*60*4)));
    if ~isempty(ind)
        height_sorted=sort(heave_height(ind));
        heave_height_5min(j)=nanmean(height_sorted(ceil(length(height_sorted)*2/3):length(height_sorted)));
        heave_height_5min_datenum(j)=nanmean(heave_height_datenum(ind));
    else
        heave_height_5min(j)=nan;
        heave_height_5min_datenum(j)=nan;
    end
    j=j+1;
end





%corrected wind adjusted to 10 m height (from saildrone_wind_adjust2022.m)
tday=coare35_infile_saildrone20hz_clean(ws,datenum_20hz,wind_hgt,at,rh,sst,bp,par,datenum_1min,sdnum,tc_name,pth_adj);
saildrone_wind_adjust20hz_clean(sdnum,tc_name,pth_adj);

load([pth_adj fname1min_adj]);
ws_10m_tilt=ws_10m;                 %with stability correction (U10N)
ws_10m_tilt_nostab=ws_10m_nostab;   %no stability correction (U10)
datenum_1min_adj=datenum_1min;


%20-Hz data
%<u'w'>, <v'w'>

j=1;
jw=1;

u_proc=nan*ones(size(u));
v_proc=u_proc;
w_proc=u_proc;
ws_proc=u_proc;
for i=1:20*60*(ave_int-ave_overlap):length(ws)-20*60*ave_int

    ui=u(i:i+20*60*ave_int);
    vi=v(i:i+20*60*ave_int);
    wi=w(i:i+20*60*ave_int);
    si=ws(i:i+20*60*ave_int);
    wdiri=wdir(i:i+20*60*ave_int);
    
    
    sdlati=sdlat(i:i+20*60*ave_int);
    sdloni=sdlon(i:i+20*60*ave_int);
    
    sdtimei=datenum_20hz(i:i+20*60*ave_int);
    
    
    sd_acc_verti=sd_acc_vert(i:i+20*60*ave_int);
    sd_vel_verti=sd_vel_vert(i:i+20*60*ave_int);
    sd_vel_northi=sd_vel_north(i:i+20*60*ave_int);
    sd_vel_easti=sd_vel_east(i:i+20*60*ave_int);
    heavei=heave_interp(i:i+20*60*ave_int);
    
    sd_vel_diri=sd_vel_dir(i:i+20*60*ave_int);
    sd_vel_speedi=sd_vel_speed(i:i+20*60*ave_int);
    
    
    %calculate wave direction
    %Thomson et al., JAOT, 2018 https://doi.org/10.1175/JTECH-D-17-0091.1
    indg=find(~isnan(sd_acc_verti)&~isnan(sd_vel_northi)&~isnan(sd_vel_easti));
    %band-pass first for swell (6-20 sec)
    [sd_acc_verti,~,~,~,~]=lanczosfilter(lanczosfilter(sd_acc_verti-nanmean(sd_acc_verti),1,1/(20*20),1000,'high'),1,1/(20*6),1000,'low');
    [sd_vel_easti,~,~,~,~]=lanczosfilter(lanczosfilter(sd_vel_easti-nanmean(sd_vel_easti),1,1/(20*20),1000,'high'),1,1/(20*6),1000,'low');
    [sd_vel_northi,~,~,~,~]=lanczosfilter(lanczosfilter(sd_vel_northi-nanmean(sd_vel_northi),1,1/(20*20),1000,'high'),1,1/(20*6),1000,'low');
    if length(indg)>ave_int*20*60 *0.75
        cc1=corrcoef(sd_acc_verti(indg),sd_vel_northi(indg));
        cc2=corrcoef(sd_acc_verti(indg),sd_vel_easti(indg));
        
        corr_vw(j)=cc1(2);
        corr_uw(j)=cc2(2);
        %calculate wave direction for each 2-min segment within 20-min for
        %calculating stdev
        for iww=1:2*60*20:length(sd_acc_verti)-2*60*20
            subv=sd_acc_verti(iww:iww+2*60*20-1);
            sube=sd_vel_easti(iww:iww+2*60*20-1);
            subn=sd_vel_northi(iww:iww+2*60*20-1);
            indg=find(~isnan(subv)&~isnan(subn)&~isnan(sube));
            if ~isempty(indg)
                cc1=corrcoef(subv(indg),subn(indg));
                cc2=corrcoef(subv(indg),sube(indg));
                corr_vw_2min(jw)=cc1(2);
                corr_uw_2min(jw)=cc2(2);
                jw=jw+1;
            else
                corr_vw_2min(jw)=nan;
                corr_uw_2min(jw)=nan;
                jw=jw+1;
            end
        end
    else
        corr_vw(j)=nan;
        corr_uw(j)=nan;
        for iww=1:2*60*20:length(sd_acc_verti)-2*60*20
            corr_vw_2min(jw)=nan;
            corr_uw_2min(jw)=nan;
            jw=jw+1;
        end
    end
    %band-pass for wind waves (<8 sec)
    [sd_acc_verti,~,~,~,~]=lanczosfilter(lanczosfilter(sd_acc_verti-nanmean(sd_acc_verti),1,1/(20*8),1000,'high'),1,1/(20*3),1000,'low');
    [sd_vel_easti,~,~,~,~]=lanczosfilter(lanczosfilter(sd_vel_easti-nanmean(sd_vel_easti),1,1/(20*8),1000,'high'),1,1/(20*3),1000,'low');
    [sd_vel_northi,~,~,~,~]=lanczosfilter(lanczosfilter(sd_vel_northi-nanmean(sd_vel_northi),1,1/(20*8),1000,'high'),1,1/(20*3),1000,'low');
    if length(indg)>ave_int*20*60 *0.75
        cc1=corrcoef(sd_acc_verti(indg),sd_vel_northi(indg));
        cc2=corrcoef(sd_acc_verti(indg),sd_vel_easti(indg));
        corr_vw2(j)=cc1(2);
        corr_uw2(j)=cc2(2);
    else
        corr_vw2(j)=nan;
        corr_uw2(j)=nan;
    end
    
    if length(find(~isnan(si)&~isnan(wi)))>ave_int*60*20*0.4

        if stream_wind  %calculate flux based only on perturbations along mean wind direction
            uave=nanmean(ui);
            vave=nanmean(vi);
            %angle of mean wind (deg. from)
            wdir_ave=atan(vave/uave)*180/pi;
            if uave<0
                wdir_ave=90-wdir_ave;
            elseif uave>=0 && vave>0
                wdir_ave=180+90-wdir_ave;
            elseif uave>=0 && vave<0
                wdir_ave=270-wdir_ave;
            end
            wdir20=atan(reshape(vi,size(ui))./ui)*180/pi;
            wdir20(ui<0)=90-wdir20(ui<0);
            wdir20(ui>=0&reshape(vi,size(ui))>0)=180+90-wdir20(ui>=0&reshape(vi,size(ui))>0);
            wdir20(ui>=0&reshape(vi,size(ui))<0)=270-wdir20(ui>=0&reshape(vi,size(ui))<0);
            dtheta=wdir20-wdir_ave;
            dtheta(dtheta<-180)=dtheta(dtheta<-180)+360;
            dtheta(dtheta>180)=dtheta(dtheta>180)-(360-dtheta(dtheta>180));
            
            %wind components paralell and perpendicular to direction of
            %mean wind
            ws_parai=si.*cos(reshape(dtheta,size(si))*pi/180);
            ws_perpi=si.*sin(reshape(dtheta,size(si))*pi/180);

            fluxu(j)=nanmean((reshape(ws_parai,size(wi))-nanmean(ws_parai)).*(wi-nanmean(wi)));
            if stream_cross %only cross-wind component
                fluxv(j)=nanmean((reshape(ws_perpi,size(wi))-nanmean(ws_perpi)).*(wi-nanmean(wi)));
                fluxu(j)=0;
            elseif stream_along %only along-wind component
                fluxv(j)=0;
            else    %both components
                fluxv(j)=nanmean((reshape(ws_perpi,size(wi))-nanmean(ws_perpi)).*(wi-nanmean(wi)));
            end
            fluxs(j)=nanmean((si-nanmean(si)).*(wi-nanmean(wi)));
            uave_all(j)=uave;                       %20-min mean zonal wind
            vave_all(j)=vave;
            su_ave_all(j)=nanmean(abs(ws_parai));   %20-min ave wind speed along mean wind direction
            sv_ave_all(j)=nanmean(abs(ws_perpi));
            wdir_ave_all(j)=wdir_ave;
            
        else
            fluxu(j)=nanmean((ui-nanmean(ui)).*(wi-nanmean(wi)));
            fluxv(j)=nanmean((vi-nanmean(vi)).*(wi-nanmean(wi)));
            fluxs(j)=nanmean((si-nanmean(si)).*(wi-nanmean(wi)));
            
        end
        
        flux_time(j)=mean(datenum_20hz(i:i+20*60*ave_int)); 
        
    else
        fluxu(j)=nan;
        fluxv(j)=nan;
        fluxs(j)=nan;
        
        flux_time(j)=mean(datenum_20hz(i:i+20*60*ave_int));
        
    end
    
    
    u_mean(j)=nanmean(ui);
    v_mean(j)=nanmean(vi);
    s_mean(j)=nanmean(si);
    wdir_mean(j)=nanmean(wdiri);

    wi_std(j)=nanstd(wi);
    
    sdlat_mean(j)=nanmean(sdlati);
    sdlon_mean(j)=nanmean(sdloni);
    sdtime_mean(j)=nanmean(sdtimei);
    

    if wi_std(j)<0.05
        fluxs(j)=nan;   %obviously bad data with no vertical wind variations
        fluxu(j)=nan;
        fluxv(j)=nan;
    end
    
    j=j+1;
end


%TC quandrant and distance, angle from TC-------------------
[tc_quad,tc_sd_dist,tc_sd_angle,tc_rmw]=tc_quadrant_clean(sdlat_mean,sdlon_mean,sdtime_mean, tc_name,sdyear,fname_ibtracs);
tc_sd_angle(tc_sd_angle<0)=tc_sd_angle(tc_sd_angle<0)+360;
%-------------------------------------------


%adjust wave angle convention so it's direction from, clockwise from north
%(same as wind)
sd_wave_angle=180/pi*atan(corr_vw./corr_uw);
sd_wave_angle(corr_vw<0 & corr_uw<0)=sd_wave_angle(corr_vw<0 & corr_uw<0)+180;
sd_wave_angle(corr_vw>0 & corr_uw<0)=sd_wave_angle(corr_vw>0 & corr_uw<0)+180;
sd_wave_angle(sd_wave_angle<0)=sd_wave_angle(sd_wave_angle<0)+360;%+180;
%convert to direction to, starting from north and going clockwise
sd_wave_angle=90-sd_wave_angle;
sd_wave_angle(sd_wave_angle<0)=sd_wave_angle(sd_wave_angle<0)+360;
%convert to direction from, starting from north
sd_wave_angle(sd_wave_angle<180)=sd_wave_angle(sd_wave_angle<180)+360;
sd_wave_angle(sd_wave_angle>180)=sd_wave_angle(sd_wave_angle>180)-180;

%2 min wave direction
sd_wave_angle_2min=180/pi*atan(corr_vw_2min./corr_uw_2min);
sd_wave_angle_2min(corr_vw_2min<0 & corr_uw_2min<0)=sd_wave_angle_2min(corr_vw_2min<0 & corr_uw_2min<0)+180;
sd_wave_angle_2min(corr_vw_2min>0 & corr_uw_2min<0)=sd_wave_angle_2min(corr_vw_2min>0 & corr_uw_2min<0)+180;
sd_wave_angle_2min(sd_wave_angle_2min<0)=sd_wave_angle_2min(sd_wave_angle_2min<0)+360;%+180;
%convert to direction to, starting from north and going clockwise
sd_wave_angle_2min=90-sd_wave_angle_2min;
sd_wave_angle_2min(sd_wave_angle_2min<0)=sd_wave_angle_2min(sd_wave_angle_2min<0)+360;
%convert to direction from, starting from north
sd_wave_angle_2min(sd_wave_angle_2min<180)=sd_wave_angle_2min(sd_wave_angle_2min<180)+360;
sd_wave_angle_2min(sd_wave_angle_2min>180)=sd_wave_angle_2min(sd_wave_angle_2min>180)-180;
%mean and stdev of 2-min wave direction every 20 min
wwx=1;
for wd=1:20/2/2:length(sd_wave_angle_2min)-9
    wdsub=sd_wave_angle_2min(wd:wd+9);
    wdmax=max(wdsub)-min(wdsub);
    if ~isnan(wdmax) && wdmax<=180
        sd_wave_angle_std(wwx)=nanstd(wdsub);
        sd_wave_angle_mean_from_2min(wwx)=nanmean(wdsub);
    else
        sd_wave_angle_std(wwx)=nan;
        sd_wave_angle_mean_from_2min(wwx)=nan;
    end
    wwx=wwx+1;
end

sd_windwave_angle=180/pi*atan(corr_vw2./corr_uw2);
sd_windwave_angle(corr_vw2<0 & corr_uw2<0)=sd_windwave_angle(corr_vw2<0 & corr_uw2<0)+180;
sd_windwave_angle(corr_vw2>0 & corr_uw2<0)=sd_windwave_angle(corr_vw2>0 & corr_uw2<0)+180;
sd_windwave_angle(sd_windwave_angle<0)=sd_windwave_angle(sd_windwave_angle<0)+360;%+180;
%convert to direction to, starting from north and going clockwise
sd_windwave_angle=90-sd_windwave_angle;
sd_windwave_angle(sd_windwave_angle<0)=sd_windwave_angle(sd_windwave_angle<0)+360;
%convert to direction from, starting from north
sd_windwave_angle(sd_windwave_angle<180)=sd_windwave_angle(sd_windwave_angle<180)+360;
sd_windwave_angle(sd_windwave_angle>180)=sd_windwave_angle(sd_windwave_angle>180)-180;


ustar2=sqrt(fluxu.^2+reshape(fluxv,size(fluxu)).^2);


%average adjusted wind, wave period to ave_int
j=1;
for i=1:(ave_int-ave_overlap):length(ws_10m_tilt)-ave_int+1    
    ws_10m_tilt_interp(j)=nanmean(ws_10m_tilt(i:i+ave_int-1));
    ws_10m_tilt_nostab_interp(j)=nanmean(ws_10m_tilt_nostab(i:i+ave_int-1));
    datenum_1min_adj_interp(j)=nanmean(tday(i:i+ave_int-1));
   j=j+1;
end

%interpolate to flux times
indg2=find(~isnan(flux_time));
indg=find(~isnan(heave_height_5min)&~isnan(heave_height_5min_datenum));
heave_height_5min_interp=interp1(heave_height_5min_datenum(indg),heave_height_5min(indg),flux_time(indg2));


%times of momentum flux, wind, wind direction
datenum_flux=flux_time(indg2);

%interpolate 1-min data to times of fluxes
indg=find(~isnan(at)&~isnan(datenum_1min));
indg2=find(~isnan(flux_time));
at_interp=interp1(datenum_1min(indg),at(indg),flux_time(indg2));
indg=find(~isnan(rh)&~isnan(datenum_1min));
rh_interp=interp1(datenum_1min(indg),rh(indg),flux_time(indg2));
indg=find(~isnan(bp)&~isnan(datenum_1min));
bp_interp=interp1(datenum_1min(indg),bp(indg),flux_time(indg2));
indg=find(~isnan(swh_filled')&~isnan(datenum_1min));
swh_interp=interp1(datenum_1min(indg),swh_filled(indg)',flux_time(indg2));


%air density
Rd=287.058;
Rv=461.495;
p1=6.1078 * 10.^(7.5*at_interp./(at_interp+237.3));
pv=100*p1.*rh_interp/100;
pd=(100*bp_interp - pv);
air_dens=(pd./(Rd*(at_interp+273.15))) + (pv./(Rv*(at_interp+273.15)));
air_dens(isnan(air_dens))=nanmean(air_dens);

%adjust wind to 10 m height using observed ustar
z0=3.4./exp(s_mean*0.4./sqrt(ustar2));
u10_from_ustar=sqrt(ustar2)./0.4.*log(10./z0);
u10_from_ustar(z0==0|s_mean==0)=nan;
%stability adjustment from COARE3.5 to make it U10N instead of U10
ws_10m_tilt_nostab_interp(ws_10m_tilt_nostab_interp<1)=ws_10m_tilt_interp(ws_10m_tilt_nostab_interp<1);
u10_from_ustar=u10_from_ustar.*reshape(ws_10m_tilt_interp,size(u10_from_ustar))./reshape(ws_10m_tilt_nostab_interp,size(u10_from_ustar));


%drag coefficient
cd=ustar2./u10_from_ustar./u10_from_ustar;
ws_10m_tilt_interp=u10_from_ustar;  %for writing to file

taux=air_dens.*fluxu;
tauy=air_dens.*fluxv;

cd(cd>0.1)=nan;
cdu=fluxu./u_mean./s_mean;
cdv=fluxv./v_mean./s_mean;

turb_flux=air_dens.*ustar2;

if save_file
    
    save([pth_save fname_save '.mat'],'ustar2','datenum_flux','cd','turb_flux',...
        'ws_10m_tilt_interp','u_mean','v_mean','taux','tauy','sd_wave_angle','sd_windwave_angle','heave_height_5min_interp','wdir_mean','tc_quad',...
        'tc_sd_dist','tc_sd_angle','tc_rmw');
    
end



function [tc_quad,sd_milton_dist,angle_diff,milton_rmw_interp]=tc_quadrant_clean(lat1083,lon1083,time1083, tc_name,tc_year,pth)
%Input:     saildrone lat, lon, time
%           TC name, year
%Output:    TC quadrant that saildrone was in at each time
%           (1=front-right, 2=rear-right, 3=rear-left, 4=front-left
%           Saildrone distance from TC center
%           Angle saildrone position makes with TC track, clockwise from
%           north
%           TC radius of max wind

tc_names=ncread(pth,'name');
tc_years=ncread(pth,'season');
tc_lon=ncread(pth,'lon');
tc_lat=ncread(pth,'lat');
tc_wnd=ncread(pth,'usa_wind');
tc_rmw=ncread(pth,'usa_rmw');
tc_time=ncread(pth,'time')+datenum(1858,11,17);

cc=1;
for isd=1:length(tc_names(1,:))
    if strcmp(tc_names(1:length(tc_name),isd)',tc_name) && tc_years(isd)==tc_year
        milton_ind(cc)=isd;
        cc=cc+1;
    end
end
milton_lat=tc_lat(:,milton_ind);
milton_lon=tc_lon(:,milton_ind);
milton_wnd=tc_wnd(:,milton_ind)/1.944;  %given in kt, convert to m/s
milton_rmw=tc_rmw(:,milton_ind)*1.852;  %given in nm, convert to km
milton_time=tc_time(:,milton_ind);

%interpolate TC times and locations to saildrone times and locations
indg=find(~isnan(milton_time));
milton_lon_interp=interp1(milton_time(indg),milton_lon(indg),time1083);
milton_lat_interp=interp1(milton_time(indg),milton_lat(indg),time1083);
%interpolate TC RMW to saildrone times, locations
milton_rmw_interp=interp1(milton_time(indg),milton_rmw(indg),time1083);
milton_rmw_interp=milton_rmw_interp';

%angle of TC track at each SD time, clockwise from north
milton_dx=nan*ones(1,length(milton_lon_interp));
milton_dy=nan*ones(1,length(milton_lon_interp));

milton_dx(2:end-1)=(milton_lon_interp(3:end)-milton_lon_interp(1:end-2))*pi/180.*cos(milton_lat_interp(2:end-1)*pi/180)*6.67e6;
milton_dy(2:end-1)=(milton_lat_interp(3:end)-milton_lat_interp(1:end-2))*pi/180*6.67e6;

milton_angle=atan(milton_dx./milton_dy)*180/pi;
milton_angle(milton_dx>0 & milton_dy<0)=180+milton_angle(milton_dx>0 & milton_dy<0);
milton_angle(milton_dx<0 & milton_dy<0)=180+milton_angle(milton_dx<0 & milton_dy<0);
milton_angle(milton_dx<0 & milton_dy>=0)=360+milton_angle(milton_dx<0 & milton_dy>=0);

%angle of saildrone location relative to TC location, clockwise from north
sd_milton_dx=nan*ones(1,length(lon1083));
sd_milton_dy=nan*ones(1,length(lon1083));

sd_milton_dx=(lon1083-milton_lon_interp)*pi/180.*cos(milton_lat_interp*pi/180)*6.67e6;
sd_milton_dy=(lat1083-milton_lat_interp)*pi/180*6.67e6;
sd_milton_dist=sqrt(sd_milton_dx.^2+sd_milton_dy.^2);

sd_milton_angle=atan(sd_milton_dx./sd_milton_dy)*180/pi;
sd_milton_angle(sd_milton_dx>0 & sd_milton_dy<0)=180+sd_milton_angle(sd_milton_dx>0 & sd_milton_dy<0);
sd_milton_angle(sd_milton_dx<0 & sd_milton_dy<0)=180+sd_milton_angle(sd_milton_dx<0 & sd_milton_dy<0);
sd_milton_angle(sd_milton_dx<0 & sd_milton_dy>=0)=360+sd_milton_angle(sd_milton_dx<0 & sd_milton_dy>=0);

%calculate TC quadrant
tc_quad=nan*ones(1,length(milton_angle));
angle_diff=sd_milton_angle-milton_angle;
tc_quad(angle_diff>=0 & angle_diff<90)=1;
tc_quad(angle_diff>=90 & angle_diff<180)=2;
tc_quad(angle_diff>=180 & angle_diff<270)=3;
tc_quad(angle_diff>=270 & angle_diff<360)=4;

tc_quad(angle_diff<=0 & angle_diff>-90)=4;
tc_quad(angle_diff<=-90 & angle_diff>-180)=3;
tc_quad(angle_diff<=-180 & angle_diff>-270)=2;
tc_quad(angle_diff<=-270 & angle_diff>=-360)=1;




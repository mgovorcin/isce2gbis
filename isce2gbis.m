function GBISinput = ISCE2GBIS(isce_dir, output_dir, seed_ref, crop)
%Function to generate GBIS INSAR input file from ISCE generated interferogram
%
% Usage: GBISinput = ISCE2GBIS(isce_dir, seed_ref, crop)
% isce_dir  = '/dir/' - define the directory with filt_topophase.unw.geo
% file
%
% output_dir = 'dir' -define the output dir for GBIS insar input file
%
% seed_ref = [Lon,Lat] - define the unwrapping seed point (default left
% upper corner), move it closer to deformation area but not on the area
% affected by deformation (OPTIONAL)
%
% crop= [lon_min, lon_max, lat_min, lat_max] - Define the window for
% cropping (OPTIONAL)
% =========================================================================
% This function is created by: Marin Govorcin (22.08.2019)
% Address: Faculty of Geodesy, University of Zagreb
% Email: mgovorcin@geof.hr
% TO DO LIST: Script works with StripmapApp.py- update it to work with
% TopsApp.py and InsarApp.py (code lines: 99, 102, 138)

tic
%% EDIT THE ENVIRONMENTAL PARAMETERS !!!
 %define the ISCE installation folder, example
 %/dir/isce/isce-2.2.0/install/isce
 
 isce_install_dir = '/home/mgovorcin/Installation_dir/isce/isce-2.2.0/install/isce/';
 
 %OPTIONAL - Define the directory of Water shape file for cliping the interferogram with coastal
 %boundaries, usefull for coastal areas or islands
 %can be downloaded at https://osmdata.openstreetmap.de/data/water-polygons.html
 
 water_dir='/home/mgovorcin/Working_dir/SCRIPTS/ISCE2GBIS/water-polygons-split-4326';
 
 plot_flag =1; %Debug /plot the GBIS input
 
%% Setup the input information
if ~exist('isce_dir')
    fprintf('ERROR - Define the ISCE input directory\n')
    return
end

if ~exist('output_dir')
    fprintf('ERROR - Define the ISCE2GBIS output directory\n')
    return
end

%in case of StripmapAp.py where filt_topophase.unw.geo is in folder
%"interferogram"
if strfind(isce_dir,'interferogram')~=0
    isce_dir = regexprep(isce_dir,'interferogram','','ignorecase');
    stripmap_flag = 1;
end

if ~exist('isce_install_dir')
    fprintf('ERROR - Define the ISCE installation directory\n /dir/isce/isce-2.2.0/install\n')
    return
end

if exist('water_dir')
    water_flag=1;
else
    water_flag=0;
end

if exist('seed_ref')
    seed_flag=1;
else
    seed_flag=0;
end

if exist('crop')
    crop_flag=1;
else
    crop_flag=0;
end

 %SET NAME OF OUTPUT
 %Define the name of ISCE processing folder as SAT_DIRECTION_MASTER_SLAVE;
 %example ERS_DESC_199612237_20010907
 filepath =  fileparts(strcat(isce_dir));
 [filepath, output_name] = fileparts(filepath);


%% EXPORT AND CONVERSION
%Define ISCE export to VRT command
if stripmap_flag==1
isce_export = sprintf('cd %s; python %s/applications/isce2gis.py vrt -i ./interferogram/filt_topophase.unw.geo',isce_dir,isce_install_dir);
isce_export_mask = sprintf('cd %s; python %s/applications/isce2gis.py vrt -i ./interferogram/filt_topophase.unw.conncomp.geo',isce_dir,isce_install_dir);
isce_base_dir =isce_dir;
isce_dir = sprintf('%sinterferogram',isce_dir);

%get wavelength
[ans,wavelength] = system(sprintf('grep -Eo wavelength.{20} %s/stripmapProc.xml | awk ''NR==1 {print $1}'' | sed ''s/wavelength>//g''',isce_base_dir));
wavelength = str2num(strcat(wavelength));

else
isce_export = sprintf('cd %s; python %s/applications/isce2gis.py vrt -i ./filt_topophase.unw.geo',isce_dir,isce_install_dir);
isce_base_dir =isce_dir;

%[ans,wavelength] = system(sprintf('grep -Eo wavelength.{20} %s/stripmapProc.xml | awk ''NR==1 {print $1}'' | sed ''s/wavelength>//g''',isce_base_dir));
%wavelength = str2num(strcat(wavelength));
end

%Export to VRT file
system(isce_export);
system(isce_export_mask);

%Get VRT info
[ans,x_min] = system(sprintf('gdalinfo %s/filt_topophase.unw.geo.vrt | grep -Eo "Lower Left.{15}" | awk ''{print $4}''',isce_dir));
[ans,y_min] = system(sprintf('gdalinfo %s/filt_topophase.unw.geo.vrt | grep -Eo "Lower Left.{28}" | awk ''{print $5}''',isce_dir));
[ans,x_max] = system(sprintf('gdalinfo %s/filt_topophase.unw.geo.vrt | grep -Eo "Upper Right.{14}" | awk ''{print $4}''',isce_dir));
[ans,y_max] = system(sprintf('gdalinfo %s/filt_topophase.unw.geo.vrt | grep -Eo "Upper Right.{27}" | awk ''{print $5}''',isce_dir));

if water_flag == 1
%crop water shape file for your area of interest
water_crop_area=sprintf('ogr2ogr -overwrite -wrapdateline -clipdst %s %s %s %s %s/water.shp %s/water_polygons.shp',strcat(x_min),strcat(y_min),strcat(x_max),strcat(y_max), isce_dir, water_dir);
system(water_crop_area);

crop_water=sprintf('gdal_rasterize -b 1 -b 2 -burn 0 -burn 0 %s/water.shp %s/filt_topophase.unw.geo.vrt',isce_dir,isce_dir);
crop_water_mask=sprintf('gdal_rasterize -b 1 -burn 0  %s/water.shp %s/filt_topophase.unw.conncomp.geo.vrt',isce_dir,isce_dir);

system(crop_water);
system(crop_water_mask);
end

%prepare isce phase output (Lon Lat Phase)
vrt2xyz = sprintf('gdal_translate -of Gtiff -b 2 %s/filt_topophase.unw.geo.vrt %s/unwrap_phase.tif; gdal_translate -of XYZ -b 2 %s/filt_topophase.unw.geo.vrt %s/gbis_lonlatphase.xyz',isce_dir,isce_dir,isce_dir,isce_dir);
system(vrt2xyz);
mask2xyz= sprintf('gdal_translate -of XYZ -b 1 %s/filt_topophase.unw.conncomp.geo.vrt %s/gbis_unw_conncomp.xyz',isce_dir,isce_dir);
system(mask2xyz);

%prepare isce Inc and heading output (Inc Heading)
if stripmap_flag==1
inc_head = sprintf('gdalwarp -overwrite  -multi -te %s %s %s %s %s/geometry/los.rdr.geo.vrt %s/los.rdr.geo.clip.vrt; gdal_translate -of XYZ -b 1 %s/los.rdr.geo.clip.vrt %s/gbis_inc.xyz; gdal_translate -of XYZ -b 2 %s/los.rdr.geo.clip.vrt %s/gbis_heading.xyz',strcat(x_min),strcat(y_min),strcat(x_max),strcat(y_max), isce_base_dir, isce_dir, isce_dir,isce_dir,isce_dir,isce_dir);
system(inc_head);
else
inc_head = sprintf('gdalwarp -overwrite  -multi -te %s %s %s %s %s/los.rdr.geo.vrt %s/los.rdr.geo.clip.vrt; gdal_translate -of XYZ -b 1 %s/los.rdr.geo.clip.vrt %s/gbis_inc.xyz; gdal_translate -of XYZ -b 2 %s/los.rdr.geo.clip.vrt %s/gbis_heading.xyz',strcat(x_min),strcat(y_min),strcat(x_max),strcat(y_max), isce_base_dir, isce_dir, isce_dir,isce_dir,isce_dir,isce_dir);
system(inc_head);
end 
%% LOAD DATA 
INC_INPUT = importdata([isce_dir '/gbis_inc.xyz']);
PHASE_INPUT = importdata([isce_dir '/gbis_lonlatphase.xyz']);
HEADING_INPUT = importdata([isce_dir '/gbis_heading.xyz']);
CONN_PARTS = importdata([isce_dir '/gbis_unw_conncomp.xyz']);
GBIS_NAME = output_name;

%% PREPARE GBIS INPUTFILE
%define satellite heading from pixel heading
HEADING_INPUT=HEADING_INPUT(:,3)*-1+90;
c=find(HEADING_INPUT(:,1)== 90);
HEADING_INPUT(c,:)=0;

%CREATE MATRIX
GBIS=[PHASE_INPUT INC_INPUT(:,3) HEADING_INPUT];

%mask phase with connected parts
GBIS(CONN_PARTS(:,3)==0,:)=[];

%delete zero phase values
zero=find(GBIS(:,3)==0);
GBIS(zero,:) =[];

%correct phase for seed point
if seed_flag ==1
ref=find(GBIS(:,1)>seed_ref(:,1)-0.0005 & GBIS(:,1)<seed_ref(:,1)+0.0005 & GBIS(:,2)>seed_ref(:,2)-0.0005 & GBIS(:,2)<seed_ref(:,2)+0.0005);
seed=GBIS(ref,3);
PHASE=GBIS(:,3)-mean(seed);

GBIS(:,3)=PHASE;

%delete zero phase values
zero2=find(GBIS(:,4)==0);
GBIS(zero2,:) = [];
end

%crop area
if crop_flag==1
crop_lim=find(GBIS(:,1)>crop(:,1) & GBIS(:,1)<crop(:,2) & GBIS(:,2)>crop(:,3) & GBIS(:,2)<crop(:,4));
GBIS=GBIS(crop_lim,:);
crop_lim=find(CONN_PARTS(:,1)>crop(:,1) & CONN_PARTS(:,1)<crop(:,2) & CONN_PARTS(:,2)>crop(:,3) & CONN_PARTS(:,2)<crop(:,4));
CONN_PARTS=CONN_PARTS(crop_lim,:);
end

for i=1:2
    if i==1
%define GBIS input variables
l=length(GBIS);
%Heading=ones(l,1)*HEADING_INPUT;
Heading=double(GBIS(:,5));
Lon=double(GBIS(:,1));
Lat=double(GBIS(:,2));
Phase=double(GBIS(:,3));
Inc=double(GBIS(:,4));
Conn_part=double(CONN_PARTS(CONN_PARTS(:,3)~=0,3));

if ~exist('output_name')
    output_name = 'insar_gbis_input'
end

mat_name=strcat(output_dir,'/',output_name,'.mat');
save(mat_name,'Lon','Lat','Phase','Inc','Heading','Conn_part')
fprintf('Saving GBIS INSAR inputfile: %s \n',mat_name);
%rm Heading Lon Lat Phase Inc mat_name

    elseif i==2
%decimate GBIS input variables until less then 1 million points are
%achieved
for j=1:100
    if length(GBIS)>1000000
        GBIS=GBIS(1:2:end,:);
    elseif length(GBIS)<1000000
         break
    else
    end
end

l=length(GBIS);
%Heading=ones(l,1)*HEADING_INPUT;
Heading=double(GBIS(:,5));
Lon=double(GBIS(:,1));
Lat=double(GBIS(:,2));
Phase=double(GBIS(:,3));
Inc=double(GBIS(:,4));
Conn_part=double(CONN_PARTS(CONN_PARTS(:,3)~=0,3));
mat_name=strcat(output_dir,'/',output_name,'_sparse.mat');
save(mat_name,'Lon','Lat','Phase','Inc','Heading','Conn_part')
fprintf('Saving decimated GBIS INSAR inputfile: %s \n',mat_name);
    end  
end
 if plot_flag == 1
 % Plot wrapped dataset
    figure
    convertedPhase = (Phase / (4*pi)) * wavelength;   % Convert phase from radians to m
    los = single(-convertedPhase);                    % Convert to Line-of-sigth displacement in m
    scatter(Lon, Lat, [], mod(los,wavelength/2),'.');
    hold on
    a=plot(seed_ref(1,1),seed_ref(1,2),'ks','MarkerFaceColor','black','MarkerSize',10);
    hold off
    legend(a,'Reference Point')
    cmapSeismo = colormap_cpt('GMT_seis.cpt', 100);
    colormap(cmapSeismo)
    caxis([0 wavelength/2])
    axis tight
    axis equal
    title('Wrapped Interferogram')
    xlabel('Longitude (degrees)')
    ylabel('Latitude (degrees)')
    colorbar
    fig_name=strcat(output_dir,'/',output_name,'_isce_wrapp_ifg.png');   
    saveas(gcf,fig_name)
    figure
    scatter(CONN_PARTS(CONN_PARTS(:,3)~=0,1),CONN_PARTS(CONN_PARTS(:,3)~=0,2),[],CONN_PARTS(CONN_PARTS(:,3)~=0,3),'.')
    colorbar
    axis tight
    axis equal
    title('Unwrapped Interferogram components')
    xlabel('Longitude (degrees)')
    ylabel('Latitude (degrees)')
    fig_name=strcat(output_dir,'/',output_name,'_isce_unw_comm.png');
    saveas(gcf,fig_name)
 end         
 toc 








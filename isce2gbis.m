function GBISinput = ISCE2GBIS(isce_dir, output_dir, seed_ref, crop)
%Function to generate GBIS INSAR input file from ISCE generated interferogram
%
% Usage: GBISinput = ISCE2GBIS(isce_dir, seed_ref, crop)
% isce_dir  = '/dir/' - define the ISCE processing directory
% 
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
 
 
 stripmap_flag =0;
 tops_flag=0;
 
%% Setup the input information
% Check the inputs
if ~exist('isce_dir')
    fprintf('ERROR - Define the ISCE input directory\n')
    return
end

if ~exist('output_dir')
    fprintf('ERROR - Define the ISCE2GBIS output directory\n')
    return
else
    %Check if output_dir exists
    if ~exist(output_dir,'dir')
        % Create if does not exist
        mkdir(output_dir)
    end
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


%% Find the ISCE interferogram

%Look for the filt_*.unw.geo in the input folder, in matlab version after
%2015 this can be done with 'dir **/*unw.geo'

[out,ifg_dir] = system(['ls ',isce_dir,'/**/*unw.geo']);

ifg_dir = strsplit(ifg_dir);

for n = 1 : (length(ifg_dir)-1) %Skip the last empty cell
    ifg_dir_new{n} = ifg_dir{n};
end


%Check if unwrapping is done in 2 stages and use only unwrap_2_stage

    for n = 1:length(ifg_dir_new)
        ck = strfind(ifg_dir_new{n},'2stage');
        if ~isempty(ck)
            ifg = ifg_dir_new{n};
            ifg_dir = fileparts(ifg_dir_new{n});
            stage2_flag = 1;
        else
            ifg = ifg_dir_new{n};
            ifg_dir = fileparts(ifg_dir_new{n});
        end
    end
    
    clear ifg_dir_new ck n

    
%Look for the connected part of unw filt_*.unw,conncomp.geo in the input folder, in matlab version after
%2015 this can be done with 'dir **/*unw.conncomp.geo'

if ~isempty(dir([ifg_dir,'/*.unw.conncomp.geo']))
   [out,ifg_cp] = system(['ls ',isce_dir,'/**/*unw.conncomp.geo']);
   cp_flag = 1;
end

%Look for the los.rdr.geo in the input folder, in matlab version after
%2015 this can be done with 'dir **/los.rdr.geo'
[out,los_file] = system(['ls ',isce_dir,'/**/los.rdr.geo']);
    los_new=strsplit(los_file);
    los_file = los_new{1};
    clear los_new
    
if strfind(los_file,'No such file')
    fprintf('los.rdr.geo file is missing, put it on geocoding list and reprocess\n')
    return
end
%Get Proc log file

if isempty(dir([isce_dir,'/*Proc.xml']))
    fprintf('Processing was done with Stack, define wavelength manually\n')
    w_flag =1;
else
    proc = dir([isce_dir,'/*Proc.xml']);
    proc_file = [isce_dir,'/',proc.name];
end

%%

 %SET NAME OF OUTPUT
 %Define the name of ISCE processing folder as SAT_DIRECTION_MASTER_SLAVE;
 %example ERS_DESC_199612237_20010907
 [filepath, output_name] = fileparts(strcat(isce_dir));
 display(output_name)

 %% EXPORT AND CONVERSION
isce_export = sprintf('cd %s; python %s/applications/isce2gis.py vrt -i %s',isce_dir,isce_install_dir,ifg);
system(isce_export);
%create conncomp vrt file 
if exist('cp_flag')
    isce_export_mask = sprintf('cd %s; python %s/applications/isce2gis.py vrt -i %s',isce_dir,isce_install_dir,ifg_cp);
    system(isce_export_mask);
end 

%get wavelength
%When interfergoram is obtain though stack processing
if exist('w_flag')
    wavelength=input('Define wavelength [m]\n');
else
    %check if the processing is done with tops
    if ~isempty(strfind(proc_file,'tops'))
        wavelength = 0.05546576; %for Sentinel-1, wavelenth is not stored in topsProc.xml
    else
       [ans,wavelength] = system(sprintf('grep -Eo wavelength.{20} %s | awk ''NR==1 {print $1}'' | sed ''s/wavelength>//g''',proc_file));
       wavelength = str2num(strcat(wavelength));
    end
end

%%CROP IFG BEFORE OTHER PROCESSING

%% EXPORT AND CONVERSION
%Define ISCE export to VRT command
% if stripmap_flag==1
% isce_export = sprintf('cd %s; python %s/applications/isce2gis.py vrt -i ./interferogram/filt_topophase.unw.geo',isce_dir,isce_install_dir);
% isce_export_mask = sprintf('cd %s; python %s/applications/isce2gis.py vrt -i ./interferogram/filt_topophase.unw.conncomp.geo',isce_dir,isce_install_dir);
% isce_base_dir =isce_dir;
% isce_dir = sprintf('%sinterferogram',isce_dir);
% 
% %get wavelength
% [ans,wavelength] = system(sprintf('grep -Eo wavelength.{20} %s/stripmapProc.xml | awk ''NR==1 {print $1}'' | sed ''s/wavelength>//g''',isce_base_dir));
% wavelength = str2num(strcat(wavelength));
% 
% elseif tops_flag==1
% isce_export = sprintf('cd %s; python %s/applications/isce2gis.py vrt -i ./merged/filt_topophase.unw.geo',isce_dir,isce_install_dir);
% isce_export_mask = sprintf('cd %s; python %s/applications/isce2gis.py vrt -i ./merged/filt_topophase.unw.conncomp.geo',isce_dir,isce_install_dir);
% isce_base_dir =isce_dir;
% isce_dir = sprintf('%sinterferogram',isce_dir);
% 
% %get wavelength
% [ans,wavelength] = system(sprintf('grep -Eo wavelength.{20} %s/topsProc.xml | awk ''NR==1 {print $1}'' | sed ''s/wavelength>//g''',isce_base_dir));
% wavelength = str2num(strcat(wavelength));
% 
% else
% isce_export = sprintf('cd %s; python %s/applications/isce2gis.py vrt -i ./filt_topophase.unw.geo',isce_dir,isce_install_dir);
% isce_base_dir =isce_dir;
% 
% %[ans,wavelength] = system(sprintf('grep -Eo wavelength.{20} %s/stripmapProc.xml | awk ''NR==1 {print $1}'' | sed ''s/wavelength>//g''',isce_base_dir));
% %wavelength = str2num(strcat(wavelength));
% end

% %Export to VRT file
% system(isce_export);
% system(isce_export_mask);

%Get VRT info
[ans,x_min] = system(sprintf('gdalinfo %s.vrt | grep -Eo "Lower Left.{15}" | awk ''{print $4}''',ifg));
[ans,y_min] = system(sprintf('gdalinfo %s.vrt | grep -Eo "Lower Left.{28}" | awk ''{print $5}''',ifg));
[ans,x_max] = system(sprintf('gdalinfo %s.vrt | grep -Eo "Upper Right.{14}" | awk ''{print $4}''',ifg));
[ans,y_max] = system(sprintf('gdalinfo %s.vrt | grep -Eo "Upper Right.{27}" | awk ''{print $5}''',ifg));

% Create temp folder to store files through the code processing
temp_dir = [isce_dir,'/temp'];
if ~exist(temp_dir,'dir')
    mkdir(temp_dir);
end


if water_flag == 1
%crop water shape file for your area of interest
water_crop_area=sprintf('ogr2ogr -overwrite -wrapdateline -clipdst %s %s %s %s %s/water.shp %s/water_polygons.shp',strcat(x_min),strcat(y_min),strcat(x_max),strcat(y_max), temp_dir, water_dir);
fprintf('Cropping water shp file\n')
system(water_crop_area);

fprintf('Use water shp to mask water areas in interferogram\n')
crop_water=sprintf('gdal_rasterize -b 1 -b 2 -burn 0 -burn 0 %s/water.shp %s.vrt',temp_dir,ifg);
system(crop_water);

if exist('cp_flag')
    fprintf('Use water shp to mask water areas in interferogram conncomp file\n')
    crop_water_mask=sprintf('gdal_rasterize -b 1 -burn 0  %s/water.shp %s/filt_topophase.vrt',temp_dir,ifg_cp);
    system(crop_water_mask);
end
end

%prepare isce phase output (Lon Lat Phase)
fprintf('Prepare isce phase output (Lon Lat Phase)\n')
    vrt2xyz = sprintf('gdal_translate -of Gtiff -b 2 %s.vrt %s/unwrap_phase.tif; gdal_translate -of XYZ -b 2 %s.vrt %s/gbis_lonlatphase.xyz',ifg,temp_dir,ifg,temp_dir);
    system(vrt2xyz);
if exist('cp_flag')
    fprintf('Prepare isce concomp output\n')
    mask2xyz= sprintf('gdal_translate -of XYZ -b 1 %s.vrt %s/gbis_unw_conncomp.xyz',ifg_cp,temp_dir);
    system(mask2xyz);
end


%prepare isce Inc and heading output (Inc Heading)
% if stripmap_flag==1
fprintf('Prepare isce inc and heading output (Lon Lat Inc Heading)\n')
inc_head = sprintf('gdalwarp -overwrite  -multi -te %s %s %s %s %s.vrt %s/los.rdr.geo.clip.vrt; gdal_translate -of XYZ -b 1 %s/los.rdr.geo.clip.vrt %s/gbis_inc.xyz; gdal_translate -of XYZ -b 2 %s/los.rdr.geo.clip.vrt %s/gbis_heading.xyz',strcat(x_min),strcat(y_min),strcat(x_max),strcat(y_max), los_file, temp_dir, temp_dir,temp_dir,temp_dir,temp_dir);
system(inc_head);
% elseif tops_flag==1
% inc_head = sprintf('gdalwarp -overwrite  -multi -te %s %s %s %s %s.vrt %s/los.rdr.geo.clip.vrt; gdal_translate -of XYZ -b 1 %s/los.rdr.geo.clip.vrt %s/gbis_inc.xyz; gdal_translate -of XYZ -b 2 %s/los.rdr.geo.clip.vrt %s/gbis_heading.xyz',strcat(x_min),strcat(y_min),strcat(x_max),strcat(y_max), isce_base_dir, isce_dir, isce_dir,isce_dir,isce_dir,isce_dir);
% system(inc_head);    
% else
% inc_head = sprintf('gdalwarp -overwrite  -multi -te %s %s %s %s %s/los.rdr.geo.vrt %s/los.rdr.geo.clip.vrt; gdal_translate -of XYZ -b 1 %s/los.rdr.geo.clip.vrt %s/gbis_inc.xyz; gdal_translate -of XYZ -b 2 %s/los.rdr.geo.clip.vrt %s/gbis_heading.xyz',strcat(x_min),strcat(y_min),strcat(x_max),strcat(y_max), isce_base_dir, isce_dir, isce_dir,isce_dir,isce_dir,isce_dir);
% system(inc_head);
% end 
%% LOAD DATA 
INC_INPUT = importdata([temp_dir '/gbis_inc.xyz']);
PHASE_INPUT = importdata([temp_dir '/gbis_lonlatphase.xyz']);
HEADING_INPUT = importdata([temp_dir '/gbis_heading.xyz']);
if exist('cp_flag')
    CONN_PARTS = importdata([temp_dir '/gbis_unw_conncomp.xyz']);
end

GBIS_NAME = output_name;

%% PREPARE GBIS INPUTFILE
%define satellite heading from pixel heading
HEADING_INPUT=HEADING_INPUT(:,3)*-1+90;
c=find(HEADING_INPUT(:,1)== 90);
HEADING_INPUT(c,:)=0;

%CREATE MATRIX
GBIS=[PHASE_INPUT INC_INPUT(:,3) HEADING_INPUT];

if exist('cp_flag')
%mask phase with connected parts
GBIS(CONN_PARTS(:,3)==0,:)=[];
end

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
if exist('cp_flag')    
    crop_lim=find(CONN_PARTS(:,1)>crop(:,1) & CONN_PARTS(:,1)<crop(:,2) & CONN_PARTS(:,2)>crop(:,3) & CONN_PARTS(:,2)<crop(:,4));
    CONN_PARTS=CONN_PARTS(crop_lim,:);
end
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
if exist('cp_flag')
Conn_part=double(CONN_PARTS(CONN_PARTS(:,3)~=0,3));
end
if ~exist('output_name')
    output_name = 'insar_gbis_input'
end

mat_name=strcat(output_dir,'/',output_name,'.mat');
if exist('cp_flag')
    save(mat_name,'Lon','Lat','Phase','Inc','Heading','Conn_part')
else
    save(mat_name,'Lon','Lat','Phase','Inc','Heading')
end
    fprintf('Saving GBIS INSAR input file: %s \n',mat_name);
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
if exist('cp_flag')
    Conn_part=double(CONN_PARTS(CONN_PARTS(:,3)~=0,3));
end
mat_name=strcat(output_dir,'/',output_name,'_sparse.mat');

if exist('cp_flag')
    save(mat_name,'Lon','Lat','Phase','Inc','Heading','Conn_part')
else
    save(mat_name,'Lon','Lat','Phase','Inc','Heading')
end

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
    if seed_flag==1
    a=plot(seed_ref(1,1),seed_ref(1,2),'ks','MarkerFaceColor','black','MarkerSize',10);
    legend(a,'Reference Point')
    end
    hold off    
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
    if exist('cp_flag')
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
 end  
 rmdir(temp_dir,'s');
 toc 

%% Hydra Mouth Opening Script

%%Modified: March 17th,2022, Tapan Goel: Updated all analysis to use pixels
%%and frames. Only convert to seconds and um at the very end to save data

% Calculates the area of hypostome for each frame in sequence.
% Manually enter the first frame of hypostome opening in stack (num)
% Manually enter the frame rate (interval)
% Manually enter the pixel conversion (conversion)  
% Run program and manually select folder containing images to be analyzed
% Manually click points around hypostome in first frame of opening
% Script will calculate hypostome area in each successive image in folder
% Output is an Matlab cell array in the folder containing the images 

% NOTE: If you are using a Mac you have to replace all of the backslashes
% with forward slashes


clear all; close all;

%% Select folder containing movie
folder = uigetdir('','Select Image Folder');
folder = [folder '\'];

if isempty(dir([folder '*.tiff'])) ~= 1
     file_list = dir([folder '*.tiff']);
else
     file_list = dir([folder '*.tif']);
end

%% Get movie details from user
[first_frame, last_frame, gap, framerate, pixpermm] = MovieDetails();

%% Analyze first image

% Size theshold for registering opening of the hypostome
cell_diam = 15; % microns
cell_area = pi*((cell_diam/2)^2); % area in microns
sz_th = cell_area/4; % minimum area for the hypostome opening to be counted
th_rng = 0.005;


% Cell-size structuring element
cell_pix = round(cell_diam/(pixpermm*2));


%% Intensity threshold for initial opening
h = figure;
im = im2double(imread([folder file_list(first_frame).name])); %read first image
im_filt = imgaussfilt(im,0.5,'FilterSize',3,'Padding','replicate'); %Blur image with a 3x3 gauss filter of sigma = .5
init_bw = roipoly(im_filt); % Mark the region of interest - mouth, in the image;

in_bound = init_bw-imerode(init_bw, strel('disk',2, 0)); %create a 2 pix thick mask of the mouth boundary from inside. the sturcturing element is a disk of radius 2. creates rounded corners in the new mask.
out_bound = imdilate(init_bw, strel('disk',2, 0))-init_bw; % create a 2 px thick mask of the mouth boundary from outside.

mean_intensity_in = sum(in_bound.*im_filt,'all')/sum(in_bound,'all'); %Mean intensity of the in boundary pixels
mean_intensity_out = sum(out_bound.*im_filt,'all')/sum(out_bound,'all'); %Mean intensity of the out boundary pixels

init_threshold = mean_intensity_in; %%Set the binarizing threshold for the first frame as the mean intensity of the inner boundary

init_stats = regionprops(init_bw, 'Centroid', 'Area');
mouth_x = init_stats(1).Centroid(1);
mouth_y = init_stats(1).Centroid(2);

close(h);


%% Initialize storage table


%% Loop over all frames.
image_number = 0;

for i = first_frame:gap:last_frame
    image_number = image_number+1
    
    %% Read and smooth image
    im = im2double(imread(strcat(folder,'\',file_list(i).name)));
    im_filt = imgaussfilt(im,0.5,'FilterSize',3,'Padding','replicate'); %Blur image with a 3x3 gauss filter of sigma = .5
    
    %% Find threshold for binarizing, search on either side of initial threshold
    %%For first frame, use the mean inside intensity.
    %%For all othe frames, scan over range of threshold values (+/- .005 of the initial threshold) and pick the
    %%one that matches the best centroid and average contrast.
    if i == first_frame
        binarythreshold = init_threshold;
    else        
        binarythreshold = ThresholdFinder(im_filt,init_threshold-0.005:.001:init_threshold+0.005,mouth_x,mouth_y);        
    end
    
    %% Now that threshold has been found, process the image to get mouth area and other mouth features:    
    b = im2bw(im_filt,binarythreshold);
    bc_dil = imdilate(imopen(imfill(1.-b,'holes'),strel('disk',2, 0)),strel('disk',2, 0));    
    ConnectedRegions = regionprops(bwconncomp(bc_dil), 'Centroid');
    
    if(length(ConnectedRegions) == 0)
        return;
    else %% if there are multiple regions, identified, pick the one closest to the previous mouth
        centroids = [ConnectedRegions.Centroid];
        centroids_x = centroids(1:2:end);
        centroids_y = centroids(2:2:end);
        centroid_distance = sqrt( (centroids_x-mouth_x).^2 + (centroids_y-mouth_y).^2);         
        [mindist, mouth_centroid_idx] = min(centroid_distance);        
        mouth_centroid = ConnectedRegions(mouth_centroid_idx).Centroid;
        bc_dil = bwselect(bc_dil,mouth_centroid(1),mouth_centroid(2)); %%isolate the region of the mouth
        
      % Remove "pockets" (uncomment this if little parts of cells at the
      % edge of the opening are being picked up.)
      
        if sum(bc_dil,'all') > 6*cell_area            
            bc_dil = imerode(bc_dil,strel('disk', cell_radius, 0));
            bc_dil = imdilate(bc_dil,strel('disk', cell_radius, 0));            
        end
              
        if length(regionprops(bc_dil,'Area')) > 1
             
                for J = 1:length(regionprops(bc_dil,'Area'))
                    
                    abc(J) = area(J).Area;
                    disp('Multiple regions detected')
                    
                end 
                
                bc_dil = bwareaopen(bc_dil, max(abc));
                mouthstats = regionprops('table',bc_dil, 'Area','Centroid','Perimeter',...
                       'Orientation','MinorAxisLength','MajorAxisLength');
                
        elseif length(area) == 1
             
             mouthstats = regionprops('table',bc_dil, 'Area','Centroid','Perimeter',...
                       'Orientation','MinorAxisLength','MajorAxisLength');
             
        else
             
             disp('Ending sequence: region of zero area detected')
             break;
             
        end 
                              
    end

    %% If the area of the mouth is greater than threshold, update mouth centroid
    if sum(bc_dil,'all') > mouthareathreshold
                
        bc_erd = imerode(bc_dil,strel('disk', cell_radius, 0));
        border = bc_dil-bc_erd;
        mouth_x = mouth_centroid(1);
        mouth_y = mouth_centroid(2);
        MouthProperties = [MouthProperties; mouthstats];
        
    else
        
       break;
        
    end
     
    
    displayim(:,:,1) = im;
    displayim(:,:,2) = im + border;
    displayim(:,:,3) = im;
    imshow(displayim,'InitialMagnification','fit')
    title(sprintf('Image %03i of %03i',i,lastframe))
    pause(0.3);
end
      
        
    



% If the script messes up it will still save the data up to when it broke. 
if (i < last_frame)
    ending_frame = i-1;
else
    ending_frame = i;   
end

%% Save workspace
MouthProperties.Frame = first_frame:gap:ending_frame;

cd(folder);
save_folder = sprintf('HM_Data_Images_%i_Through_%i-t2',num,ending_frame);
mkdir(save_folder)
cd(save_folder)
save_file = sprintf('Data_Images_%i_Through_%i-t2.mat',num,ending_frame);
save(save_file,'DataCell');
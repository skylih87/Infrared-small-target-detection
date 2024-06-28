clear;close all;
input_folder    = 'image_dataset/SNU';  % source folder
ext             = {'*.jpg','*.png','*.bmp'};
images          = [];
for i = 1:length(ext)
    images = [images;dir(fullfile(input_folder, ext{i}))];
end

for i = 1:numel(images)
    [~, name, exte] = fileparts(images(i).name);
    I   =   imread( fullfile(input_folder,images(i).name) ) ;
    [height,width,~]=size(I);
    I=imcrop(I,[1 1 width height]);
    [height,width,chennal]=size(I);
    if chennal > 1
        I = rgb2gray(I);
    end

    Zs = stdfilt(I);
    Z = double(I);
    Zs=Zs.*Z;
    x=1:width;
    y=1:height;
    [X,Y]=meshgrid(x,y);

    figure(1)
    mesh(X,-Y,Zs);hold on;
    xlabel('x Pixel');ylabel('y Pixel');zlabel('Intensity')
    hold off;

%% roughly finding
    rate = 0.5;  % gradient theshold
    max_rate = 0.9;  % maximum theshold
    mean_std=mean2(Zs);
    max_std=max(Zs(:));
    
    [col_candi,row_candi]=find(Zs > rate*max_std); 
    points_candi = [row_candi,col_candi];  %candidated pixels

    epsilon=10;
    MinPts=1;
    IDX_candi=DBSCAN(points_candi,epsilon,MinPts); % clustering
    k_candi=max(IDX_candi); % No. clustering
    mean_points_candi=zeros(k_candi,2);
    retina_pixel_FOV = 0.1; % mrad / pixel
    small_rect = (1.2/retina_pixel_FOV/2);
    large_rect = round(3/retina_pixel_FOV/2);

    g=figure;
    subplot(3,2,1)
    imshow(I);hold on;
    PlotClusterinResult(points_candi, IDX_candi)
    subplot(3,2,2)
    imshow(I);hold on;
    Colors=hsv(k_candi);
    Legends = {};
    cumulative_points=[];
    min_width_large=zeros(k_candi,1);
    max_width_large=zeros(k_candi,1);
    min_heigh_large=zeros(k_candi,1);
    max_heigh_large=zeros(k_candi,1);
    for j=1:k_candi
        if size(points_candi(IDX_candi==j,:),1) == 1
            mean_points_candi(j,:)=points_candi(IDX_candi==j,:);
        else
            mean_points_candi(j,:)=mean(points_candi(IDX_candi==j,:));
        end
        min_width_large(j) = round(min(mean_points_candi(j,1)))-large_rect; % row
        max_width_large(j) = round(max(mean_points_candi(j,1)))+large_rect;
        min_heigh_large(j) = round(min(mean_points_candi(j,2)))-large_rect; % column
        max_heigh_large(j) = round(max(mean_points_candi(j,2)))+large_rect;

        if min_width_large(j) < 1
            min_width_large(j) = 1;
        end
        if max_width_large(j) > width
            max_width_large(j) = width;
        end
        if min_heigh_large(j) < 1
            min_heigh_large(j) = 1;
        end
        if max_heigh_large(j) > height
            max_heigh_large(j) = height; 
        end
        mask_candi = zeros(height,width);
        mask_candi(min_heigh_large(j):max_heigh_large(j),min_width_large(j):max_width_large(j)) = 1;
        meanz=mean2(mask_candi.*Z);
        maxz=max(max(mask_candi.*Z));
        [col{j},row{j}]=find(mask_candi.*Z > meanz + max_rate*(maxz-meanz));

        Color = Colors(j,:);
        Legends{end+1} = ['Candidate #' num2str(j)];
        plot(mean_points_candi(j,1),mean_points_candi(j,2),'+','Color',Color)
        rectangle('Position',[min_width_large(j) min_heigh_large(j) 2*large_rect 2*large_rect],'EdgeColor',Color,'LineWidth',1)

        points{j}=[row{j} col{j}];
        cumulative_points=[cumulative_points;points{j}];
    end
hold off;
legend(Legends);

%
subplot(3,2,3)
imshow(I);hold on;

count=0;
cum_IDX=DBSCAN(cumulative_points,epsilon/2,MinPts*2);
k_obj=max(cum_IDX);
new_points=[];
PlotClusterinResult(cumulative_points,cum_IDX)
if k_obj == 0
    new_points = cumulative_points;
    k_obj = size(new_points,1);
end
Colors=hsv(k_obj);
Legends = {};


subplot(3,2,4)
imshow(I);hold on;

for j=1:k_obj
    if cum_IDX(j) ~= 0 || max(cum_IDX) ~= 0
        new_points(j,:)=mean(cumulative_points(cum_IDX==j,:));
        new_points_min(j,1)=min(cumulative_points(cum_IDX==j,1));
        new_points_min(j,2)=min(cumulative_points(cum_IDX==j,2));
        new_points_max(j,1)=max(cumulative_points(cum_IDX==j,1));
        new_points_max(j,2)=max(cumulative_points(cum_IDX==j,2));
    else
        new_points_min(j,1)=min(new_points(:,1));
        new_points_min(j,2)=min(new_points(:,2));
        new_points_max(j,1)=max(new_points(:,1));
        new_points_max(j,2)=max(new_points(:,2)); 
    end   

    Color = Colors(j,:);
    Legends{end+1} = ['Object #' num2str(j)];
    plot(new_points(j,1),new_points(j,2),'.','Color',Color)    
    rectangle('position',[new_points_min(j,:) new_points_max(j,:)-new_points_min(j,:)],'EdgeColor',Color,'LineWidth',2)

end
legend(Legends);
Legends = {};

n=0;
target_points=zeros(k_obj,2);
for j=1:k_obj
    if new_points_max(j,1)-new_points_min(j,1)<=small_rect*3 && new_points_max(j,2)-new_points_min(j,2)<=small_rect*3
        n=n+1;
        target_points(n,:)=new_points(j,:);
    end
    if j==k_obj && n==0
        n=k_obj;
        target_points(n,:)=new_points(n,:);
    end
end

min_width=zeros(n,1);
max_width=zeros(n,1);
min_heigh=zeros(n,1);
max_heigh=zeros(n,1);

patch = zeros(round(2*small_rect), round(2*small_rect));
number=0;
len = 4;
nr = 3;
nc = 3;
leny = len*nr;
lenx=  len*nc;
op = zeros(leny, lenx, nr*nc);
for ii = 1:nr*nc
    temp = zeros(len*nr,  len*nc);
    [r, c]  = ind2sub([nr, nc], ii);
    temp((r-1)*len + 1:r*len, (c-1)*len+1:c*len) = 1;
    temp = temp';
    op(:, :, ii) = temp;
end
subplot(3,2,5)
imshow(I);hold on;
intensity = zeros(n,nr*nc);

for k=1:n
    min_width(k) = round(target_points(k,1))-small_rect; % row [x]
    max_width(k) = round(target_points(k,1))+small_rect;
    min_heigh(k) = round(target_points(k,2))-small_rect; % column [y]
    max_heigh(k) = round(target_points(k,2))+small_rect;

    if min_width(k) < 1
        min_width(k) = 1;
    end
    if min_width(k) > width-small_rect*2+1
        min_width(k) = width-small_rect*2+1;
    end
    if min_heigh(k) < 1
        min_heigh(k) = 1;
    end
    if min_heigh(k) > height-small_rect*2+1
        min_heigh(k) = height-small_rect*2+1; 
    end
    patch = imcrop(I,[min_width(k) min_heigh(k) 2*small_rect-1 2*small_rect-1]);
    
    for ii = 1:nr*nc
        gimg(:, :, ii) = double(patch).*op(:,:,ii);
        intensity(k,ii)=max(max(gimg(:,:,ii)));
    end

    if intensity(k,5) - max(intensity(k,[1:4,6:end])) >= 5
        number=number+1;
        target_position{i}(number,:)=round(target_points(k,:));

    end
    if number == 0 && n==1
        number=1;
        target_position{i}(number,:)=round(target_points(n,:));
    elseif number == 0
        target_position{i}=[];
    end

    Color = Colors(k,:);
    Legends{end+1} = ['Target #' num2str(k)];
    plot(target_points(k,1),target_points(k,2),'+','Color',Color)    
    rectangle('Position',[min_width(k) min_heigh(k) 2*small_rect 2*small_rect],'EdgeColor',Color,'LineWidth',1)
end

legend(Legends);
hold off;

subplot(3,2,6)
imshow(I);hold on;
for check=1:size(target_position{i},1)
    Color = Colors(check,:);
    rectangle('Position',[target_position{i}(check,1)-small_rect target_position{i}(check,2)-small_rect 2*small_rect 2*small_rect],'EdgeColor',Color,'LineWidth',1)
end
pause
end
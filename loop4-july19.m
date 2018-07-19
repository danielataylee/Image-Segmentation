
img_mask   = imread('testmask_flip_horiz.tiff','tiff');
img   = imread('test_flip_hor.tiff','tiff');
s=1;
o=[];
j=0;
BW = im2bw(img_mask,0.8); % convert image to a binary image
preproc={};
preproc{1,1} = BW;
preprocim{1,1} = img;

while 1
    s=s+1;
    
    if s==100
        break
    end
    j = j+1;
    if j > length(preproc)  % check if preprocessed list of daughter object have all been gone through
        break               % if they have then break out of while loop
    end
    [B,L,N] = bwboundaries(preproc{j,1},'holes');
    if s== 40
        w=w;
    end
    LENGTH=[];
    %choose largest object when multiple objects detected
    if N>1
        for i =1:N
            sz = size(L);
            %check if the object has boundary element on the edge of image. If it
            %does it is not the object wanted
            if isempty(find(B{i} == 0 )) == 0 || isempty(find(B{i}(:,1) == sz(1))) == 0 || isempty(find(B{i}(:,2) == sz(2))) == 0
            else
                LENGTH = [LENGTH;[i,length(B{i})]];
            end
        end
         ind = find(LENGTH(:,2) == max(LENGTH(:,2)));
         MAX_L_ind = LENGTH(ind,1);
    else
        MAX_L_ind=1;
    end
    b=B{MAX_L_ind};
    
    %pause here for checks
    if s== 84
        w=w;
    end
    
    figure
    imshow(preproc{j,1})
    hold on
    plot(b(:,2),b(:,1),'g');
    %hold off
    
    %     figure
    %     imshow(preprocim{j,1})
    %     hold on
    %%%%%%%%%%%===================================> check cell size!!!!!
    %count the number of elements with the label of the object we want
    %which is equivalent to its size
    num = sum(L(:) == MAX_L_ind);
    if num < 1500
        continue
    end
    
    [ANGLES,X1,Y1,INFS,Po]=GetAng(b,5,preproc{j,1}); %FUNCTION TO GET THE ANGLES WHERE INPUT 1 IS BOUNDARY ARRAY AND 2ND IS THE REACH
    if length(Po) <=2
        continue;
    end
    
    for i = 1:length(Po)
        plot(b(Po(:,1),2),b(Po(:,1),1),'r.');
    end
    figure
    imshow(preprocim{j,1})
    hold on
    [P_D_val_x,P_D_val_y] = GetDnP(b,Po);
    if isempty(P_D_val_x) == 1
        continue
    end
    
    [LabB,nL] = bwlabel(preproc{j,1}) ;%label objects
    
    
    while 1
        for i =1:length(P_D_val_x)
            [BW2,STAT] = GetBI(P_D_val_x(i),P_D_val_y(i),preproc{j,1},b); %cut new objects
            if STAT==0
                continue %if line drawn crosses outside of object then go back up
            end
            [B2,L2,N2] = bwboundaries(BW2,4,'holes');
            %make sure the objects selected are from the same original
            %object
            if nL>=2
                W=[];
                for k = 1:N2
                    a= B2{k};
                    [xm,I]= max(a(:,1));
                    ym =a(I,2);
                    w=LabB(xm,ym);
                    W=[W;w];
                end
                inds=[];
                for l=1:length(W)-1
                    for m=l+1 : length(W)
                        if W(l)==W(m)
                            inds = [l,m];
                            break
                        end
                    end
                    if isempty(inds) == 0
                        break
                    end
                end
            else %if only 2 object or less are found inds array is automatically 1 and 2
                inds = [1,2];
            end
            
            %%%%%%%%%%%%%%%%%%% Test for eccentricity of the 2 new objects
            b21= B2{inds(1)}; % boundary of object 1
            num21 = sum(L2(:) == inds(1)); %size of object 1
            b22 = B2{inds(2)}; % boundary of object 2
            num22 = sum(L2(:) == inds(2)); % size of object 2
            
            
            if num < 20000 && num > 5000  %check size to do circularity test since large objects are predicted to segment into non circular objects
                Ecc21 = (4*pi*num21)/(length(b21)^2); %eccentricity of obj 1
                Ecc22 = (4*pi*num22)/(length(b22)^2); %eccentricity of obj 2
  
                if Ecc21 > 0.2 || Ecc22 > 0.2
                else
                    if i ==length(P_D_val_x)
                        break
                    else
                        continue       %if it is too uncircular go back up forloop
                    end
                end
            end
                
            % check size of daughter objects
            if num21 <200 || num22 <200 %%%==================> check size of daughter objects
                if i ==length(P_D_val_x)
                    break
                else
                    continue
                end
            end
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%
            %split this is first half
            
            [bleh21]=GetcropImg(b21,BW2); %crop object mask
            [bleh21img]=GetcropImg(b21,preprocim{j,1}); %%%% crop object image
            figure
            hold on
            imshow(bleh21);
            figure
            imshow(bleh21img); %%%%
            M = 1+length(preproc);
            preproc{M,1} = bleh21;
            preprocim{M,1} = bleh21img; %%%%
            
            %2nd half
            [bleh22]=GetcropImg(b22,BW2);
            [bleh22img]=GetcropImg(b22,preprocim{j,1}); %%%%
            figure
            hold on
            imshow(bleh22);
            figure
            imshow(bleh22img); %%%%
            M = 1+length(preproc);
            preproc{M,1} = bleh22;
            preprocim{M,1} = bleh22img; %%%%
            
            
            break
            
        end
        
        break
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%FUNCTIONS
function [Bin_Im,stat] = GetBI(P_x,P_y,BLW,cc)

COORD_X1 = cc(P_x,2);
COORD_X2 = cc(P_y,2);
COORD_Y1 = cc(P_x,1);
COORD_Y2 = cc(P_y,1);
Bin_Im = BLW;
stat =1;

burnedImage = BLW;
% Create line mask, h, as an ROI object over the second image in the bottom row.
figure('visible','off');
imshow(burnedImage);
hLine = imline(gca,[COORD_X1 COORD_X2],[COORD_Y1 COORD_Y2]); % Second argument defines line endpoints.
hLine2 = imline(gca,[(COORD_X1+1) (COORD_X2+1)],[COORD_Y1 COORD_Y2]);
% Create a binary image ("mask") from the ROI object.
binaryImage2 = hLine.createMask();
binaryImage21 = hLine2.createMask();
%combine 2 lines to make it 2 pixels thick
binaryImagecomb = binaryImage2 | binaryImage21;
% Burn line into image by setting it to 0 (black) wherever the mask is true.

%check that the line does not pass outside the object
BBI = burnedImage(binaryImagecomb);
%Find elements in image on line cordinates that have value 0
%if find less than 2 pixels that are 0 (which may be the inflection points
%then line does not pass through background. If it is larger then it does
%pass background
if length(find(BBI==0)) > 2
    stat = 0;
end

burnedImage(binaryImagecomb) = 0;
% Display the image with the "burned in" line.
%subplot(2, 4, 8);
%imshow(burnedImage);

Bin_Im = burnedImage;
close(gcf);
%imshow(Bin_Im);

end
function [ANGS,XX,YY,INF,P] = GetAng(c,reach,BLW)
ANGS =[];
XX=[];
YY=[];
INF=[];
BOO=[];
%infl_p_x =[];
%infl_p_y =[];
for i = 1:length(c(:,1))
    X = c(i,2);
    XX = [XX;X];
    Y = c(i,1);
    YY = [YY;Y];
    
    if i>reach && i<=length(c(:,1))-reach %section of the boundary chain code
        
        %Calculate angle between left and right reach position with middle point
        minus_i = i-reach;
        plus_i = i+reach;
        P_left = [c(minus_i,2);c(minus_i,1)];
        P_right = [c(plus_i,2);c(plus_i,1)];
        P_mid = [c(i,2);c(i,1)];
        
        
        
        Pml = sqrt((P_mid(1)-P_left(1))^2 + (P_mid(2)-P_left(2))^2);
        Pmr = sqrt((P_mid(1)-P_right(1))^2 + (P_mid(2)-P_right(2))^2);
        Plr = sqrt((P_left(1)-P_right(1))^2 + (P_left(2)-P_right(2))^2);
        
        ANGLE = acos((Pml^2 + Pmr^2 - Plr^2)/(2*Pml*Pmr));
        ANGLE_DEG = rad2deg(ANGLE);
        
        ANGS = [ANGS;ANGLE_DEG];
        
        %%%%%%%%%%%%%%checks if angle is from a bulge or inflection point
        m_x = round(mean([P_right(1),P_left(1)]));
        m_y = round(mean([P_right(2),P_left(2)]));
        %sm = BLW(m_x,m_y);
        %co_x = c(m_x,1);
        %co_y = c(m_y,2);
        BOO= [BOO; [BLW(m_y,m_x),i,m_x,m_y]];
        
    end
end

for j =1: length(ANGS)
    
    if ANGS(j) <150 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%====> change inflection angle here!!!!!
        minus_j = j-reach;
        plus_j = j+reach;
        %record the index and the angle of ones that meet the requirement
        %angle
        if BOO(j,1) ==0 %if logic value is zero then it is an inflection point
            INF = [INF;[plus_j,ANGS(j)]];
        end
        %INF = [INF;[plus_j,ANGS(j)]];
        
    end
end
l=1;
L=1;
H=1;
t =1;
INF_filt=[];
P=[];

while l==1
    if isempty(INF)==1 || length(INF(:,1)) <2
        break
        
    else
        if t>length(INF(:,1))
            [val,ind]= sort(INF_filt);
            if length(INF_filt)==2 || length(INF_filt)==1
                P = [P;[INF_filt(1),INF_filt(2)]];
            else
                s = val(:,1);
                s=s(ind(:,2));
                P = [P;[s(1),val(1,2)]];
            end
            l=0;
            
            
        elseif INF(L,1) < INF(H,1) + 30  %%%%%%%%%%%%%%%%%=>print smallest angle within ... indices
            INF_filt = [INF_filt; [INF(L,1),INF(L,2)]];
            L=L+1;
            t=t+1;
        else
            [val,ind]= sort(INF_filt);
            if length(INF_filt)==2 || length(INF_filt)==1
                P = [P;[INF_filt(1),INF_filt(2)]];
            else
                s = val(:,1);
                s=s(ind(:,2));
                
                P = [P;[s(1),val(1,2)]];
            end
            H=t;
            INF_filt =[];
            
            
        end
    end
end

end
% function to get distance & perimeter
function [one,two] = GetDnP(c,p)
DIST =[];
PERIM=[];
coord_start=[];
coord_end=[];
P_D_val = [];
%P_P_val = [];
one =[];
two=[];
for i = 1: length(p)-1
    for j = i+1: length(p)
        %distance between inflection points using x and y coordinates
        dist = sqrt((c(p(i,1),1)-c(p(j,1),1))^2 + (c(p(i,1),2)-c(p(j,1),2))^2);
        DIST = [DIST;dist];
        
        %perimeter calculated by finding the difference between indices of
        %inflection points
        perim = sqrt((p(j,1)-p(i,1))^2);%cw perim
        perimccw = p(i,1) + (length(c)-p(j,1));
        if perimccw < perim
            PERIM = [PERIM;perimccw];
        else
            PERIM = [PERIM;perim];
        end
        coord_start = [coord_start; p(i,1)];
        coord_end = [coord_end; p(j,1)];
    end
end
[D_val,D_ind] = sort(DIST);
P_sort = PERIM(D_ind);
coord_start_sort = coord_start(D_ind);
coord_end_sort = coord_end(D_ind);

for t = 1:length(P_sort)
    if P_sort(t)>100  %%%%%%%%%%%%%%%%%%=================> change perimeter length here!!!!
        P_D_val = [P_D_val; D_val(t)]; %list of ascending distance
        %P_P_val = [P_P_val; P_sort(t)];
        one = [one;coord_start_sort(t)]; %list of coordinates sorted by ascending distance
        two = [two;coord_end_sort(t)];
    end
end
end
function [IMAGE] = GetcropImg(c,i) %input is boundary array and bachground image
max_x = max(c(:,2))+3;
max_y = max(c(:,1))+3;
min_x = min(c(:,2))-3;
min_y = min(c(:,1))-3;
rect = [min_x-1,min_y,(max_x-min_x),(max_y-min_y)];
IMAGE = imcrop(i,rect);
end
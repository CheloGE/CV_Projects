function M = ImageMorphingTriangulation(warp_frac,dissolve_frac)

if nargin < 1
    warp_frac = .5;
end

if nargin < 2
    dissolve_frac= warp_frac; 
end


% ream images
I = im2double(imread('a.png'));
J = im2double(imread('c.png'));

% load mat file with points, variables Ip,Jp
 load('points.mat');
 
% initialize ouput image (morphed)
M = zeros(size(I));

%  Triangulation (on the mean shape)
MeanShape = (1/2)*Ip+(1/2)*Jp;
TRI = delaunay(MeanShape(:,1),MeanShape(:,2));


% number of triangles
TriangleNum = size(TRI,1); 

% find coordinates in images I and J
CordInI = zeros(3,3,TriangleNum);
CordInJ = zeros(3,3,TriangleNum);

for i =1:TriangleNum
  for j=1:3
    
    CordInI(:,j,i) = [ Ip(TRI(i,j),:)'; 1];
    CordInJ(:,j,i) = [ Jp(TRI(i,j),:)'; 1]; 
    
  end
end

% create new intermediate shape according to warp_frac
Mp = (1-warp_frac)*Ip+warp_frac*Jp; 

 
% create a grid for the morphed image
[x,y] = meshgrid(1:size(M,2),1:size(M,1));

% for each element of the grid of the morphed image, find  which triangle it falls in
TM = tsearchn([Mp(:,1) Mp(:,2)],TRI,[x(:) y(:)]);



% YOUR CODE STARTS HERE
indNan=find(isnan(TM));
TM(isnan(TM)) = 1;

for i = 1:size(x,1)
    for j = 1:size(x,2)
        b=[x(i,j);y(i,j);1];
        k=sub2ind(size(x),i,j);
        

        tri_vert=TRI(TM(k),:);
        Ma = [Mp(tri_vert(1),1) Mp(tri_vert(2),1) Mp(tri_vert(3),1); 
            Mp(tri_vert(1),2) Mp(tri_vert(2),2) Mp(tri_vert(3),2);
            1 1 1];
        
        sol = Ma\b;
        
        alpha=sol(1);
        beta=sol(2);
        gamma=sol(3);
        
        xyz_I=CordInI(:,:,TM(k))*[alpha;beta;gamma];
        xy_I=round([(xyz_I(2)/xyz_I(3)),(xyz_I(1)/xyz_I(3))]);
        
        xyz_J=CordInJ(:,:,TM(k))*[alpha;beta;gamma];
        xy_J=round([(xyz_J(2)/xyz_J(3)),(xyz_J(1)/xyz_J(3))]);
       
        IndM(i,j, 1) = sub2ind(size(I),i,j, 1);

        IndM(i,j, 2) = sub2ind(size(I),i,j, 2);

        IndM(i,j, 3) = sub2ind(size(I),i,j, 3);


        for k_nan= 1:length(indNan)  
            if indNan(k_nan)==k
                xy_I(1)=1;
                xy_I(2)=1;
                xy_J(1)=1;
                xy_J(2)=1;
                break
            end
        end

        IndI(i,j, 1) = sub2ind(size(I),xy_I(1), xy_I(2) , 1);

        IndI(i,j, 2) = sub2ind(size(I),xy_I(1), xy_I(2) , 2);

        IndI(i,j, 3) = sub2ind(size(I),xy_I(1), xy_I(2) , 3);
        
        IndJ(i,j, 1) = sub2ind(size(J),xy_J(1), xy_J(2) , 1);

        IndJ(i,j, 2) = sub2ind(size(J),xy_J(1), xy_J(2) , 2);

        IndJ(i,j, 3) = sub2ind(size(J),xy_J(1), xy_J(2) , 3);
        
    end
    
end 
% YOUR CODE ENDS HERE



% cross-dissolve
M(IndM)=(1-dissolve_frac)* I(IndI)+ dissolve_frac * J(IndJ);


figure(100);
subplot(1,3,1);
imshow(I);
hold on;
triplot(TRI,Ip(:,1),Ip(:,2))
hold off;
title('First')

subplot(1,3,2);
imshow(M);
hold on;
triplot(TRI,Mp(:,1),Mp(:,2))
hold off
title('Morphed')

subplot(1,3,3);
imshow(J);
hold on;
triplot(TRI,Jp(:,1),Jp(:,2))
hold off
title('Second')

[length(find(M~=I)),length(find(M~=J))]
[length(find(Mp~=Ip)), length(find(Mp~=Jp))]
size(M)
end


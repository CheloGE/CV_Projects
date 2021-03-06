function [seamless,A] = PoissonImageEditing

% read images 
target= im2double(imread('target_3.jpg')); 
source= im2double(imread('source_3.jpg')); 

% mask image
mask=imread('mask_3.bmp');

% image offsets
row_offset=20;
col_offset=40; 


% N: Number of pixels in the mask
N=sum(mask(:)); 

% enumerating pixels in the mask
mask_id = zeros(size(mask));
mask_id(mask) = 1:N;   
    
% neighborhood size for each pixel in the mask
[ir,ic] = find(mask);

Np = zeros(N,1); 

for ib=1:N
    
    i = ir(ib);
    j = ic(ib);
    
    Np(ib)=  double((row_offset+i> 1))+ ...
             double((col_offset+j> 1))+ ...
             double((row_offset+i< size(target,1))) + ...
             double((col_offset+j< size(target,2)));
end




% compute matrix A

% We use a sparse matrix since most of the data is zero
% and we only want to track non zero values (optimizing resources)

A = sparse(1:N,1:N,Np,N,N);
for curr_pixel=1:N
    % get row and column of current pixel in the mask
    row = ir(curr_pixel);
    col = ic(curr_pixel);
    A_row = A(curr_pixel,:);
    
    if(mask(row-1,col)~=0)
        A_row_ind = mask_id(sub2ind(size(mask),row-1,col));
        A_row(A_row_ind)=-1;
    end
    
    if(mask(row+1,col)~=0)
        A_row_ind = mask_id(sub2ind(size(mask),row+1,col));
        A_row(A_row_ind)=-1;
    end
    
    if(mask(row,col-1)~=0)
        A_row_ind = mask_id(sub2ind(size(mask),row,col-1));
        A_row(A_row_ind)=-1;
    end
    
    if(mask(row,col+1)~=0)
        A_row_ind = mask_id(sub2ind(size(mask),row,col+1));
        A_row(A_row_ind)=-1;
    end
    
    A(curr_pixel,:)=A_row;
    
end




% output intialization
seamless = target; 


for color=1:3 % solve for each colorchannel

    % compute b for each color
    b=zeros(N,1);
    
    for ib=1:N
    
    i = ir(ib);
    j = ic(ib);
    
            
      if (i>1) 
          b(ib)=b(ib)+ target(row_offset+i-1,col_offset+j,color)*(1-mask(i-1,j))+...
                          source(i,j,color)-source(i-1,j,color);
      end

      if (i<size(mask,1))
          b(ib)=b(ib)+  target(row_offset+i+1,col_offset+j,color)*(1-mask(i+1,j))+ ...
                           source(i,j,color)-source(i+1,j,color);
      end

      if (j>1)
          b(ib)= b(ib) +  target(row_offset+i,col_offset+j-1,color)*(1-mask(i,j-1))+...
                           source(i,j,color)-source(i,j-1,color);
      end


      if (j<size(mask,2))
          b(ib)= b(ib)+ target(row_offset+i,col_offset+j+1,color)*(1-mask(i,j+1))+...
                         source(i,j,color)-source(i,j+1,color); 
      end     


 

    end

    
     % solve linear system A*x = b;
    % your CODE begins here

    x = A\b;
    % your CODE ends here

   
    


    
    % impaint target image
    
     for ib=1:N
           seamless(row_offset+ir(ib),col_offset+ic(ib),color) = x(ib);
     end
    
end


%figure(1), imshow(target);
%figure(3), imshow(seamless);
%figure(2), imshow(source);

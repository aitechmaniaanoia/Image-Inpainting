function imgout = poisson_blend(im_s, mask_s, im_t)
% -----Input
% im_s     source image (object)
% mask_s   mask for source image (1 meaning inside the selected region)
% im_t     target image (background)
% -----Output
% imgout the blended image

[imh, imw, nb] = size(im_s);

function [A, b] = createAbMatrix(imh,imw, mask_s, im_s, im_t)
    %TODO: initialize counter, A (sparse matrix) and b.
    %Note: A don't have to be k��k,
    %      you can add useless variables for convenience,
    %      e.g., a total of imh*imw variables
    %TODO: fill the elements in A and b, for each pixel in the image

    %fill the element in A
    i = double.empty;
    j = double.empty;
    v = double.empty;
    k = 1;

    for x = 1:imh
        for y = 1:imw
            if mask_s(x, y)==1
                i = [i k];
                j = [j k];
                v = [v 4];
                if mask_s(x-1, y) == 1 %left
                    i = [i k];
                    j = [j (k-1)];
                    v = [v -1];
                end
                if mask_s(x+1, y) == 1 %right
                    i = [i k];
                    j = [j (k+1)];
                    v = [v -1];
                end
                if mask_s(x, y-1) == 1 %up
                    distance = double(sum(double(mask_s(x-1,y:end))) + sum(double(mask_s(x,1:(y-1)))));
                    i = [i k];
                    j = [j (k-distance)];
                    v = [v -1];
                end
                if mask_s(x, y+1) == 1 %down
                    distance = double(sum(double(mask_s(x,(y+1):end))) + sum(double(mask_s(x+1,1:y))));
                    i = [i k];
                    j = [j (k+distance)];
                    v = [v -1];
                end
            k = k+1; 
            end
        end
    end

    A = sparse(i,j,v);
    %spy(A);

    %b
    b = double(zeros(imh, imw));
    % fill the element in b
    for x_m = 1:imh
        for y_m = 1:imw
            if mask_s(x_m, y_m)==1
                b(x_m,y_m) = 4*im_s(x_m,y_m) - im_s(x_m+1,y_m) - im_s(x_m-1,y_m) - im_s(x_m,y_m+1) - im_s(x_m,y_m-1);

                if mask_s(x_m-1,y_m)==0 %left
                    b(x_m,y_m) = b(x_m,y_m) + im_t(x_m-1,y_m);
                end
                if mask_s(x_m+1,y_m)==0 %right
                    b(x_m,y_m) = b(x_m,y_m) + im_t(x_m+1,y_m);
                end
                if mask_s(x_m,y_m-1)==0 %up
                    b(x_m,y_m) = b(x_m,y_m) + im_t(x_m,y_m-1);
                end
                if mask_s(x_m,y_m+1)==0 %down
                    b(x_m,y_m) = b(x_m,y_m) + im_t(x_m,y_m+1);
                end

            elseif mask_s(x_m, y_m)==0
                b(x_m,y_m) = NaN;
            end
        end
    end
    
    b1 = double.empty;
    % convert b from 2d to 1d
    for x_m = 1:imh
        for y_m = 1:imw
            b1 = [b1 b(x_m,y_m)];
        end
    end
    
    b = rmmissing(b1)';

end

%TODO: consider different channel numbers

for c = 1: nb
    ims = im_s(:,:,c);
    imt = im_t(:,:,c);
    [A,b] = createAbMatrix(imh, imw,mask_s, ims, imt);
    %TODO: solve the equation
    %use "lscov" or "\", please google the matlab documents
    solution = A\b;
    error = sum(abs(A*solution-b));
    disp(error)
    result(:,c) = double(solution); 
    
end

%TODO: copy those variable pixels to the appropriate positions
%      in the output image to obtain the blended image
n = 1;
for x_r = 1:imh
    for y_r = 1:imw
        if mask_s(x_r, y_r)==1
            for m = 1:nb
                im_t(x_r, y_r,m) = result(n,m);
            end
            n = n+1;
        end
    end
end
imgout = im_t;
end
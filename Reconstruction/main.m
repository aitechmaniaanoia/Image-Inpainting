clc;
clear;
close all;

imgin = im2double(imread('./target.jpg'));

[imh, imw, nb] = size(imgin);
assert(nb==1);
% the image is grayscale

V = zeros(imh, imw);
V(1:imh*imw) = 1:imh*imw;
%V(y,x) = (y-1)*imw + x
% use V(y,x) to represent the variable index of pixel (x,y)
% Always keep in mind that in matlab indexing starts with 1, not 0

%TODO: initialize counter, A (sparse matrix) and b.
% sparse matrix A
%A = zeros((imh*imw), (imh*imw));
% b
b = double(zeros(imh, imw));

%TODO: fill the elements in A and b, for each pixel50 in the image
% A 
%except edge and corner
%A_test = sparse(1,1,0,2500,2500)

i = double.empty;
j = double.empty;
v = double.empty;

for m = 1:(imh - 2)
    for n = 2:(imw - 1)
       % A_test = A_test + sparse(m*50+n, m*50+n, 4, 2500, 2500)
        %A_test = A_test + sparse(m*50+n, m*50+n+1, -1, 2500, 2500)
        %A_test = A_test + sparse(m*50+n, m*50+n-1, -1, 2500, 2500)
        %A_test = A_test + sparse(m*50+n, m*50+n-50, -1, 2500, 2500)
        %A_test = A_test + sparse(m*50+n, m*50+n+50, -1, 2500, 2500)
        i = [i m*imw+n];
        j = [j m*imw+n];
        v = [v 4];
        
        i = [i m*imw+n];
        j = [j m*imw+n+1];
        v = [v -1];
      
        i = [i m*imw+n];
        j = [j m*imw+n-1];
        v = [v -1];
        
        i = [i m*imw+n];
        j = [j m*imw+n-imw];
        v = [v -1];
        
        i = [i m*imw+n];
        j = [j m*imw+n+imw];
        v = [v -1];
        
    end
end


%left and right edge
for m = 1:(imh - 2)
    for n = [1, imw]
        %A_test = A_test + sparse(m*50+n, m*50+n, 2, 2500, 2500)
        %A_test = A_test + sparse(m*50+n, m*50+n-50, -1, 2500, 2500)
        %A_test = A_50test + sparse(m*50+n, m*50+n+50, -1, 2500, 2500)
        i = [i m*imw+n];
        j = [j m*imw+n];
        v = [v 2];
        
        i = [i m*imw+n];
        j = [j m*imw+n+imw];
        v = [v -1];
        
        i = [i m*imw+n];
        j = [j m*imw+n-imw];
        v = [v -1];
    end
end

%up and down edge
for m = [0, (imh - 1)]
    for n = 2:(imw - 1)
        %A_test = A_test + sparse(m*50+n, m*50+n, 2, 2500, 250500)
        %A_test = A_test + sparse(m*50+n, m*50+n-1, -1, 2500, 2500)
        %A_test = A_test + sparse(m*50+n, m*50+n+1, -1, 2500, 2500)
        i = [i m*imw+n];
        j = [j m*imw+n];
        v = [v 2];
        
        i = [i m*imw+n];
        j = [j m*imw+n+1];
        v = [v -1];        
        
        i = [i m*imw+n];
        j = [j m*imw+n-1];
        v = [v -1];
    end
end

%corner
for m = [0, (imh - 1)]
    for n = [1, imw]
        %A_test = A_test + sparse(m*50+n, m*50+n, 0, 2500, 2500)
        i = [i m*imw+n];
        j = [j m*imw+n];
        v = [v 1];
    end
end

A = sparse(i,j,v);

% b
%except edge and corner
for x = 2:(imh-1)
    for y = 2:(imw-1)
        b(x,y) = 4*imgin(x,y) - imgin(x+1,y) - imgin(x-1,y) - imgin(x,y+1) - imgin(x,y-1);
    end
end
% edgeA = sparse(i,j,v);
spy(A);
for y = 2:imw - 1
    b(1, y)= 2*imgin(1,y) - imgin(1,y-1) - imgin(1,y+1)
    b(imh, y)= 2*imgin(imh,y) - imgin(imh,y-1) - imgin(imh,y+1)
end

for x = 2:imh - 1
    b(x, 1) = 2*imgin(x,1) - imgin(x+1,1) - imgin(x+1,1)
    b(x, imw) = 2*imgin(x,1) - imgin(x+1,1) - imgin(x+1,1)
end
% corner
b(1,1) = 0;
b(1,imw) = 0;
b(imh,1) = 0;
b(imh, imw) = 0;
    
%TODO: add extra constraints
%
b(1,1) = imgin(1,1) + 0.5;
b(1,imw) = imgin(1,imw) + 0.5;
b(imh,1) = imgin(imh,1) - 0.1;
b(imh, imw) = imgin(imh,imw) - 0.1;

%b convert 2d to 1d
b = reshape(b,[],1);

%spy(A)
%TODO: solve the equation
%use "lscov" or "\", please google the matlab documents
solution = A\b;
error = sum(abs(A*solution-b));
disp(error)
imgout = reshape(solution,[imh,imw]);
imwrite(imgout,'output_.png');
figure(), hold off, imshow(imgout);


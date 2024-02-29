% Central Catadioptric Camera Calibration Using Line Image and Mirror Contour
% Inputs:
%    Ci - Line image (at least 2)
%    C0 - Mirror Contour
%
% Outputs:
%    K - 3-by-3 - Intrinsic parameters
%    l - 1-by-1 - Mirror parameter

function [est_T_K, est_T_l] = T(line_image,mirror_contour)
    %%
    %Method T needs at least 3 line images to calibration
    [~,~,Image_num] = size(line_image);
    
    if Image_num < 2
        error('Method T needs at least 2 line images to calibration.')
    end
    
    %%
    %Step 1, find vertices%
    
    l0 = zeros(Image_num,3);
    
    for i = 1:Image_num
         
        
        [V,e] = eig(line_image(:,:,i),mirror_contour);
        e = diag(e);
        if isreal(e) 
            [~,k1] = max(abs(median(sign(e))-sign(e))); 
        else 
            for j = 1:3
                if isreal(e(j))
                    k1 = j;
                end
            end
        end
        v = V(:,k1); %get vertex v
        l0(i,:) = line_image(:,:,i)*v; %get the line passing through the principal point
     end


    %%
    %Step 2, obtain the simplified ICPs%

    
    [~,~,V] = svd(l0); 
    Op = V(:,end);
    Op = Op/Op(3);
    

    Tp=[1 0 0;
        0 1 0;
        Op']';

    l = zeros(Image_num,3);
    line_image_p = zeros(3,3,Image_num);

    for i = 1:Image_num
        line_image_p(:,:,i) = Tp'*line_image(:,:,i)*Tp;
        l(i,:) = line_image_p(:,:,i)*[0 0 1]'; %get vanishing lines 
    end
    
    
    %
    %Function, find the intersections between line L[A1;B1;D1] and conic C=[a1 b1 c1 d1 e1 1]'%

    function [X,Y] = crossLandC(L,C)
    A1 = L(1);
    B1 = L(2);
    D1 = L(3);

    C = C/C(3,3);
    a1 = C(1,1);
    b1 = 2*C(1,2);
    c1 = C(2,2);
    d1 = 2*C(1,3);
    e1 = 2*C(2,3);

    X = zeros(2,1);
    Y = zeros(2,1);
    %sovle
    X(1) = -(D1 + (B1*(A1*(A1^2*conj(e1)^2 + B1^2*conj(d1)^2 + D1^2*conj(b1)^2 - 4*B1^2*conj(a1) - 4*A1^2*conj(c1) - 4*D1^2*conj(a1)*conj(c1) + 4*A1*B1*conj(b1) - 2*A1*B1*conj(d1)*conj(e1) - 2*A1*D1*conj(b1)*conj(e1) + 4*A1*D1*conj(c1)*conj(d1) + 4*B1*D1*conj(a1)*conj(e1) - 2*B1*D1*conj(b1)*conj(d1))^(1/2) - A1^2*conj(e1) + A1*B1*conj(d1) + A1*D1*conj(b1) - 2*B1*D1*conj(a1)))/(2*(conj(c1)*A1^2 - conj(b1)*A1*B1 + conj(a1)*B1^2)))/A1;

    X(2) = -(D1 - (B1*(A1^2*conj(e1) + A1*(A1^2*conj(e1)^2 + B1^2*conj(d1)^2 + D1^2*conj(b1)^2 - 4*B1^2*conj(a1) - 4*A1^2*conj(c1) - 4*D1^2*conj(a1)*conj(c1) + 4*A1*B1*conj(b1) - 2*A1*B1*conj(d1)*conj(e1) - 2*A1*D1*conj(b1)*conj(e1) + 4*A1*D1*conj(c1)*conj(d1) + 4*B1*D1*conj(a1)*conj(e1) - 2*B1*D1*conj(b1)*conj(d1))^(1/2) - A1*B1*conj(d1) - A1*D1*conj(b1) + 2*B1*D1*conj(a1)))/(2*(conj(c1)*A1^2 - conj(b1)*A1*B1 + conj(a1)*B1^2)))/A1;

    Y(1) = (A1*(A1^2*conj(e1)^2 + B1^2*conj(d1)^2 + D1^2*conj(b1)^2 - 4*B1^2*conj(a1) - 4*A1^2*conj(c1) - 4*D1^2*conj(a1)*conj(c1) + 4*A1*B1*conj(b1) - 2*A1*B1*conj(d1)*conj(e1) - 2*A1*D1*conj(b1)*conj(e1) + 4*A1*D1*conj(c1)*conj(d1) + 4*B1*D1*conj(a1)*conj(e1) - 2*B1*D1*conj(b1)*conj(d1))^(1/2) - A1^2*conj(e1) + A1*B1*conj(d1) + A1*D1*conj(b1) - 2*B1*D1*conj(a1))/(2*(conj(c1)*A1^2 - conj(b1)*A1*B1 + conj(a1)*B1^2));

    Y(2) = -(A1^2*conj(e1) + A1*(A1^2*conj(e1)^2 + B1^2*conj(d1)^2 + D1^2*conj(b1)^2 - 4*B1^2*conj(a1) - 4*A1^2*conj(c1) - 4*D1^2*conj(a1)*conj(c1) + 4*A1*B1*conj(b1) - 2*A1*B1*conj(d1)*conj(e1) - 2*A1*D1*conj(b1)*conj(e1) + 4*A1*D1*conj(c1)*conj(d1) + 4*B1*D1*conj(a1)*conj(e1) - 2*B1*D1*conj(b1)*conj(d1))^(1/2) - A1*B1*conj(d1) - A1*D1*conj(b1) + 2*B1*D1*conj(a1))/(2*(conj(c1)*A1^2 - conj(b1)*A1*B1 + conj(a1)*B1^2));

    end

    mi = [];
    mj = [];

    for i = 1:Image_num
        [X,Y] = crossLandC(l(i,:),line_image_p(:,:,i));
        mi(:,i) = [X(1) Y(1) 1]';
        mj(:,i) = [X(2) Y(2) 1]';    
    end
    

    %%
    %Step 3, obtain the intrinsic and mirror parameters%

    A = [];
    for i = 1:Image_num
        A = [A;real(mi(1,i)^2) real(mi(1,i)*mi(2,i)) real(mi(2,i)^2) real(mi(1,i)) real(mi(2,i)) 1;
            imag(mi(1,i)^2) imag(mi(1,i)*mi(2,i)) imag(mi(2,i)^2) imag(mi(1,i)) imag(mi(2,i)) 0];
    end
    
    
    [~,~,V] = svd(A); 
    c = V(:,end);
    c = c/c(6);
    C = [c(1)   c(2)/2 c(4)/2;
       c(2)/2 c(3)   c(5)/2;
       c(4)/2 c(5)/2 c(6)  ];
    K1 = inv(chol(C));
    K1 = K1/K1(3,3);
    est_T_K = Tp*K1; % get K
    

    
    est_T_l = 0;
    
    for i = 1:Image_num
        Qsi = est_T_K'*line_image(:,:,i)*est_T_K; 
        Qsi = Qsi/Qsi(3,3);
        li = sqrt(1-(Qsi(3,3)*Qsi(2,1))/(Qsi(3,2)*Qsi(3,1))); 
        est_T_l = est_T_l + li;
    end
    
    est_T_l = est_T_l/Image_num; % get l

end



%% Generate Random Central Catadioptric Mirror Contour
% Inputs:
%    d0 - Distance between mirror contour and sphere center
%    K - 3-by-3 - Intrinsic parameters
%    l - 1-by-1 - Mirror parameter
%
% Outputs:
%    mirror_contour - 3-by-3 - Mirror contour

function [mirror_contour] = generate_mirror_contour(d0,K,l)

  if d0 >= 1 && d0 < 0
        error('The range of d0 should be 0 <= d0 < 1.')
  end
    
  line_image = zeros(3,3);
  
  x = [0;0;1];
  n = x/norm(x);
  d = d0;
  Q = [(l^2-1)*n(1)^2+(d+l*n(3))^2  (l^2-1)*n(1)*n(2)                -(l*d+n(3))*n(1);
      (l^2-1)*n(1)*n(2)               (l^2-1)*(n(2))^2+(d+l*n(3))^2   -(l*d+n(3))*n(2);
      -(l*d+n(3))*n(1)                -(l*d+n(3))*n(2)                 d^2-n(3)^2      ];
  C = inv(K')*Q*inv(K);
  C = C/C(3,3);

  mirror_contour = C;
  
end

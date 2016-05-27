function feat = Feature( A, B )
%FEATURE Summary of this function goes here
%   Detailed explanation goes here
   feat = mean([corr2(A(:,:,1),B(:,:,1)),corr2(A(:,:,2),B(:,:,2)),corr2(A(:,:,3),B(:,:,3))]);
end


function xc2mat=xcorr2g(A,B)%(mat1,mat2)

[mA,nA] = size(A);
[mB,nB] = size(B);%
xc2mat = nan(mA+mB-1,nA+nB-1);
%isize = size(B,1); jsize = size(B,2);
[workmat,npadi,npadj] = pad_edges(A,B,1); %returns A padded by size of B with nans 
%for each i and j, choose the correct sub-matrix of A (size of B) to multiply with B
for i = 1:size(xc2mat,1)
    for j = 1:size(xc2mat,2) %get correct submatrix       
        subA = workmat(npadi+i-floor(mB):npadi+i-1 , npadj+j-floor(nB):npadj+j-1 );
        
        subAtimesB=subA .* B;
        
        notnanids = find(~isnan(subAtimesB));%normalized to the number of nontnan components (average)        
        
        n=length(notnanids);
        
        if n<20; xc2mat(i,j) = NaN; continue; end
        
        SAdotB =    sum( subAtimesB(notnanids)    ); 
        sumSA =     sum(       subA(notnanids)    );
        sumB =      sum(          B(notnanids)    );
        sumSAsqrd = sum(       subA(notnanids).^2 );
        sumBsqrd =  sum(          B(notnanids).^2 );
        
        xc2mat(i,j) = (     n*SAdotB - sumSA.*sumB) ./ ...
                      ( sqrt( n*sumSAsqrd - sumSA.^2 ) .* ...
                        sqrt( n* sumBsqrd - sumB .^2 ) ); 
    end % for j
end % for i

end

function out_mat=PointCorr(Rxx,RxxR)
    nan_mat=Rxx .* RxxR;
    notnan_inds = find(~isnan(nan_mat));  %normalized to the number of nontnan components (average)
    n=length(notnan_inds);
    if n < 2
        out_mat = NaN;
    end
    sigma_x_y =sum(nan_mat(notnan_inds));
    sigma_x =      sum(Rxx(notnan_inds));
    sigma_y =      sum(RxxR(notnan_inds));
    sigma_x2 = sum(Rxx(notnan_inds).^2);
    sigma_y2 = sum(RxxR(notnan_inds).^2);

    out_mat = (n*sigma_x_y - sigma_x.*sigma_y) ./ ...
        sqrt(n*sigma_x2-sigma_x.^2) ./ ...
        sqrt(n*sigma_y2-sigma_y.^2);
end

%         xc2mat(i,j) = (     n*sumSAdotB - sumSA.*sumB) ./ ...
%                       sqrt( n*sumSAsqrd - sumSA.^2   ) ./ ...
%                       sqrt( n* sumBsqrd - sumB .^2   )  ;

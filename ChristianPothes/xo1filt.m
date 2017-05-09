% S. Alex Kandel, 07/31/2006
% function im_out = xfilt(im_in,mask,filt_width)
% 
% im_in: 	input image, two-dimensional matrix (any size) of floating-point heights
% mask:	input mask(s); either:
%		1) a two-dimensional matrix of the same size of the image, with values of 1 for
%		   points to be included in fitting, and values of 0 everywhere else; or
%		2) a cell array, where each element is a mask, as defined in (1)
% filt_width:	frequency cut-off of the line-by-line high-pass filter
%		filt_width=0 for DC, filt_width=1 for (pi*x/N) frequency components, etc.
% additional:	'force' makes xfilt proceed even if the data is not flat.
% arguments 	(this is sometimes but not usually a good idea.)

function im_out = xo1filt(im_in)

    disp ("We enter xo1filt.m");

    filt_width = 1 ;
    im_out = [];

    num_masks = length(im_in);
						% this is used to define which points in the image will be fit
    [m,n] = size(im_in);
    disp ("The value of m*n is: "), disp (m*n);

    % FIXME
    mem = (2*filt_width+1)*m+num_masks-1 ;
    % disp ("The value of mem is: "), disp (mem);

    basis_mat = spalloc(m, n, m*n);  % sparse matrices needed
    disp ("size first spalloc: "), disp(length(basis_mat));

    basis_mat_mask = spalloc(m, n, m*n);    
    disp ("size second spalloc: "), disp(length(basis_mat_mask));

    disp ("enter for loop");
    for i=1:m
        line_basis_mat = ones(n,1);    			% always fit each scan line with a constant 
									% value    
        for j=1:filt_width					% add sin and cos components according to filt_width
            line_basis_mat = [line_basis_mat (sin((0:(n-1))/(n-1)*pi*j))' (cos((0:(n-1))/(n-1)*pi*j))'];
        end;    
        disp ("exit internal for loop");
        line_basis_mat_mask = line_basis_mat .* (ones(1,2*filt_width+1));   
        start_col = num_masks+(i-1)*(2*filt_width+1);
        start_row = (i-1)*n+1;
        basis_mat(start_row:(start_row+n-1),start_col:(start_col+2*filt_width)) = line_basis_mat;
        basis_mat_mask(start_row:(start_row+n-1),start_col:(start_col+2*filt_width)) = line_basis_mat_mask;
    end;
    disp ("exit for loop");
    total_fit = (basis_mat_mask\reshape(im_in',[],1))'*basis_mat';	% fit to the masked test matrix, and
										% multiply the result by the unmasked matrix
    im_out = im_in - reshape(total_fit',n,m)';
    disp ("The length of im_out is:"), disp (length(im_out))
    
return;

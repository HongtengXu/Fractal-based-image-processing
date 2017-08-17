%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% image superresolution and hallucination based on local fractal analysis
% on gradient domain
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;

for N=1
    
    filename=sprintf('0 (%d).png',N);
    
    in=imread(filename);   
    ori=in;
    
    [R,C,L]=size(ori);

    

    % parameters
    %up scaling
    uprate=4;%2,4,8; 
    % block size
    bs=max([4,2*uprate]);
    % blurring kernel
    gs=4*log2(uprate)+1;
    ker = fspecial('gaussian', [gs,gs], 0.4*uprate);
    % scaling for fractal analysis
    rad=3; 
    % weight of regularizer of gradient
    lambda1=0.5;
    % initial interpolation method
    inter='icbi'; % or 'bicubic'
    
    % for superresolution or enhancement
    fun='SR';
%     fun='EN';
    
     
    
    
    Lori=padarray(ori,[1,1],'symmetric','post');
    YUVori=rgb2ycbcr(Lori);
    if strcmp(inter,'icbi')==1
        H_YUVori =icbi(YUVori, log2(uprate), 8, 1, true);
        H_YUVori = H_YUVori(1:end-1,1:end-1, :);
    else
        H_YUVori=imresize(YUVori(1:end-1,1:end-1, :), uprate);
    end
    H_rgbori = ycbcr2rgb(H_YUVori);
    Yori=im2double(YUVori(:,:,1));
    
    
    
    
    tic;
    [outE, gout, gin, Alpha, Amply]=FractalGradientEnhance(Yori, uprate, ...
        ker, lambda1, bs, rad, inter, fun);
    


    H_YUVori(:,:,1)=im2uint8(outE);
    HoutE=ycbcr2rgb(H_YUVori);

    t6=toc;
   
    out1=imresize(ori,uprate,'nearest');
    
 
    filename2=sprintf('X%d_%d_%d.png'...
        ,uprate,N,round(t6));
    filename2=strcat(fun,'_',filename2);
    imwrite(im2uint8(HoutE),filename2);
   
end
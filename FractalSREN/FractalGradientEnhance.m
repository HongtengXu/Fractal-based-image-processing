%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Gradient Enhancement by Fractal analysis
%
% Potential application: superresolution, detail enhancement
%
%                                                   Hongteng Xu
%
% 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [out, gout, gin, Alpha, Amply]=FractalGradientEnhance(in, uprate, ker, ...
    lambda1, bs, rad, inter, fun)

% uprate:       >=2 superresolution(SR)
% ker:          default gaussian kernel for SR, size of kernel: (2*ks+1)^2
% lambda1:      the weight of regularization term of gradient
% bs:           block size
% rad:          the number of sample to fit fractal dimension
% inter:        the interpolation method: 'bilinear', 'bicubic', 'icbi'.
% fun:          'SR' or 'EN' which correspond to preserving invariance of
%               fractal in the gradient domain.



ks=(size(ker,1)-1)/2;


lowS=padarray(in,[2*ks+1,2*ks+1],'symmetric');
LS=padarray(lowS,[1,1],'symmetric','post');
        
DSx=LS(1:end-1,1:end-1)-LS(1:end-1,2:end);
DSy=LS(1:end-1,1:end-1)-LS(2:end,1:end-1);        
gin=sqrt(DSx.^2+DSy.^2);

if strcmp(inter,'icbi')==1        
    in =im2double( icbi(im2uint8(LS), log2(uprate), 8, 1, true) );
    in = in(1:end-1,1:end-1);
else
    in = im2double( imresize(im2uint8(lowS), uprate, inter) );
end
        
in_L=padarray(in,[1,1],'symmetric','post');
DLx=in_L(1:end-1,1:end-1)-in_L(1:end-1,2:end);
DLy=in_L(1:end-1,1:end-1)-in_L(2:end,1:end-1);
        
XX=DLx;
XX(DLx==0)=eps;
ratio=abs((DLy)./(XX));
clear XX;
        
Mx=2.*double(DLx>=0)-1;
My=2.*double(DLy>=0)-1;
gout=sqrt(DLx.^2+DLy.^2);



[dim_H,L_H]=fractalanalysis(gin, rad, 1 );
[dim_L,L_L]=fractalanalysis(gout, rad, uprate);
        

[R,C]=size(dim_L);
bound=ceil(bs);
dim_H=imresize(dim_H,[R,C]);
L_H=imresize(L_H,[R,C]);
        
Alpha=log2(uprate).*dim_H./dim_L;        
Amply=exp(L_H-Alpha.*L_L);


        
clear dim_H;
clear dim_L;
clear L_L;
clear L_H;

        
gout_L=padarray(gout,[bs*ceil(R/bs)-R, bs*ceil(C/bs)-C],'symmetric','post');
clear gout;
gout_L=padarray(gout_L, [bound,bound], 'symmetric');
        

        
a_L=padarray(Alpha,[bs*ceil(R/bs)-R, bs*ceil(C/bs)-C],'symmetric','post');
clear Alpha;
a_L=padarray(a_L, [bound,bound], 'symmetric');
        
a1_L=padarray(Amply,[bs*ceil(R/bs)-R, bs*ceil(C/bs)-C],'symmetric','post');
clear Amply;
a1_L=padarray(a1_L, [bound,bound], 'symmetric');
        
index=zeros(size(gout_L));
goutF=index;

        
for i=1:ceil(R/bs)
    for j=1:ceil(C/bs)

        pg=gout_L(1+(i-1)*bs:i*bs+2*bound, 1+(j-1)*bs:j*bs+2*bound);
        pa=a_L(1+(i-1)*bs:i*bs+2*bound, 1+(j-1)*bs:j*bs+2*bound);
        pa1=a1_L(1+(i-1)*bs:i*bs+2*bound, 1+(j-1)*bs:j*bs+2*bound);

                
        tmp=pg.^mean(pa(:)); 
        AA=norm(pg(:),2)/norm(tmp(:),2);
         
        if strcmp(fun,'EN')==1
            tmp=AA.*mean(pa1(:)).*tmp;
        else
            tmp=AA.*tmp;
        end
        tmp(isnan(tmp))=pg(isnan(tmp));
                
        goutF(1+(i-1)*bs:i*bs+2*bound, 1+(j-1)*bs:j*bs+2*bound)...
                        =goutF(1+(i-1)*bs:i*bs+2*bound, 1+(j-1)*bs:j*bs+2*bound)+tmp;
        index(1+(i-1)*bs:i*bs+2*bound, 1+(j-1)*bs:j*bs+2*bound)...
                        =index(1+(i-1)*bs:i*bs+2*bound, 1+(j-1)*bs:j*bs+2*bound)+1;

    end
end

gout=goutF./index;
        
gout=gout(bound+1:bound+R, bound+1:bound+C);
        
Alpha=a_L(bound+1:bound+R, bound+1:bound+C);
        
Amply=a1_L(bound+1:bound+R, bound+1:bound+C);


clear tmp;
clear goutF;
clear gout_L;
clear index;
clear a_L;
clear a1_L;
        
   
X=(gout.^2)./(ratio.^2+1);
Y=(gout.^2)-X;
DHx=Mx.*sqrt(X);
DHy=My.*sqrt(Y);
        
        
        
        
clear ratio;
clear Mx;
clear My
clear DSx;
clear DSy;
clear DLx;
clear DLy;
        
out=gradientCon_SR_fast(in, ker, DHx, DHy, lambda1);
        
clear DHx;
clear DHy;
out = out(1+uprate*(2*ks+1):end-uprate*(2*ks+2),1+uprate*(2*ks+1):end-uprate*(2*ks+2));
gout=gradient(out);
    


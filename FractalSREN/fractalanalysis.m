function [out1,out2]=fractalanalysis(in, rad, step)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% out1 is fractal dimension
% out2 is fractal length
%
% Hongteng Xu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

r=step.*(1:rad);
logr=log((2*r+1));
logr=logr(:);
A=[logr,ones(length(logr),1)];
out1=in;
out2=in;


tmp=zeros(size(in,1),size(in,2),rad);
for i=1:rad
    

    G=fspecial('gaussian',[1+2*r(i),1+2*r(i)],0.5*r(i));
    x=padarray(in,[r(i),r(i)],'symmetric');
    x=log(length(G(:))*conv2(255*sqrt(2).*x,G));
    tmp(:,:,i)=x(2*r(i)+1:end-2*r(i),2*r(i)+1:end-2*r(i));
    
end



for r=1:size(in,1)
    for c=1:size(in,2)
        logB=tmp(r,c,:);
        logB=logB(:);
        H=A\logB;

        out1(r,c)=H(1);
        out2(r,c)=H(2);
    end
end
clear tmp;

M=isnan(out1);
tmp=out1(~M);
out1(M)=mean(tmp(:));
out1=abs(out1);


M=isnan(out2);
tmp=out2(~M);
out2(M)=mean(tmp(:));
out2=abs(out2);


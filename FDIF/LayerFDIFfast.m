function [out, filters] = LayerFDIFfast( in, model, filters )




Sf = ones(2*model.wr+1)./((2*model.wr+1)^2);
in1 = padarray(in, [model.wr, model.wr], 'symmetric');
in1 = conv2(in1, Sf, 'valid');
in2 = in.^model.scale;
in3 = padarray(in2, [model.wr, model.wr], 'symmetric');
in3 = conv2(in3, Sf, 'valid');
in3 = 1./in3;
in3(isinf(in3))=0;
in = in1.*in2.*in3;
patches = im2col(padarray(in,[model.wr,model.wr],'symmetric'),...
    [2*model.wr+1, 2*model.wr+1], 'sliding');
out = size(in);

if isempty(filters)    
    filters = zeros(size(patches));
    
    filter = model.weight(:,:,1);
    imgblur = padarray(in, [1,1], 'symmetric');
    gf = fspecial('gaussian', [3,3], 0.5);
    imgblur = conv2(imgblur, gf, 'valid');
    [dx, dy]=gradient(imgblur);

    dx = padarray(dx, [model.wr, model.wr], 'symmetric');
    dy = padarray(dy, [model.wr, model.wr], 'symmetric');

    for c = 1:size(in,2)
        for r = 1:size(in,1)
            Dx = dx(r:r+2*model.wr, c:c+2*model.wr);
            Dy = dy(r:r+2*model.wr, c:c+2*model.wr);

            mat = [Dx(:),Dy(:)];
            [V, ~] = eig(mat'*mat);

            deg = atand(V(1,2)/V(2,2));

            flt = imrotate(filter, deg, 'nearest', 'crop');
            flt = flt./sum(flt(:));

            filters(:,(c-1)*size(in,1)+r) = flt(:);
            out(r,c) = flt(:)'*patches(:, (c-1)*size(in,1)+r);
        end
    end
else
    out = sum( filters.*patches );
    out = reshape(out, [size(in,1), size(in,2)]);    
end
out=out./max(out(:));

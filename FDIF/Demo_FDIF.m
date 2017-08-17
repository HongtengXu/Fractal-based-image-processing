clear

% parameters of model
model.wr = 4; % radius of patch
model.ln = 5; % number of layer
model.fn = 30; % number of filter
model.scale = 2; % intensity of fractal enhancement ?for curve detection?
filter = zeros(model.wr*2+1, model.wr*2+1, model.fn);
filter(model.wr+1,:,1) = 1./(model.wr*2+1);
% for i=1:model.fn-1
%     deg = i*(180/model.fn);
%     if deg==45
%         tmp = fliplr(eye(model.wr*2+1));
%     else if deg==135
%             tmp = eye(model.wr*2+1);
%         else       
%             tmp = imrotate(filter(:,:,1), deg, 'bilinear', 'crop');
%         end
%     end
%     tmp = tmp./sum(tmp(:));
%     filter(:,:,i+1)=tmp;
% end
model.weight = filter; % structural filter bank

for i = 1:9
    ori = imread(sprintf('ori_%d.png', i));
    if size(ori,3)>1
        ori = rgb2gray(ori);
    end
    ori = imresize(ori, 250/size(ori,2));
    ori = im2double(ori);
    
    
    M = mean(ori(:));
    tmp = ori;
    tic
    figure
    subplot(1,model.ln+1,1)
    imshow(ori)
    for l=1:model.ln      
        if l<=1
            [tmp, filters] = LayerFDIFfast( tmp, model, [] );
        else
            [tmp, filters] = LayerFDIFfast( tmp, model, filters);
        end
        tmp = M/mean(tmp(:))*tmp;
        
        subplot(1, model.ln+1, l+1)
        imshow(tmp)
        imwrite(tmp, sprintf('new%d_%d.png', i,l))
        toc;
    end
  
end
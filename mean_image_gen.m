function [mean_image] = mean_image_gen(series,dim,n,file)

    mean_image = zeros(dim(1),dim(2));

    % create mean image from average of 'series'
    for j = 1:dim(2)
        for i = 1:dim(1)
            mean_image(i,j) = nanmean(series(i,j,:));
            if isnan(mean_image(i,j))
                mean_image(i,j) = 0;
            else if (mean_image(i,j) < 0)
                    mean_image(i,j) = 0;
                end
            end
        end
    end
    figure(222)
    %imagesc(log10(mean_image))
    h=imagesc(mean_image);
    colorbar
    
    filename = strrep(file,'_lpd-targ.fits','');
    filename = strrep(filename,'ktwo','outputs/pipeout_ktwo');
    %filename = strcat(filename,'_',ms_str);
    %filename = strcat(filename,'_target', num2str(targetNum))
    filename = strcat(filename,'_image.jpg');
    saveas(h,filename,'jpg');
        
    end
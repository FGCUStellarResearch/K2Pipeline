function bkg = fit_background2(apdim,series)

    for count=1:apdim
        temp = series(:,:,count);
        temp = temp(:);
        temp = temp(temp>0);
        temp = sort(temp);
        bkg(count) = median(temp(1:ceil(numel(temp)/5))); %median of faintest 20% of pixels
    end  
    
end
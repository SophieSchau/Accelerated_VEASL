function [ spat_corr ] = calculate_spat_corr( GT, Recon, temporal_concat, mask,show_fig )
%CALCULATE_SPAT_CORR
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<5
    show_fig =0;
end
if nargin<4
    mask = ones(size(GT,1), size(GT,2), size(GT,3), size(GT,5));
end

if show_fig
    figure
end

if size(GT,5) == 2 %if nonVE
    temp(:,:,:,:,1:3) = repmat(GT(:,:,:,:,1),1,1,1,1,3);
    temp(:,:,:,:,4) = GT(:,:,:,:,2);
    
    GT = temp;
    
    temp(:,:,:,:,1:3) = repmat(Recon(:,:,:,:,1),1,1,1,1,3);
    temp(:,:,:,:,4) = Recon(:,:,:,:,2);
    
    Recon = temp;
end

colors = [1 1 0; 0 1 1; 1 0 1; 1 1 1];


if temporal_concat % average across temporal dimension
    spat_corr = zeros(size(Recon, 5),1);
    for v = 1:size(Recon, 5) %vessel
                Reconlist = Recon(:,:,:,:,v);
                Reconlist(repmat(~mask(:,:,:,v),1,1,1,size(Recon,4))) = [];
                Reconlist = Reconlist(:);
                
                GTlist = GT(:,:,:,:,v);
                GTlist(repmat(~mask(:,:,:,v),1,1,1,size(Recon,4))) = [];
                GTlist = GTlist(:);
                
                spat_corr(v) = corr(abs(Reconlist), abs(GTlist));
                
                if show_fig
                    if v < size(Recon, 5)
                        scatter(abs(GTlist),abs(Reconlist),20,colors(v,:))
                        hold on
                    end
                    legend('RICA', 'LICA', 'BA')
                end
                clear('Reconlist', 'GTlist')

    end
else % calculate for each time frame separately
    spat_corr = zeros(size(Recon, 5),size(Recon, 4));
    
    for v = 1:size(Recon, 5) %vessel
        for t = 1:size(Recon, 4) %time
                Reconlist = Recon(:,:,:,t,v);
                Reconlist(~mask(:,:,:,v)) = [];
                Reconlist = Reconlist(:);
                
                GTlist = GT(:,:,:,t,v);
                GTlist(~mask(:,:,:,v)) = [];
                GTlist = GTlist(:);
                
                spat_corr(v,t) = corr(abs(Reconlist), abs(GTlist));
                
                if show_fig
                    if v < size(Recon, 5)
                        scatter(abs(GTlist),abs(Reconlist),20,colors(v,:))
                        hold on 
                    end
                    legend('RICA', 'LICA', 'BA')
                end

        end
    end
end

spat_corr(isnan(spat_corr)) = 0;

end

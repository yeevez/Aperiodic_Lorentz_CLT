function [validation_vec,validation_vec_2] = validate(positions,scatterers,generators,radius)
    %This function will take the positions and
    %scatterer centers from a trial and calculate the distance between each
    %pair. it counts how many bounces were not within the radius (plus or
    %minus a suitable error margin) from their corresponding center. Then,
    %it calculates if the point that we concluded to be the center of the
    %scatterer we bounce off of is actually one of the scatterer centers we
    %originally calculated - it does this rather inefficiently, so don't
    %run it on a large sample
    %we assume a single trial as the input for simplicity
    nsamples = size(positions,2);
    validation_vec = zeros(1,nsamples);
    for i=1:nsamples
        position = positions(:,i);
        scatterer = scatterers(:,i);
        norm(position-scatterer,2)
        if abs(norm(position-scatterer,2) - radius) > 10e-6
            validation_vec(1,i) = 1;
            sprintf("bounce %d is too far away from its scatterer",i)
        end
    end

    sprintf("%d out of %d total bounces were too far away from their scatterer", sum(validation_vec),nsamples)
    
    validation_vec_2 = zeros(1,nsamples);
    matrix = [1,2,1,1,2;2,5,4,3,5;1,4,6,5,6;1,3,5,7,9;2,5,6,9,14];
    [v,~]=eig(matrix);
    outdim = 3;
    r=rotation_matrix(v,outdim);
    pgrid=20;
    window=0.5;
    for i = 1:nsamples
        scatterer = scatterers(:,i);
        for j=1:size(generators,2)
            generator = generators(:,j);
            [sp,~,~] = scatterer_positions(r,window,pgrid,generator,outdim,false,radius);
            for k=1:size(sp,2)
                center = sp(:,k);
                if abs(norm(center-scatterer)) <= 10e-6
                    validation_vec_2(1,i) = 1;
                end
            end
        end
    end

    sprintf("%d out of %d total scatterer centers were reproduced during validation", sum(validation_vec_2),nsamples)
    end
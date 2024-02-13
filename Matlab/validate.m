function validation_vec = validate(positions,scatterers,radius)
    %This function will take the positions and
    %scatterer centers from a trial and calculate the distance between each
    %pair. it counts how many bounces were not within the radius (plus or
    %minus a suitable error margin) from their corresponding center. Then,
    %it calculates if the point that we concluded to be the center of the 
    %we assume a single trial as the input for simplicity
    nsamples = size(positions,2);
    validation_vec = zeros(1,nsamples);
    for i=1:nsamples
        position = positions(:,i);
        scatterer = scatterers(:,i);
        norm(position-scatterer,2)
        if abs(norm(position-scatterer,2) - radius) > 10e-5
            validation_vec(1,i) = 1;
        else
            sprintf("bounce %d is too far away from its scatterer",i)
        end
    end

    sprintf("%d out of %d total bounces were too far away from their scatterer", sum(validation_vec),nsamples)
function [boolean, posOfColl, newV] = inCube(posE, centerCube, vel, side)
    newCenter = centerCube;
    x1 = newCenter(0) - (side/2);
    x2 = newCenter(0) + (side/2);
    y1 = newCenter(1) - (side/2);
    y2 = newCenter(1) + (side/2);
    if ( posE(0) < x1 )
        if (vel(0) > 0)
            possibleCollisionPoint = vel(1)*( x1 - posE(0))/vel(0);
            if (possibleCollisionPoint  + posE(1)> y1 && possibleCollisionPoint  + posE(1)< y2)
                boolean = true;
                newV = [-1*vel(0), vel(1)];
                posOfColl = [x1, possibleCollisionPoint + posE(1)];
                return
            end
        end
    end

    if ( posE(0) > x2 )
        if (vel(0) < 0)
            possibleCollisionPoint = vel(1)*( posE(0) - x2)/vel(0);
            if (possibleCollisionPoint + posE(1)> y1 && possibleCollisionPoint + posE(1)< y2)
                boolean = true;
                newV = [-1*vel(0), vel(1)];
                newpos = [ x2, possibleCollisionPoint + posE(1) ];
                return
            end
        end
    end

    if ( posE(1) < y1 )
        if(vel(1) > 0)
            possibleCollisionPoint = vel(0)*( y1 - posE(1))/vel(1);
            if (possibleCollisionPoint + posE(0)> x1 && possibleCollisionPoint + posE(0) < x2)
                boolean = true;
                newV = [vel(0), -1*vel(1)];
                newpos = [ possibleCollisionPoint + posE(0), y1 ];
                return
            end
        end
    end

    if ( posE(1) > y2)
        if (vel(1) < 0)
            possibleCollisionPoint = vel(0)*( posE(1) - y2)/vel(1);
            if (possibleCollisionPoint + posE(0) > x1 && possibleCollisionPoint + posE(0) < x2)
                boolean = true;
                newV = [vel(0), -1*vel(1)];
                newpos = [ possibleCollisionPoint + posE(0), y2 ];
                return
            end
        end
    end
    boolean = false;
    newV = NULL;
    newpos = NULL;
end
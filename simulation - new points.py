"""
This template .py file contains all the helper functions necessary to run the simulation, and can be run
as is to generate collisions under the same parameters that we have been doing so far, but is designed to 
be copied and some parameters changed. For example, you could change the matrix that generates the aperiodic
points, the radius of the points, the shape of the scatterers, the number of collisions, the number of
collisions whose positions we actually save in each sample, or more broadly generalize or modify the simulation
 in other ways. Make sure to adjust the names of the samples being saved to reflect the specific parameters
of this simulation run
"""

import numpy as np

##code to generate irrational planes from hyperbolic matrices.
def hyperbolic_matrix(dim,rand = False):
    #returns a random, dim dimensional symmetric hyperbolic matrix
    U = np.eye(dim)
    if not rand:
        for i in range(dim):
            for j in range(dim):
                if i < j:
                    U[i][j] = 1

    else:
        for i in range(dim):
            for j in range(dim):
                if i < j:
                    U[i][j] = np.random.randint(1,high=5)
    L = U.T
    A = L @ U
    return A

def irrational_plane(A, dim):
    #returns a dim-1 dimensional eigen subspace of a dim dimensional hyperbolic matrix
    subspace_dim = dim - 1
    eigenvals,eigenvecs = np.linalg.eig(A)
    eigenvectors = []
    for i in range(dim):
        eigenvectors.append(eigenvecs[:,i])
    return eigenvectors[:subspace_dim]

##code to project aperiodic points onto a plane
def point_projection(vectors,point):
    #returns the 2 dimensional representation of the point projection,
    #ie the coefficients multiplied by the eigenvectors to create the ortho projection
    projection = np.zeros(len(vectors))
    for i in range(len(vectors)):
        projection[i] = point.dot(vectors[i]) #coefficient
    return projection


def aperiodic_points(R,dim,eig,borderSize,position):
    #x1<x2, y1<y1
    #input a range of integers to iterate over R, a dimension, a hyperbolic matrix, and a window size - get in return a set of aperiodic points within the range
    points = []
    patch_center = sum([position[i]*eig[i] for i in range(dim-1)]) #3d representation of position on plane
    patch_center = np.floor(patch_center).astype(int)

    r = np.eye(dim)
    x = patch_center[0]
    y = patch_center[1]
    
    for i in range(x-R,x+R):
        for j in range(y-R,y+R):
            b = (j*eig[0][0] - i*eig[0][1])/(eig[0][0]*eig[1][1]-eig[1][0]*eig[0][1])
            a = (i-b*eig[1][0])/eig[0][0]
            z = a*eig[0][2] + b*eig[1][2]
            #z = (i*eig[0] + j*eig[1])[2]
            borderFloor = np.floor((z-borderSize))
            borderRoof = np.floor((z+borderSize+1))
            for point in range( borderFloor.astype(int), borderRoof.astype(int)):
                if point > (z-borderSize) and point < (z+borderSize+1):
                    projected = point_projection(eig, np.array([i,j,point]))
                    points.append(projected)
    return points

#Return closest point with respect to x to thing in xPosition
def indexClosestPoint(listOfVectors, xPosition):
    size = len(listOfVectors)
    return auxIndexClosestPoint(listOfVectors, 0, size-1, xPosition)
def auxIndexClosestPoint(listOfVectors, nstart, nend, xPosition):
#    if (mend==nend and mstart==nstart):
#        print(listOfVectors[nstart],listOfVectors[nend], nstart, nend, xPosition, "\n")
    size = nend - nstart + 1
    if (size == 1):
        return nend
    if (size == 2):
        if ( (listOfVectors[nend][0]-xPosition)*(listOfVectors[nend][0]-xPosition) < (listOfVectors[nstart][0]-xPosition)*(listOfVectors[nstart][0]-xPosition)):
            return nend
        else:
            return nstart
    halfSize = (nstart + np.floor(size/2)).astype(int)
    if((listOfVectors[halfSize][0]) == xPosition):
        return halfSize
    if((listOfVectors[halfSize][0]) < xPosition):
        return auxIndexClosestPoint(listOfVectors, halfSize, nend, xPosition)
    return auxIndexClosestPoint(listOfVectors, nstart, halfSize, xPosition)


#Returns list of centers with which a colision occurs and corresponding list of projections of the center to a line spaned by the velocity
def posHitInArea(listOfVectors, position, velocity, radius,p):
    cp = indexClosestPoint(listOfVectors, position[0])
    it = 0
    length = len(listOfVectors)
    collisionList = []
    projectionList = []
    proj_matrix = np.outer(velocity,velocity)
    iter_direction = np.sign(velocity[0]).astype(int)
    radiusSquared = radius**2
    radiusDoubled = 2*radius

    if abs(velocity[1] + 1) < 0.02:
        k = 2
    else:
        k = 1

    if (velocity[0] == 0):
        itB = 0
        itA = 0
        before = listOfVectors[cp+itB]
        after = listOfVectors[cp+itA]
        while(before[0]-position[0]>-radius):
            if (cp+itB==length):
                break
            before = listOfVectors[cp+itB]
            currentVector = before
            newPosCenterSphere =  currentVector - position
            #We check that the sphere is close enough to collide and if so append
            projectionOfCenter = proj_matrix @ newPosCenterSphere
            centerToProjection = projectionOfCenter - newPosCenterSphere
            CTPLengthSQ = centerToProjection.dot(centerToProjection)
            if (CTPLengthSQ <= radiusSquared):
                collisionList.append( currentVector)
                projectionList.append(projectionOfCenter+position)
            itB+=1
        while(after[0]-position[0]<radius):
            if (cp+itA<0):
                break
            after = listOfVectors[cp+itA]
            currentVector = after
            newPosCenterSphere =  currentVector - position
            #We check that the sphere is close enough to collide and if so append
            projectionOfCenter = proj_matrix @ newPosCenterSphere
            centerToProjection = projectionOfCenter - newPosCenterSphere
            CTPLengthSQ = centerToProjection.dot(centerToProjection)
            if (CTPLengthSQ <= radiusSquared):
                collisionList.append( currentVector)
                projectionList.append(projectionOfCenter+position)
            itA-=1
            after = listOfVectors[cp+itA]
        return collisionList,projectionList
    while (cp+it<length and cp+it>=0):
        currentVector = listOfVectors[cp+it]
        #We move the plane so the electrton is in the center and this moves the circles to their according position
        newPosCenterSphere =  currentVector - position
        #We check that the sphere is close enough to collide and if so append
        projectionOfCenter = proj_matrix @ newPosCenterSphere
        centerToProjection = projectionOfCenter - newPosCenterSphere
        norm = np.linalg.norm(centerToProjection,ord=p)
        if (norm < radius):
            collisionList.append( currentVector)
            projectionList.append(projectionOfCenter+position)
        if (len(projectionList)>k):
            if iter_direction*(currentVector[0]-projectionList[k][0])>radiusDoubled:
                return collisionList,projectionList
        it += iter_direction
    return collisionList,projectionList

##code to simulate the path of a particle through an aperiodic medium
def simulate_diffusion(radius,max_bounces,R,bordersize,p=2):
    #primary function to simulate particle collisions in media. input parameters: list of aperiodic points,
    #radius of obstacles, number of collisions to simulate to, and range of values to initialize position and
    #velocity within. This function continuously checks for the next particle collision until the stopping
    #criteria have been met, and returns a list of the position and velocity of the particle at each collision
    #it relies on a number of helper functions:
    n_bounces = 0
    dim = 3
    A = hyperbolic_matrix(dim)
    eigenvecs = irrational_plane(A,dim)
    results = np.zeros(((max_bounces//1000)+ 1,2))
    position = np.random.uniform(-1000,1000,2)
    velocity = np.random.uniform(-R,R,2)
    velocity = velocity/np.linalg.norm(velocity)
    largest_flight_time = 0

    points = sortVecX(aperiodic_points(R,dim,eigenvecs,bordersize,position))
    length = len(points)
    
    for point in points:
        #make sure we don't initialize position inside a scatterer
        if np.linalg.norm((position-point)) < radius:
            diff = (position - point)
            scalar = radius - np.linalg.norm(diff) + 0.01
            position += scalar*diff
            break

    results[0] = position
    
    while n_bounces < max_bounces:
        n_bounces += 1
        if len(points) > length:
            points = sortVecX(aperiodic_points(R,dim,eigenvecs,bordersize,position))
        (position,old_position,velocity,results,points) = next_collision(position,velocity,radius,results,points,R,dim,eigenvecs,bordersize,p)
        if np.linalg.norm(position-old_position) > largest_flight_time:
            largest_flight_time = np.linalg.norm(position-old_position)
        if n_bounces%1000 == 0:
            results[n_bounces/1000] = position
        #print("num bounces: ", n_bounces)
    return (results,largest_flight_time)

def next_collision(position,velocity,radius,results,points,R,dim,eig,bordersize,p):
        #given a position, velocity and number of bounces, finds the next collision in the simulation
        collision_point,obj_center = find_collision_point(position,velocity,points,radius,p)
        
        i = 0
        while len(collision_point) == 0:
            points = sortVecX(aperiodic_points(R+i,dim,eig,bordersize,position))
            collision_point,obj_center = find_collision_point(position,velocity,points,radius,p)
            i += 2
        '''
        tan = np.array([collision_point[0]**(p-1),collision_point[1]**(p-1)])
        tan /= np.linalg.norm(tan,ord=p)
        #perp = (obj_center-collision_point)/radius
        if tan[0] != 0:
            perp = np.array([-(tan[1]/tan[0]),1])
        else:
            perp = np.array([1,-(tan[0]/tan[1])])
        perp /= np.linalg.norm(perp, ord=p)


        velocity_new = -(velocity.dot(tan)*tan - velocity.dot(perp)*perp)'''
        perp = (obj_center-collision_point)/radius
        if perp[0] != 0:
            ortho_perp = np.array([-(perp[1]/perp[0]),1])
        else:
            ortho_perp = np.array([1,-(perp[0]/perp[1])])
        ortho_perp /= np.linalg.norm(ortho_perp)

        velocity_new = -(velocity.dot(perp)*perp - velocity.dot(ortho_perp)*ortho_perp)
        position_new = collision_point

        velocity_new = velocity_new/np.linalg.norm(velocity_new)
        
        return(position_new,position,velocity_new,results,points)

def sortVecX(vectors):
    #sorts list of points by x value in ascending order
    return sorted(vectors, key=lambda v: v[0])

def sortVecY(vectors):
    #sorts list of points by y value in ascending order
    return sorted(vectors, key=lambda v: v[1])

def find_collision_point(position, velocity,points,radius,p):
    #finds the next point of collision for the particle. calls the posHitInArea function to
    #generate list of possible object centers resulting in collisions, calculates collision points,
    #and returns the closest one to the current position.
    obj_centers,projections = posHitInArea(points,position,velocity,radius,p)
    centers_new = []
    possible_collisions = []
    radius_squared = radius**2
    for i in range(len(obj_centers)):
        center = obj_centers[i]
        projection = projections[i]
        a = np.sqrt(np.abs(radius**2-(center-projection).dot((center-projection))))
        collision_point_1 = projection + a*velocity
        collision_point_2 = projection - a*velocity
        if not (np.isclose(position,collision_point_1)).all() and not (np.isclose(position,collision_point_2)).all():
            if np.sign(collision_point_1[0] - position[0]) == np.sign(velocity[0]) and np.sign(collision_point_1[1] - position[1]) == np.sign(velocity[1]):
                centers_new.append(center)
                centers_new.append(center)
                possible_collisions.append(collision_point_1)
                possible_collisions.append(collision_point_2)
                
        
    if not possible_collisions:
        collision = np.array([])
        center = np.array([])
        return(collision,center)
    collision = possible_collisions[0]
    center = centers_new[0]
    for i in range(len(possible_collisions)):
        if (position-possible_collisions[i]).dot((position-possible_collisions[i])) < (position-collision).dot((position-collision)):
            collision = possible_collisions[i]
            center = centers_new[i]
    return (collision,center)


#CODE TO ACTUALLY COLLECT DATA: SET N_SAMPLES = NUMBER OF DESIRED SAMPLES, n_existing_samples TO NUMBER OF SAMPLES YOU ALREADY HAVE,
#max_bounces TO NUMBER OF COLLISIONS YOU WISH TO SIMULATE. leave radius, R, bordersize parameters unchanged.

N_SAMPLES = 10000
n_existing_samples = 0
radius = 0.3
max_bounces = 5000000
R = 3
bordersize = 1.5244586697611553


for i in range(n_existing_samples+1,N_SAMPLES):
    filename = '2d_circular_obstacle_5milcollisions_matrixA_sample_' + str(i)
    flight_filename = '2d_circular_obstacle_5milcollisions_matrixA_sample_' + str(i) + '_flight_time'
    results,flight_time = simulate_diffusion(radius,max_bounces,R,bordersize)
    np.save(filename,results,allow_pickle = True)
    np.save(flight_filename,np.array([flight_time]))
    print("sample " + str(i) + " completed")
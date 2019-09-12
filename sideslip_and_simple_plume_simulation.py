from __future__ import print_function
import numpy as np
import os, sys
import matplotlib
import matplotlib.pyplot as plt
import random
import json
from numpy import *
import operator

trap_number = 10
trap_angle_list = np.linspace(0,np.pi*2, trap_number, endpoint=False)
trap_names = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
distance_to_trap = 1000 # in meters
trap_point_list = np.array([[np.cos(a)*distance_to_trap, np.sin(a)*distance_to_trap] for a in trap_angle_list])

tracking_prob_cutoff = 0.1

plume_angular_width = 0*np.pi/180. # in radians

fly_number = 360

preferred_groundspeed_along_body_axis = 1.0 #meters per second.
forward_airspeed_limit = 1.8 #meters per second.
reverse_airspeed_limit = -0.2 #meters per second.

wind_direction_list = [3*np.pi/2.]
wind_speed_list = [0.0]


weight_of_perpendicular_slip = 1.0
fly_headings = list(np.linspace(0, np.pi*2, fly_number, endpoint = False))

def calc_plume_track_prob(distance_from_plume, space_constant = 100):
    """
    Need to extensively discuss this with Annie; she did a lot of work to reduce the dimensionality of the plume tracking problem and I think she ended up with something that was sort of polynomial function relative to distance from plume source.
    """
    return 1.0*(np.exp(-1*distance_from_plume/float(space_constant)))

def scalar_projection(a,b): #projects a onto b, yields magnitude in direction of b
    return np.dot(a,b) / np.linalg.norm(b)

def vector_projection(a,b): #projects a onto b, yields scaled vector in direction of b
    return b * np.dot(a,b) / np.dot(b,b)

def calculate_trajectory_vector (wind_speed, wind_dir,
                            fly_heading, preferred_groundspeed_along_body_axis,
                            forward_airspeed_limit, reverse_airspeed_limit,
                            weight_of_perpendicular_slip, ax_handle):
        wind_vector = np.array([np.cos(wind_dir)*wind_speed, np.sin(wind_dir)*wind_speed]) #length of this vector is in units of m/s
        fly_heading_unit_vector = np.array([np.cos(fly_heading), np.sin(fly_heading)]) # length of this vector is not meaningful.
        parallel_wind_vector = vector_projection(wind_vector, fly_heading_unit_vector) #length of this vector is in units of m/s, negative values allowed

        perpendicular_wind_vector = wind_vector - parallel_wind_vector #length of this vector is in units of m/s

        attempted_groundspeed_vector_along_body_axis = fly_heading_unit_vector * preferred_groundspeed_along_body_axis
        attempted_airspeed_vector_along_body_axis = attempted_groundspeed_vector_along_body_axis - parallel_wind_vector
        a = scalar_projection(attempted_airspeed_vector_along_body_axis,fly_heading_unit_vector)

        if reverse_airspeed_limit <= a <= forward_airspeed_limit:
            groundspeed_vector_along_body_axis = attempted_groundspeed_vector_along_body_axis
        elif a > forward_airspeed_limit: # fly cannot thrust enough to achieve preferred groundspeed along body axis
            actual_airspeed_vector_along_body_axis = forward_airspeed_limit*fly_heading_unit_vector
            groundspeed_vector_along_body_axis = actual_airspeed_vector_along_body_axis + parallel_wind_vector
            if scalar_projection(groundspeed_vector_along_body_axis, fly_heading_unit_vector) <0: #fly would be pushed backwards along her body axis
                print ('maintaining this heading is not possible in this wind')
                return (None, None)
        elif a < reverse_airspeed_limit: # fly cannot brake enough to achieve preferred groundspeed along body axis
            print ('fly cannot brake enough')
            actual_airspeed_vector_along_body_axis = reverse_airspeed_limit*fly_heading_unit_vector
            groundspeed_vector_along_body_axis = actual_airspeed_vector_along_body_axis + parallel_wind_vector
        trajectory_vector = groundspeed_vector_along_body_axis + perpendicular_wind_vector*weight_of_perpendicular_slip
        return [wind_vector, trajectory_vector, perpendicular_wind_vector, parallel_wind_vector, fly_heading_unit_vector]


def check_if_point_is_on_ray (point, ray): # ray will look like [[x,y],[x2,y2]]
    ray_to_test = np.subtract(point,ray[0])
    ray_from_origin = np.subtract(ray[1], ray[0])
    #magnitude = scalar_projection(ray_to_test,ray_from_origin)
    sign = np.sign(np.dot(ray_to_test,ray_from_origin))
    # print ()
    # print (ray_from_origin)
    # print (ray_to_test)
    # print (magnitude)
    if sign <0:
        return False
    else:
        return True

def get_intersect(a1, a2, b1, b2):
    """
    Returns the point of intersection of the lines passing through a2,a1 and b2,b1.
    a1: [x, y] a point on the first line
    a2: [x, y] another point on the first line
    b1: [x, y] a point on the second line
    b2: [x, y] another point on the second line
    """
    s = np.vstack([a1,a2,b1,b2])        # s for stacked
    h = np.hstack((s, np.ones((4, 1)))) # h for homogeneous
    l1 = np.cross(h[0], h[1])           # get first line
    l2 = np.cross(h[2], h[3])           # get second line
    x, y, z = np.cross(l1, l2)          # point of intersection

    if z == 0:                          # lines are parallel
        return (float('inf'), float('inf'))
    trajectory_ray = np.array([a1,a2])
    plume_ray = np.array([b1,b2])
    intersection_point = np.array([x,y])
    # print ('trajectory ray: '+ str(trajectory_ray))
    # print ('plume ray: ' +str(plume_ray))
    # print ('intersection point: '+str(intersection_point))
    # print ()
    trajectory_bool = check_if_point_is_on_ray(intersection_point,trajectory_ray)
    plume_bool = check_if_point_is_on_ray(intersection_point, plume_ray)
    if trajectory_bool and plume_bool:
        print ('both booleans true')
        return (x/z, y/z)
    else:
        print ('at least one of the booleans was false')
        return (float('inf'),float('inf'))

def ask_if_fly_intersects_a_plume(trajectory, trap_point,plume_angular_width, wind_direction):
    """
    trajectory is a single (x,y) point; defined with respect to release site at origin
    trap_point is a single (x,y) point
    plume_angular_width is in radians
    wind_direction is in radians
    """
    plume_angle_1 = (wind_direction + plume_angular_width/2. + 2*np.pi) % (2*np.pi)
    plume_angle_2 = (wind_direction - plume_angular_width/2. + 2*np.pi) % (2*np.pi)
    plume_unit_vector_1 = [np.cos(plume_angle_1)+trap_point[0], np.sin(plume_angle_1)+trap_point[1]]
    plume_unit_vector_2 = [np.cos(plume_angle_2)+trap_point[0], np.sin(plume_angle_2)+trap_point[1]]

    # print (trap_point)
    # print ('plume angular width: ' +str(plume_angular_width*180/np.pi))
    # print ('wind direction:      ' +str(wind_direction*180/np.pi))
    # print ('plume angle 1:       ' + str(plume_angle_1*180/np.pi))
    # print ('plume angle 2:       ' + str(plume_angle_2*180/np.pi))
    # print ('plume unit vector 1: ' +str(plume_unit_vector_1))
    # print ('plume unit vector 2: ' +str(plume_unit_vector_2))
    # print ()

    intersection_1 =  get_intersect((0,0) , trajectory, trap_point, plume_unit_vector_1 )
    distance_1 = np.sqrt((intersection_1[0]-trap_point[0])**2 + (intersection_1[1]-trap_point[1])**2)
    print (distance_1)
    intersection_2 =  get_intersect((0,0) , trajectory, trap_point, plume_unit_vector_2 )
    distance_2 = np.sqrt((intersection_2[0]-trap_point[0])**2 + (intersection_2[1]-trap_point[1])**2)
    print (distance_2)
    return np.min([distance_1,distance_2])
    #get_intersect((0,0) , trajectory, trap_point, plume_unit_vector_2 )

fig = plt.figure(figsize=(8,8))
ax_handle = plt.subplot(1,1,1)

dictionary = {}
for wind_speed in wind_speed_list:
    dictionary[str(wind_speed)] = {}
    for wind_dir in wind_direction_list:
        experiment_dictionary = {'trajectories':[],'wind vectors': [], 'perp_vectors':[], 'par_vectors':[], 'heading_unit_vectors':[],'groundspeeds along trajectories':[], 'windspeeds along trajectories':[], 'trajectory_vectors':[]}

        for fly_heading in fly_headings:
            return_list = calculate_trajectory_vector(wind_speed, wind_dir,
                                        fly_heading, preferred_groundspeed_along_body_axis,
                                        forward_airspeed_limit, reverse_airspeed_limit,
                                        weight_of_perpendicular_slip, ax_handle)
            if return_list[0]==None:
                continue
            wind_vector = return_list[0]
            trajectory_vector = return_list[1]

            trajectory_angle = np.arctan2(trajectory_vector[1],trajectory_vector[0]) # trajectory angle in radians
            trajectory_angle = (trajectory_angle+2*np.pi)%(2*np.pi) #recast from 0 to 2pi
            windspeed_along_trajectory = scalar_projection(wind_vector, trajectory_vector)
            groundspeed_along_trajectory = np.linalg.norm(trajectory_vector)
            experiment_dictionary['trajectories'].append(trajectory_angle)
            experiment_dictionary['groundspeeds along trajectories'].append(groundspeed_along_trajectory)
            experiment_dictionary['windspeeds along trajectories'].append(windspeed_along_trajectory)
            experiment_dictionary['wind vectors'].append(list(wind_vector))
            experiment_dictionary['perp_vectors'].append(list(return_list[2]))
            experiment_dictionary['par_vectors'].append(list(return_list[3]))
            experiment_dictionary['heading_unit_vectors'].append(list(return_list[4]))
            experiment_dictionary['trajectory_vectors'].append(list(trajectory_vector))

            #now asking if the fly, taking this trajectory, would get trapped
            list_of_tracking_probs = []
            list_of_distances_from_source = []
            print ()
            for trap_point in trap_point_list:
                distance_from_source = ask_if_fly_intersects_a_plume(trajectory = trajectory_vector,
                                                            trap_point = trap_point,
                                                            plume_angular_width = plume_angular_width,
                                                            wind_direction = wind_dir)
                if distance_from_source == inf:
                    print ('distance_from_source is infinite')
                    tracking_prob = 0
                else:
                    tracking_prob = calc_plume_track_prob(distance_from_source)
                list_of_tracking_probs.append(tracking_prob)
                list_of_distances_from_source.append(distance_from_source)
            print (list_of_tracking_probs)
            index, tracking_prob_value = max(enumerate(list_of_tracking_probs), key=operator.itemgetter(1))
            if tracking_prob_value > tracking_prob_cutoff:


                trap_point = trap_point_list[index]

                # print ('Fly heading:           ' +str(fly_heading*180/np.pi))
                print ('Fly trajectory:        ' + str(trajectory_angle*180/np.pi))
                print ('Wind direction:        '+ str(wind_dir*180/np.pi))
                trap_angle = np.arctan2(trap_point[1],trap_point[0])
                trap_angle = (trap_angle+2*np.pi)%(2*np.pi)
                print ('Trapped at trap angle: ' +str(trap_angle*180/np.pi))
                print ('Trapped at trap index: '  +str(index))
                print ('Tracking prob value:   ' +str(tracking_prob_value))
                print ('Distance from source:  ' +str(list_of_distances_from_source[index]))

        dictionary[str(wind_speed)][str(wind_dir)] = experiment_dictionary


with open('./sideslip_and_simple_plume_trap_counts.json', mode = 'w') as f:
    json.dump(dictionary,f, indent = 1)

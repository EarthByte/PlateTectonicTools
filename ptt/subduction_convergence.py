
"""
    Copyright (C) 2016 The University of Sydney, Australia
    
    This program is free software; you can redistribute it and/or modify it under
    the terms of the GNU General Public License, version 2, as published by
    the Free Software Foundation.
    
    This program is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.
    
    You should have received a copy of the GNU General Public License along
    with this program; if not, write to Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
"""


#############################################################
# Find the convergence rate of trenches (subduction zones). #
#############################################################


from __future__ import print_function
import math
import pygplates
import warnings


# Required pygplates version.
PYGPLATES_VERSION_REQUIRED = pygplates.Version(12)

# PyGPlates version 22 can handle topological lines (can get their sub-sub-segment plate IDs).
USING_PYGPLATES_VERSION_GREATER_EQUAL_22 = (hasattr(pygplates, 'Version') and pygplates.Version.get_imported_version() >= pygplates.Version(22))

# PyGPlates version 23 has a method to get overriding and subducting plates.
USING_PYGPLATES_VERSION_GREATER_EQUAL_23 = (hasattr(pygplates, 'Version') and pygplates.Version.get_imported_version() >= pygplates.Version(23))

# The default threshold sampling distance along trenches (subduction zones).
DEFAULT_THRESHOLD_SAMPLING_DISTANCE_DEGREES = 0.5
DEFAULT_THRESHOLD_SAMPLING_DISTANCE_KMS = math.radians(DEFAULT_THRESHOLD_SAMPLING_DISTANCE_DEGREES) * pygplates.Earth.equatorial_radius_in_kms

DEFAULT_TIME_RANGE_YOUNG_TIME = 0
DEFAULT_TIME_RANGE_OLD_TIME = 200
DEFAULT_TIME_INCREMENT = 1

DEFAULT_VELOCITY_DELTA_TIME = 1


# Determine the overriding and subducting plates of the subduction shared sub-segment.
#
# Note: This is now a method in PyGPlates version 23 called pygplates.ResolvedTopologicalSharedSubSegment.get_overriding_and_subducting_plates().
def find_overriding_and_subducting_plates(subduction_shared_sub_segment):
    
    # Get the subduction polarity of the nearest subducting line.
    subduction_polarity = subduction_shared_sub_segment.get_feature().get_enumeration(pygplates.PropertyName.gpml_subduction_polarity)
    if (not subduction_polarity) or (subduction_polarity == 'Unknown'):
        return

    # There should be two sharing topologies - one is the overriding plate and the other the subducting plate.
    sharing_resolved_topologies = subduction_shared_sub_segment.get_sharing_resolved_topologies()
    if len(sharing_resolved_topologies) != 2:
        return

    overriding_plate = None
    subducting_plate = None
    
    geometry_reversal_flags = subduction_shared_sub_segment.get_sharing_resolved_topology_geometry_reversal_flags()
    for index in range(2):

        sharing_resolved_topology = sharing_resolved_topologies[index]
        geometry_reversal_flag = geometry_reversal_flags[index]

        if sharing_resolved_topology.get_resolved_boundary().get_orientation() == pygplates.PolygonOnSphere.Orientation.clockwise:
            # The current topology sharing the subducting line has clockwise orientation (when viewed from above the Earth).
            # If the overriding plate is to the 'left' of the subducting line (when following its vertices in order) and
            # the subducting line is reversed when contributing to the topology then that topology is the overriding plate.
            # A similar test applies to the 'right' but with the subducting line not reversed in the topology.
            if ((subduction_polarity == 'Left' and geometry_reversal_flag) or
                (subduction_polarity == 'Right' and not geometry_reversal_flag)):
                overriding_plate = sharing_resolved_topology
            else:
                subducting_plate = sharing_resolved_topology
        else:
            # The current topology sharing the subducting line has counter-clockwise orientation (when viewed from above the Earth).
            # If the overriding plate is to the 'left' of the subducting line (when following its vertices in order) and
            # the subducting line is not reversed when contributing to the topology then that topology is the overriding plate.
            # A similar test applies to the 'right' but with the subducting line reversed in the topology.
            if ((subduction_polarity == 'Left' and not geometry_reversal_flag) or
                (subduction_polarity == 'Right' and geometry_reversal_flag)):
                overriding_plate = sharing_resolved_topology
            else:
                subducting_plate = sharing_resolved_topology
    
    if overriding_plate is None:
        return
    
    if subducting_plate is None:
        return
    
    return (overriding_plate, subducting_plate, subduction_polarity)


def subduction_convergence(
        rotation_features_or_model,
        topology_features,
        threshold_sampling_distance_radians,
        time,
        velocity_delta_time = 1.0,
        anchor_plate_id = 0):
    # Docstring in numpydoc format...
    """Find the convergence and absolute velocities sampled along trenches (subduction zones) at a particular geological time.
    
    Each sampled point along trench returns the following information:
    
    0 - longitude of sample point
    1 - latitude of sample point
    2 - subducting convergence (relative to trench) velocity magnitude (in cm/yr)
    3 - subducting convergence velocity obliquity angle (angle between trench normal vector and convergence velocity vector)
    4 - trench absolute (relative to anchor plate) velocity magnitude (in cm/yr)
    5 - trench absolute velocity obliquity angle (angle between trench normal vector and trench absolute velocity vector)
    6 - length of arc segment (in degrees) that current point is on
    7 - trench normal azimuth angle (clockwise starting at North, ie, 0 to 360 degrees) at current point
    8 - subducting plate ID
    9 - trench plate ID
    
    The obliquity angles are in the range (-180 180). The range (0, 180) goes clockwise (when viewed from above the Earth) from the
    trench normal direction to the velocity vector. The range (0, -180) goes counter-clockwise.
    You can change the range (-180, 180) to the range (0, 360) by adding 360 to negative angles.
    The trench normal is perpendicular to the trench and pointing toward the overriding plate.
    
    Note that the convergence velocity magnitude is negative if the plates are diverging (if convergence obliquity angle
    is greater than 90 or less than -90). And note that the absolute velocity magnitude is negative if the trench (subduction zone)
    is moving towards the overriding plate (if absolute obliquity angle is less than 90 or greater than -90) - note that this
    ignores the kinematics of the subducting plate.
    
    Parameters
    ----------
    rotation_features_or_model : pygplates.RotationModel, or any combination of str, pygplates.FeatureCollection, pygplates.Feature
        The rotation model can be specified as a RotationModel. Or it can be specified as a rotation feature collection,
        or rotation filename, or rotation feature, or sequence of rotation features, or a sequence (eg, list or tuple) of any combination
        of those four types.
    topology_features: any combination of str, pygplates.FeatureCollection, pygplates.Feature
        The topological boundary and network features and the topological section features they reference (regular and topological lines).
        Can be specified as a feature collection, or filename, or feature, or sequence of features, or a sequence (eg, list or tuple)
        of any combination of those four types.
    threshold_sampling_distance_radians: float
        Threshold sampling distance along trench (in radians).
    time: float
        The reconstruction time at which to query subduction convergence.
    velocity_delta_time: float, optional
        The delta time interval used for velocity calculations. Defaults to 1My.
    anchor_plate_id: int, optional
        The anchor plate of the rotation model. Defaults to zero.
    
    Returns
    -------
    list of tuples
        The results for all points sampled along trench.
        The size of the returned list is equal to the number of sampled points.
        Each tuple in the list corresponds to a point and has the following tuple items:
        
        * longitude of sample point
        * latitude of sample point
        * subducting convergence (relative to trench) velocity magnitude (in cm/yr)
        * subducting convergence velocity obliquity angle (angle between trench normal vector and convergence velocity vector)
        * trench absolute (relative to anchor plate) velocity magnitude (in cm/yr)
        * trench absolute velocity obliquity angle (angle between trench normal vector and trench absolute velocity vector)
        * length of arc segment (in degrees) that current point is on
        * trench normal azimuth angle (clockwise starting at North, ie, 0 to 360 degrees) at current point
        * subducting plate ID
        * trench plate ID
    
    Notes
    -----
    Each point in the output is the midpoint of a great circle arc between two adjacent points in the trench polyline.
    The trench normal vector used in the obliquity calculations is perpendicular to the great circle arc of each point (arc midpoint)
    and pointing towards the overriding plate (rather than away from it).
    
    Each trench is sampled at approximately uniform intervals along its length (specified via a threshold sampling distance).
    The sampling along the entire length of a trench is not exactly uniform. Each segment along a trench is sampled
    such that the samples have a uniform spacing that is less than or equal to the threshold sampling distance. However each segment
    in a trench might have a slightly different spacing distance (since segment lengths are not integer multiples of
    the threshold sampling distance).
    
    The trench normal (at each arc segment mid-point) always points *towards* the overriding plate.
    """
    
    # Check the imported pygplates version.
    if pygplates.Version.get_imported_version() < PYGPLATES_VERSION_REQUIRED:
        raise RuntimeError('Using pygplates version {0} but version {1} or greater is required'.format(
                pygplates.Version.get_imported_version(), PYGPLATES_VERSION_REQUIRED))
    
    # Turn rotation data into a RotationModel (if not already).
    rotation_model = pygplates.RotationModel(rotation_features_or_model)
    
    # Turn topology data into a list of features (if not already).
    topology_features = pygplates.FeaturesFunctionArgument(topology_features)
    
    # Resolve our topological plate polygons (and deforming networks) to the current 'time'.
    # We generate both the resolved topology boundaries and the boundary sections between them.
    resolved_topologies = []
    shared_boundary_sections = []
    pygplates.resolve_topologies(topology_features.get_features(), rotation_model, resolved_topologies, time, shared_boundary_sections, anchor_plate_id)
    
    # List of tesselated subduction zone (trench) shared subsegment points and associated convergence parameters
    # for the current 'time'.
    output_data = []
    
    # Iterate over the shared boundary sections of all resolved topologies.
    for shared_boundary_section in shared_boundary_sections:
    
        # Skip sections that are not subduction zones (trenches).
        if shared_boundary_section.get_feature().get_feature_type() != pygplates.FeatureType.gpml_subduction_zone:
            continue
        
        # Iterate over the shared sub-segments of the current subducting line.
        # These are the parts of the subducting line that actually contribute to topological boundaries.
        for shared_sub_segment in shared_boundary_section.get_shared_sub_segments():
        
            # Find the overriding and subducting plates on either side of the shared sub-segment.
            if USING_PYGPLATES_VERSION_GREATER_EQUAL_23:
                # PyGPlates version 23 has a method to get overriding and subducting plates.
                overriding_and_subducting_plates = shared_sub_segment.get_overriding_and_subducting_plates(return_subduction_polarity=True)
            else:
                # Otherwise call the function defined above.
                overriding_and_subducting_plates = find_overriding_and_subducting_plates(shared_sub_segment)
            if not overriding_and_subducting_plates:
                warnings.warn('Unable to find the overriding and subducting plates of the subducting sub-segment "{0}" at {1}Ma.\n'
                              '    Either the subduction polarity is not properly set or there are not exactly 2 topologies sharing the sub-segment.\n'
                              '    Ignoring current sub-segment.'.format(
                                  shared_sub_segment.get_feature().get_name(), time),
                              RuntimeWarning)
                continue
            overriding_plate, subducting_plate, subduction_polarity = overriding_and_subducting_plates
            overriding_plate_id = overriding_plate.get_feature().get_reconstruction_plate_id()
            subducting_plate_id = subducting_plate.get_feature().get_reconstruction_plate_id()
            
            # We need to reverse the trench vector direction if overriding plate is to
            # the right of the subducting line since great circle arc normal is always to the left.
            if subduction_polarity == 'Left':
                trench_normal_reversal = 1
            else:
                trench_normal_reversal = -1
            
            # The plate ID of the trench line (as opposed to the subducting plate).
            #
            # Update: The plate IDs of the trench line and overriding plate can differ
            # even in a non-deforming model due to smaller plates, not modelled by topologies, moving
            # differently than the larger topological plate being modelled - and the trench line
            # having plate IDs of the smaller plates near them. For that reason we use the plate ID
            # of the trench line whenever we can.
            #
            # If the current shared sub-segment is part of a topological line then we obtain its sub-sub-segments
            # (if we have pyGPlates version 22 or above). This is because trench lines that are
            # topological lines might actually be deforming (or intended to be deforming) and hence their
            # plate ID is not meaningful or at least we can't be sure whether it will be zero or the
            # overriding plate (or something else). In this case we look at the plate IDs of the
            # sub-sub-segments. However if we have pyGPlates version 21 or below then we cannot do this,
            # in which case (for a topological line) we'll use the overriding plate ID instead.
            #
            if USING_PYGPLATES_VERSION_GREATER_EQUAL_22:
                sub_segments_of_topological_line_sub_segment = shared_sub_segment.get_sub_segments()
                if sub_segments_of_topological_line_sub_segment:
                    # Iterate over the sub-sub-segments associated with the topological line.
                    for sub_sub_segment in sub_segments_of_topological_line_sub_segment:
                        trench_plate_id = sub_sub_segment.get_feature().get_reconstruction_plate_id()
                        sub_segment_geometry = sub_sub_segment.get_resolved_geometry()
                        _sub_segment_subduction_convergence(
                                output_data,
                                time,
                                sub_segment_geometry,
                                trench_plate_id,
                                subducting_plate_id,
                                trench_normal_reversal,
                                threshold_sampling_distance_radians,
                                velocity_delta_time,
                                rotation_model,
                                anchor_plate_id)
                else: # It's not a topological line...
                    trench_plate_id = shared_sub_segment.get_feature().get_reconstruction_plate_id()
                    sub_segment_geometry = shared_sub_segment.get_resolved_geometry()
                    _sub_segment_subduction_convergence(
                            output_data,
                            time,
                            sub_segment_geometry,
                            trench_plate_id,
                            subducting_plate_id,
                            trench_normal_reversal,
                            threshold_sampling_distance_radians,
                            velocity_delta_time,
                            rotation_model,
                            anchor_plate_id)
            else: # Cannot handle topological lines (so use overriding plate ID when one is detected)...
                if isinstance(shared_boundary_section.get_topological_section(), pygplates.ResolvedTopologicalLine):
                    trench_plate_id = overriding_plate_id
                else:
                    trench_plate_id = shared_sub_segment.get_feature().get_reconstruction_plate_id()
                sub_segment_geometry = shared_sub_segment.get_resolved_geometry()
                _sub_segment_subduction_convergence(
                        output_data,
                        time,
                        sub_segment_geometry,
                        trench_plate_id,
                        subducting_plate_id,
                        trench_normal_reversal,
                        threshold_sampling_distance_radians,
                        velocity_delta_time,
                        rotation_model,
                        anchor_plate_id)
    
    # Return data sorted since it's easier to compare results (when at least lon/lat is sorted).
    return sorted(output_data)


def _sub_segment_subduction_convergence(
        output_data,
        time,
        sub_segment_geometry,
        trench_plate_id,
        subducting_plate_id,
        trench_normal_reversal,
        threshold_sampling_distance_radians,
        velocity_delta_time,
        rotation_model,
        anchor_plate_id):
    
    # Get the rotation of the subducting plate relative to the trench line
    # from 'time + velocity_delta_time' to 'time'.
    convergence_relative_stage_rotation = rotation_model.get_rotation(
            time,
            subducting_plate_id,
            time + velocity_delta_time,
            trench_plate_id,
            anchor_plate_id=anchor_plate_id)
    #
    # In the following:
    #   * T is for Trench (subduction zone line)
    #   * S is subducting plate
    #   * A is anchor plate
    #
    # The trenches have been reconstructed using the rotation "R(0->t,A->T)":
    #
    #   reconstructed_geometry = R(0->t,A->T) * present_day_geometry
    #
    # We can write "R(0->t,A->T)" in terms of the convergence stage rotation "R(t+dt->t,T->S)" as:
    #
    #   R(0->t,A->T)  = R(0->t,A->S) * R(0->t,S->T)
    #                 = R(0->t,A->S) * inverse[R(0->t,T->S)]
    #                 = R(0->t,A->S) * inverse[R(t+dt->t,T->S) * R(0->t+dt,T->S)]
    #                 = R(0->t,A->S) * inverse[stage_rotation * R(0->t+dt,T->S)]
    #                 = R(0->t,A->S) * inverse[R(0->t+dt,T->S)] * inverse[stage_rotation]
    #                 = R(0->t,A->S) * R(0->t+dt,S->T) * inverse[stage_rotation]
    #
    # So to get the *reconstructed* subduction line geometry into the stage rotation reference frame
    # we need to rotate it by "inverse[R(0->t,A->S) * R(0->t+dt,S->T)]":
    #
    #   reconstructed_geometry = R(0->t,A->T) * present_day_geometry
    #                          = R(0->t,A->S) * R(0->t+dt,S->T) * inverse[stage_rotation] * present_day_geometry
    #   inverse[R(0->t,A->S) * R(0->t+dt,S->T)] * reconstructed_geometry = inverse[stage_rotation] * present_day_geometry
    #
    # Once we've done that we can calculate the velocities of those geometry points
    # using the stage rotation. Then the velocities need to be rotated back from the
    # stage rotation reference frame using the rotation "R(0->t,A->S) * R(0->t+dt,S->T)".
    # 
    from_convergence_stage_frame = (
        rotation_model.get_rotation(
                time,
                subducting_plate_id,
                anchor_plate_id=anchor_plate_id) *
        rotation_model.get_rotation(
                time + velocity_delta_time,
                trench_plate_id,
                fixed_plate_id=subducting_plate_id,
                anchor_plate_id=anchor_plate_id))
    to_convergence_stage_frame = from_convergence_stage_frame.get_inverse()
    
    # Get the rotation of the trench relative to the anchor plate
    # from 'time + velocity_delta_time' to 'time'.
    #
    # Note: We don't need to convert to and from the stage rotation reference frame
    # like the above convergence because...
    #
    #   R(0->t,A->T)  = R(t+dt->t,A->T) * R(0->t+dt,A->T)
    #
    #   reconstructed_geometry = R(0->t,A->T) * present_day_geometry
    #                          = R(t+dt->t,A->T) * R(0->t+dt,A->T) * present_day_geometry
    #
    # ...where *reconstructed* subduction line geometry is already in the frame of the stage rotation "R(0->t2,A->F)".
    trench_equivalent_stage_rotation = rotation_model.get_rotation(
            time,
            trench_plate_id,
            time + velocity_delta_time,
            anchor_plate_id=anchor_plate_id)
    
    # Ensure the shared sub-segment is tessellated to within the threshold sampling distance.
    tessellated_shared_sub_segment_polyline = (
            sub_segment_geometry.to_tessellated(threshold_sampling_distance_radians))
    
    # Iterate over the great circle arcs of the tessellated polyline to get the
    # arc midpoints, lengths and trench normals.
    # There is an arc between each adjacent pair of points in the polyline.
    arc_midpoints = []
    arc_lengths = []
    trench_normals = []
    for arc in tessellated_shared_sub_segment_polyline.get_segments():
        if not arc.is_zero_length():
            arc_midpoints.append(arc.get_arc_point(0.5))
            arc_lengths.append(arc.get_arc_length())
            # The normal to the trench in the direction of subduction (towards overriding plate).
            trench_normals.append(trench_normal_reversal * arc.get_great_circle_normal())
    
    # Shouldn't happen, but just in case the shared sub-segment polyline coincides with a point.
    if not arc_midpoints:
        return
    
    # The trench normals relative to North (azimuth).
    # Convert global 3D normal vectors to local (magnitude, azimuth, inclination) tuples (one tuple per point).
    trench_local_normals = pygplates.LocalCartesian.convert_from_geocentric_to_magnitude_azimuth_inclination(
            arc_midpoints, trench_normals)
    
    # Calculate the convergence velocities at the arc midpoints.
    #
    # Note; We need to convert the reconstructed geometry points into the convergence stage rotation
    # reference frame to calculate velocities and then convert the velocities using the
    # reverse transform as mentioned above.
    arc_midpoints_in_convergence_stage_frame = [
            to_convergence_stage_frame * arc_midpoint
                    for arc_midpoint in arc_midpoints]
    convergence_velocity_vectors_in_convergence_stage_frame = pygplates.calculate_velocities(
            arc_midpoints_in_convergence_stage_frame,
            convergence_relative_stage_rotation,
            velocity_delta_time,
            pygplates.VelocityUnits.cms_per_yr)
    convergence_velocity_vectors = [
            from_convergence_stage_frame * velocity
                    for velocity in convergence_velocity_vectors_in_convergence_stage_frame]
    
    # Calculate the trench absolute velocities at the arc midpoints.
    trench_absolute_velocity_vectors = pygplates.calculate_velocities(
            arc_midpoints, trench_equivalent_stage_rotation,
            velocity_delta_time, pygplates.VelocityUnits.cms_per_yr)
    
    for arc_index in range(len(arc_midpoints)):
        arc_midpoint = arc_midpoints[arc_index]
        arc_length = arc_lengths[arc_index]
        trench_normal = trench_normals[arc_index]
        trench_normal_azimuth = trench_local_normals[arc_index][1]
        lat, lon = arc_midpoint.to_lat_lon()
        
        # The direction towards which we rotate from the trench normal in a clockwise fashion.
        clockwise_direction = pygplates.Vector3D.cross(trench_normal, arc_midpoint.to_xyz())
        
        # Calculate the convergence rate parameters.
        convergence_velocity_vector = convergence_velocity_vectors[arc_index]
        if convergence_velocity_vector.is_zero_magnitude():
            convergence_velocity_magnitude = 0
            convergence_obliquity_degrees = 0
        else:
            convergence_velocity_magnitude = convergence_velocity_vector.get_magnitude()
            convergence_obliquity_degrees = math.degrees(pygplates.Vector3D.angle_between(
                    convergence_velocity_vector, trench_normal))
            # Anti-clockwise direction has range (0, -180) instead of (0, 180).
            if pygplates.Vector3D.dot(convergence_velocity_vector, clockwise_direction) < 0:
                convergence_obliquity_degrees = -convergence_obliquity_degrees
            
            # See if plates are diverging (moving away from each other).
            # If plates are diverging (moving away from each other) then make the
            # velocity magnitude negative to indicate this. This could be inferred from
            # the obliquity but it seems this is the standard way to output convergence rate.
            if math.fabs(convergence_obliquity_degrees) > 90:
                convergence_velocity_magnitude = -convergence_velocity_magnitude
        
        # Calculate the trench absolute rate parameters.
        trench_absolute_velocity_vector = trench_absolute_velocity_vectors[arc_index]
        if trench_absolute_velocity_vector.is_zero_magnitude():
            trench_absolute_velocity_magnitude = 0
            trench_absolute_obliquity_degrees = 0
        else:
            trench_absolute_velocity_magnitude = trench_absolute_velocity_vector.get_magnitude()
            trench_absolute_obliquity_degrees = math.degrees(pygplates.Vector3D.angle_between(
                    trench_absolute_velocity_vector, trench_normal))
            # Anti-clockwise direction has range (0, -180) instead of (0, 180).
            if pygplates.Vector3D.dot(trench_absolute_velocity_vector, clockwise_direction) < 0:
                trench_absolute_obliquity_degrees = -trench_absolute_obliquity_degrees
            
            # See if the trench absolute motion is heading in the direction of the
            # overriding plate. If it is then make the velocity magnitude negative to
            # indicate this. This could be inferred from the obliquity but it seems this
            # is the standard way to output trench velocity magnitude.
            #
            # Note that we are not calculating the motion of the trench
            # relative to the overriding plate - they are usually attached to each other
            # and hence wouldn't move relative to each other.
            if math.fabs(trench_absolute_obliquity_degrees) < 90:
                trench_absolute_velocity_magnitude = -trench_absolute_velocity_magnitude
        
        # The data will be output in GMT format (ie, lon first, then lat, etc).
        output_data.append((
                lon,
                lat,
                convergence_velocity_magnitude,
                convergence_obliquity_degrees,
                trench_absolute_velocity_magnitude,
                trench_absolute_obliquity_degrees,
                math.degrees(arc_length),
                math.degrees(trench_normal_azimuth),
                subducting_plate_id,
                trench_plate_id))


def write_output_file(output_filename, output_data):
    with open(output_filename, 'w') as output_file:
        for output_line in output_data:
            output_file.write(' '.join(str(item) for item in output_line) + '\n')


def create_coverage_feature_from_convergence_data(
        subduction_convergence_data,
        time):
    """Create a feature with a coverage geometry containing the calculated convergence and absolute velocity data.
    
    Parameters
    ----------
    subduction_convergence_data : list of tuples
        The subduction convergence data calculated by :func:`subduction_convergence`.
        Each tuple in the list contains the calculated data for a single sample point on a trench line.
    time: float
        The reconstruction time associated with the subduction convergence data.
    
    Returns
    -------
    pygplates.Feature
        The feature with a coverage geometry containing the calculated convergence and absolute velocity data.
    """
    
    # Convert the list of tuples (one tuple per sample point) into a tuple of lists (one list per data parameter).
    (all_lon,
     all_lat,
     all_convergence_velocity_magnitude_cm_per_yr,
     all_convergence_obliquity_degrees,
     all_trench_absolute_velocity_magnitude_cm_per_yr,
     all_trench_absolute_obliquity_degrees,
     all_subducting_length_degrees,
     all_trench_normal_azimuth_degrees,
     all_subducting_plate_id,
     all_trench_plate_id) = zip(*subduction_convergence_data)
    
    # Put all convergence data for the current reconstruction time into a single feature.
    coverage_feature = pygplates.Feature()
    
    # Make it only appear at 'time'.
    coverage_feature.set_valid_time(time + 0.5, time - 0.5)
    
    # Add each data parameter as a separate scalar coverage.
    coverage_geometry = pygplates.MultiPointOnSphere(zip(all_lat, all_lon))
    coverage_scalars = {
        pygplates.ScalarType.create_gpml('ConvergenceVelocityMagnitude') : all_convergence_velocity_magnitude_cm_per_yr,
        pygplates.ScalarType.create_gpml('ConvergenceObliquityDegrees') : all_convergence_obliquity_degrees,
        pygplates.ScalarType.create_gpml('TrenchAbsoluteVelocityMagnitude') : all_trench_absolute_velocity_magnitude_cm_per_yr,
        pygplates.ScalarType.create_gpml('TrenchAbsoluteObliquityDegrees') : all_trench_absolute_obliquity_degrees,
        pygplates.ScalarType.create_gpml('SubductingLengthDegrees') : all_subducting_length_degrees,
        pygplates.ScalarType.create_gpml('TrenchNormalAzimuthDegrees') : all_trench_normal_azimuth_degrees,
        pygplates.ScalarType.create_gpml('SubductingPlateId') : all_subducting_plate_id,
        pygplates.ScalarType.create_gpml('TrenchPlateId') : all_trench_plate_id,
    }
    coverage_feature.set_geometry((coverage_geometry, coverage_scalars))
    
    return coverage_feature


def subduction_convergence_over_time(
        output_filename_prefix,
        output_filename_extension,
        rotation_filenames,
        topology_filenames,
        threshold_sampling_distance_radians,
        time_young,
        time_old,
        time_increment,
        velocity_delta_time = 1.0,
        anchor_plate_id = 0,
        output_gpml_filename = None):
    
    # Check the imported pygplates version.
    if pygplates.Version.get_imported_version() < PYGPLATES_VERSION_REQUIRED:
        raise RuntimeError('Using pygplates version {0} but version {1} or greater is required'.format(
                pygplates.Version.get_imported_version(), PYGPLATES_VERSION_REQUIRED))
    
    if time_increment <= 0:
        raise ValueError('The time increment "{0}" is not positive and non-zero.'.format(time_increment))
    
    if time_young > time_old:
        raise ValueError('The young time {0} is older (larger) than the old time {1}.'.format(time_young, time_old))
    
    rotation_model = pygplates.RotationModel(rotation_filenames)
    
    # Read/parse the topological features once so we're not doing at each time iteration.
    topology_features = [pygplates.FeatureCollection(topology_filename)
            for topology_filename in topology_filenames]
    
    if output_gpml_filename:
        coverage_features = []
    
    # Iterate over the time rage.
    time = time_young
    while time <= pygplates.GeoTimeInstant(time_old):
        
        # print('Time {0}'.format(time))
        
        # Returns a list of tesselated trench points and associated convergence parameters
        # to write to the output file for the current 'time'.
        output_data = subduction_convergence(
                rotation_model,
                topology_features,
                threshold_sampling_distance_radians,
                time,
                velocity_delta_time,
                anchor_plate_id)
        
        if output_data:
            output_filename = '{0}_{1:0.2f}.{2}'.format(output_filename_prefix, time, output_filename_extension)
            write_output_file(output_filename, output_data)
            
            # Also keep track of convergence data if we need to write out a GPML file.
            if output_gpml_filename:
                coverage_feature = create_coverage_feature_from_convergence_data(output_data, time)
                coverage_features.append(coverage_feature)

        # Increment the time further into the past.
        time += time_increment
    
    if output_gpml_filename:
        # Write out all coverage features to a single GPML file.
        pygplates.FeatureCollection(coverage_features).write(output_gpml_filename)
    
    return 0 # Success


if __name__ == '__main__':
    
    import argparse
    import os.path
    import sys
    import traceback
    
    ########################
    # Command-line parsing #
    ########################
    
    def warning_format(message, category, filename, lineno, file=None, line=None):
        # return '{0}:{1}: {1}:{1}\n'.format(filename, lineno, category.__name__, message)
        return '{0}: {1}\n'.format(category.__name__, message)
    
    # Print the warnings without the filename and line number.
    # Users are not going to want to see that.
    warnings.formatwarning = warning_format
   
    def main():
        
        __description__ = \
    """Find the convergence rates along trenches (subduction zones) over time.
    
    For each time (over a range of times) an output xy file is generated containing the resolved trenches
    (with point locations as the first two columns x and y) and the following convergence rate parameters in subsequent columns:
    
      - subducting convergence (relative to trench) velocity magnitude (in cm/yr)
      - subducting convergence velocity obliquity angle (angle between trench normal vector and convergence velocity vector)
      - trench absolute (relative to anchor plate) velocity magnitude (in cm/yr)
      - trench absolute velocity obliquity angle (angle between trench normal vector and trench absolute velocity vector)
      - length of arc segment (in degrees) that current point is on
      - trench normal azimuth angle (clockwise starting at North, ie, 0 to 360 degrees) at current point
      - subducting plate ID
      - trench plate ID
    
    The obliquity angles are in the range (-180 180). The range (0, 180) goes clockwise (when viewed from above the Earth) from the
    trench normal direction to the velocity vector. The range (0, -180) goes counter-clockwise.
    You can change the range (-180, 180) to the range (0, 360) by adding 360 to negative angles.
    
    Note that the convergence velocity magnitude is negative if the plates are diverging (if convergence obliquity angle
    is greater than 90 or less than -90). And note that the absolute velocity magnitude is negative if the trench
    is moving towards the overriding plate (if absolute obliquity angle is less than 90 or greater than -90) - note that this
    ignores the kinematics of the subducting plate.
    
    Each point in the output is the midpoint of a great circle arc between two adjacent points in the trench polyline.
    The trench normal vector used in the obliquity calculations is perpendicular to the great circle arc of each point (arc midpoint)
    and pointing towards the overriding plate (rather than away from it).
    
    Each trench is sampled at approximately uniform intervals along its length (specified via a threshold sampling distance).
    The sampling along the entire length of a trench is not exactly uniform. Each segment along a trench is sampled
    such that the samples have a uniform spacing that is less than or equal to the threshold sampling distance. However each segment
    in a trench might have a slightly different spacing distance (since segment lengths are not integer multiples of
    the threshold sampling distance).
    
    The trench normal (at each arc segment mid-point) always points *towards* the overriding plate.

    NOTE: Separate the positional and optional arguments with '--' (workaround for bug in argparse module).
    For example...

    python %(prog)s -r rotations.rot -m topologies.gpml -t 0 200 -i 1 -v 1 -d 0.5 -e xy -- convergence
     """
        
        # The command-line parser.
        parser = argparse.ArgumentParser(description = __description__, formatter_class=argparse.RawDescriptionHelpFormatter)
        
        parser.add_argument('-r', '--rotation_filenames', type=str, nargs='+', required=True,
                metavar='rotation_filename', help='One or more rotation files.')
        parser.add_argument('-m', '--topology_filenames', type=str, nargs='+', required=True,
                metavar='topology_filename', help='One or more topology files to generate resolved subducting lines.')
        parser.add_argument('-a', '--anchor', type=int, default=0,
                dest='anchor_plate_id',
                help='Anchor plate id used for reconstructing. Defaults to zero.')
        
        # Can specify only one of '-i', '-l' or '-t'.
        threshold_sampling_distance_group = parser.add_mutually_exclusive_group()
        threshold_sampling_distance_group.add_argument('-d', '--threshold_sampling_distance_degrees', type=float,
                help='Threshold sampling distance along trenches (in degrees). '
                    'Defaults to {0} degrees.'.format(DEFAULT_THRESHOLD_SAMPLING_DISTANCE_DEGREES))
        threshold_sampling_distance_group.add_argument('-k', '--threshold_sampling_distance_kms', type=float,
                help='Threshold sampling distance along trenches (in Kms). '
                    'Defaults to {0:.2f} Kms (which is equivalent to {1} degrees).'.format(
                            DEFAULT_THRESHOLD_SAMPLING_DISTANCE_KMS,
                            DEFAULT_THRESHOLD_SAMPLING_DISTANCE_DEGREES))

        parser.add_argument('-t', '--time_range', type=float, nargs=2,
                metavar=('young_time', 'old_time'),
                default=[DEFAULT_TIME_RANGE_YOUNG_TIME, DEFAULT_TIME_RANGE_OLD_TIME],
                help='The time range (in Ma) from young time to old time. '
                    'Defaults to {0} -> {1} Ma.'.format(
                        DEFAULT_TIME_RANGE_YOUNG_TIME, DEFAULT_TIME_RANGE_OLD_TIME))
        
        def parse_positive_number(value_string):
            try:
                value = float(value_string)
            except ValueError:
                raise argparse.ArgumentTypeError("%s is not a number" % value_string)
            
            if value <= 0:
                raise argparse.ArgumentTypeError("%g is not a positive number" % value)
            
            return value
        
        parser.add_argument('-i', '--time_increment', type=parse_positive_number,
                default=DEFAULT_TIME_INCREMENT,
                help='The time increment in My. Defaults to {0} My.'.format(DEFAULT_TIME_INCREMENT))
        
        parser.add_argument('-v', '--velocity_delta_time', type=parse_positive_number,
                default=DEFAULT_VELOCITY_DELTA_TIME,
                help='The delta time interval used to calculate velocities in My. '
                    'Defaults to {0} My.'.format(DEFAULT_VELOCITY_DELTA_TIME))
        
        parser.add_argument('-g', '--output_gpml_filename', type=str,
                help='Optional GPML output filename to contain the subduction convergence data for all specified times. '
                     'This can then be loaded into GPlates to display the data as scalar coverages.')
        
        parser.add_argument('output_filename_prefix', type=str,
                help='The output filename prefix. An output file is created for each geological time in the sequence where '
                    'the filename suffix contains the time and the filename extension.')
        parser.add_argument('-e', '--output_filename_extension', type=str, default='xy',
                help='The output xy filename extension. Defaults to "xy".')
        parser.add_argument('-w', '--ignore_topology_warnings', action="store_true",
                help='If specified then topology warnings are ignored (not output). '
                     'These are the warnings about not finding the overriding and subducting plates.')
        
        
        # Parse command-line options.
        args = parser.parse_args()
        
        # Topology warnings correspond to calls to "warnings.warn(... , RuntimeWarning)".
        # Ignore them if requested.
        if args.ignore_topology_warnings:
            warnings.simplefilter('ignore', RuntimeWarning)
        
        if args.time_range[0] > args.time_range[1]:
            raise argparse.ArgumentTypeError("First (young) value in time range is greater than second (old) value")
        
        # Determine threshold sampling distance.
        if args.threshold_sampling_distance_degrees:
            threshold_sampling_distance_radians = math.radians(args.threshold_sampling_distance_degrees)
        elif args.threshold_sampling_distance_kms:
            threshold_sampling_distance_radians = args.threshold_sampling_distance_kms / pygplates.Earth.equatorial_radius_in_kms
        else: # default...
            threshold_sampling_distance_radians = math.radians(DEFAULT_THRESHOLD_SAMPLING_DISTANCE_DEGREES)
        
        return_code = subduction_convergence_over_time(
                args.output_filename_prefix,
                args.output_filename_extension,
                args.rotation_filenames,
                args.topology_filenames,
                threshold_sampling_distance_radians,
                args.time_range[0],
                args.time_range[1],
                args.time_increment,
                args.velocity_delta_time,
                args.anchor_plate_id,
                args.output_gpml_filename)
        if return_code is None:
            sys.exit(1)
            
        sys.exit(0)
    
    try:
        main()
        sys.exit(0)
    except Exception as exc:
        print('ERROR: {0}: {1}'.format(os.path.basename(__file__), exc), file=sys.stderr)
        # Uncomment this to print traceback to location of raised exception.
        # traceback.print_exc()
        sys.exit(1)

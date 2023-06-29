"""
    Copyright (C) 2023 The University of Sydney, Australia
    
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


from  collections import deque
from inspect import signature
import math
import os
import os.path
from .utils import points_in_polygons
from .utils import points_spatial_tree
from .utils import proximity_query
import pygplates
import sys
import numpy as np


#############################################################################################################
# Calculate continent contours, and fragmentation index (global perimeter-to-area ratio), at various times.
#
# TODO: Replace the internal uniform lat/lon sampling with a sampling that’s uniform across the sphere
#       (so that the contour resolution doesn’t favour contours near the North/South poles).
#############################################################################################################
   

class ContouredContinent(object):
    """
    Class to represent the contour around overlapping/abutting continental blocks.
    """
    
    def __init__(self):
        self._polygons_including_continent = []
        self._polygons_excluding_continent = []
    

    def _add_landmass(self, landmass_polygon, landmass_polygon_interior_contains_continent, contour_polylines):
        """
        Add a landmass to this contoured continent.

        This is a polygon representing the landmass region/area, and one or more polylines representing the landmass boundary(s).
        
        If 'landmass_polygon_interior_contains_continent' is *true* then the interior of the polygon represents continental crust.
        In this case the polygon can have interior rings (holes) which represent oceanic crust.
        And it's also possible for a landmass polygon to be inside the interior hole of another landmass polygon, and so on.
        For example, an continental island inside an ocean basin that, in turn, is inside a larger continent.
        
        If 'landmass_polygon_interior_contains_continent' is *false* then the landmass is actually a continent that covers
        the entire globe except for a few oceanic holes. And the specified 'polygon' is actually one of those oceanic holes.
        This is a special case because you can't have a single global polygon with only interior holes (and no exterior ring).
        So intead we allow for multiple exterior ring polygons to represent these oceanic holes (and treat them specially).
        """
        if landmass_polygon_interior_contains_continent:
            self._polygons_including_continent.append(landmass_polygon)
        else:
            self._polygons_excluding_continent.append(landmass_polygon)
        
        self._polylines = contour_polylines


    def get_contours(self):
        """The polyline contours representing the boundaries of continental crust."""
        return self._polylines
    

    def are_points_inside(self, points, points_spatial_tree=None):
        """Returns a numpy 1D boolean array with same length as 'points' (and in same order) containing True for each point inside this contoured continent."""

        # A special (unlikely) case is a single continent covering the entire globe (area 4*pi).
        # This happens when there are no inclusive polygons and no exclusive polygons.
        if not self._polygons_including_continent and not self._polygons_excluding_continent:
            # All points *include* continent.
            return np.full(len(points), True)

        # Improve efficiency by re-using spatial tree of points if caller provides it (otherwise create our own).
        if not points_spatial_tree:
            points_spatial_tree = points_spatial_tree.PointsSpatialTree(points)
        
        # If we have any polygons that *exclude* continent then it means the continent landmass covers the entire globe
        # except these excluding polygons.
        if self._polygons_excluding_continent:
            # By default all points are considered *inside* this contoured continent unless proven *outside*.
            points_inside = np.full(len(points), True)

            # See if the points are inside any of the exclusive polygons.
            exclusive_polygons_containing_points = points_in_polygons.find_polygons_using_points_spatial_tree(
                    points,
                    points_spatial_tree,
                    self._polygons_excluding_continent)
            
            # Any points *inside* an exclusive polygon are considered to be *outside* this contoured continent.
            #
            # Note: If there are any inclusive polygons (which must be inside these exclusive polygons) then later
            #       they will include some of these points that we just excluded.
            for point_index in range(len(points)):
                if exclusive_polygons_containing_points[point_index]:
                    points_inside[point_index] = False
            
        # else all polygons *include* continent (ie, none exclude continent)...
        else:
            # By default all points are considered *outside* this contoured continent unless proven *inside*.
            points_inside = np.full(len(points), False)

        if self._polygons_including_continent:
            # See if the points are inside any of the inclusive polygons.
            inclusive_polygons_containing_points = points_in_polygons.find_polygons_using_points_spatial_tree(
                    points,
                    points_spatial_tree,
                    self._polygons_including_continent)
            
            # Any points *inside* an inclusive polygon are considered to be *inside* this contoured continent.
            for point_index in range(len(points)):
                if inclusive_polygons_containing_points[point_index]:
                    points_inside[point_index] = True

        return points_inside


    def get_perimeter(self):
        """Sum of the length of the contour boundaries of this contoured continent (in radians)."""
        return math.fsum(polyline.get_arc_length() for polyline in self._polylines)


    def get_area(self):
        """The area of this contoured continent (in square radians, as known as steradians)."""

        # A special (unlikely) case is a single continent covering the entire globe (area 4*pi).
        # This happens when there are no inclusive polygons and no exclusive polygons.
        if not self._polygons_including_continent and not self._polygons_excluding_continent:
            return 4 * math.pi
        
        area = 0.0
        
        # Note that we can get one or more polygons that *exclude* continent.
        # This can happen when contouring a large landmass such that there is no contour that is an exterior ring.
        # In this case the landmass covers the entire globe except for a few oceanic holes.
        # However you can't have a single polygon with only interior holes (and no exterior ring).
        # So intead we allow for multiple exterior ring polygons to represent these holes (and treat this as a special case).
        if self._polygons_excluding_continent:
            # And since we can't have an exterior ring covering the entire globe we need to explicitly add the area of the entire globe.
            area += 4 * math.pi

            # Subtract the area of the exclusive holes.
            for polygon in self._polygons_excluding_continent:
                    area -= polygon.get_area()

        # Add the areas of polygons that include continent.
        for polygon in self._polygons_including_continent:
                area += polygon.get_area()

        return area


    def get_perimeter_area_ratio(self):
        """The perimeter divided by the area (in units of 1/radians)."""
        return self.get_perimeter() / self.get_area()


class ContinentContouring(object):
    """
    Class to calculate continent mask, contours, and fragmentation index (global perimeter-to-area ratio), at various times.
    """
    
    def __init__(
            self,
            rotaton_model_or_features,
            continent_features,  # regular features (not topologies)

            # The grid spacing (in degrees) between points in the grid used for contouring/aggregrating blocks of continental polygons.
            continent_contouring_point_spacing_degrees,

            # Optional parameter specifying area threshold (in square radians) for contoured continents.
            #
            # Contoured continents with area smaller than this threshold will be excluded.
            # If this parameter is not specified then no area threshold is applied.
            #
            # Can also be a function (accepting time in Ma) and returning the area threshold.
            #
            # Note: Units here are for normalised sphere (ie, steradians or square radians) so full Earth area is 4*pi.
            #       So 0.1 covers an area of approximately 4,000,000 km^2 (ie, 0.1 * 6371^2, where Earth radius is 6371km).
            #       Conversely 4,000,000 km^2 is equivalent to (4,000,000 / 6371^2) steradians.
            continent_contouring_area_threshold_steradians = None,

            # Optional parameter specifying a distance (in radians) to expand contours ocean-ward - this also
            # ensures small gaps between continents are ignored during contouring.
            #
            # The continent(s) will be expanded by a buffer of this distance (in radians) when contouring/aggregrating blocks of continental polygons.
            # If this parameter is not specified then buffer expansion is not applied.
            #
            # This parameter can also be a function (that returns the distance).
            # The function can have a single function argument: (1) accepting time (in Ma).
            # Or it can have two function arguments: (1) the first accepting time (in Ma) and
            # (2) the second accepting the contoured continent (a 'ContouredContinent' object)
            # of the (unexpanded) contoured continent that the buffer/gap distance will apply to.
            # Hence a function with *two* arguments means a different buffer/gap distance can be specified for each contoured continent (eg, based on its area).
            #
            # Note: Units here are for normalised sphere (ie, radians).
            #       So 1.0 radian is approximately 6371 km (where Earth radius is 6371 km).
            #       Also 1.0 degree is approximately 110 km.
            continent_contouring_buffer_and_gap_distance_radians = None,

            # Optional parameter specifying area threshold (in square radians) when creating continent contours.
            #
            # Polygon contours that exclude continental crust and have an area smaller than this threshold will be excluded
            # (meaning they will now *include* continental crust, thus removing the contour).
            # This is useful for removing small holes inside continents.
            # If this parameter is not specified then no area threshold is applied.
            #
            # Can also be a function (accepting time in Ma) and returning the area threshold.
            #
            # Note: Units here are for normalised sphere (ie, steradians or square radians) so full Earth area is 4*pi.
            #       So 0.1 covers an area of approximately 4,000,000 km^2 (ie, 0.1 * 6371^2, where Earth radius is 6371km).
            #       Conversely 4,000,000 km^2 is equivalent to (4,000,000 / 6371^2) steradians.
            continent_exclusion_area_threshold_steradians = None):
        
        # Make sure pygplates has support for interior rings in polygons.
        if pygplates.Version.get_imported_version() < pygplates.Version(0, 36):
            raise RuntimeError(
                "Using pygplates version {} but version 0.36 or greater is required to support interior rings in polygons".format(
                    pygplates.Version.get_imported_version()))
        
        self.rotation_model = pygplates.RotationModel(rotaton_model_or_features)
        self.continent_features = continent_features

        if continent_contouring_area_threshold_steradians:
            # Convert area threshold to a function of time, if not already a function.
            if callable(continent_contouring_area_threshold_steradians):
                # We can call the specified function directly.
                self.continent_contouring_area_threshold_steradians_function = continent_contouring_area_threshold_steradians
            else:
                # Use a delegate function that returns the specified parameter.
                def continent_contouring_area_threshold_steradians_function(age):
                    return continent_contouring_area_threshold_steradians
                self.continent_contouring_area_threshold_steradians_function = continent_contouring_area_threshold_steradians_function
        else:  # no area threshold (specified either None or zero)
            self.continent_contouring_area_threshold_steradians_function = None
        
        if continent_contouring_buffer_and_gap_distance_radians:
            # Convert buffer/gap threshold to a function of time, if not already a function.
            if callable(continent_contouring_buffer_and_gap_distance_radians):
                callable_signature = signature(continent_contouring_buffer_and_gap_distance_radians)
                callable_num_args = len(callable_signature.parameters)
                if not (callable_num_args == 1 or callable_num_args == 2):
                    raise TypeError('Buffer/gap distance is a callable but does not have 1 or 2 arguments')
                if callable_num_args == 2:
                    # We can call the specified function directly.
                    self.continent_contouring_buffer_and_gap_distance_radians_function = continent_contouring_buffer_and_gap_distance_radians
                else:  # callable_num_args == 1
                    # The specified function only accepts age (not area).
                    # So use a delegate function that calls it and ignores area.
                    def continent_contouring_buffer_and_gap_distance_radians_function(age, area):
                        return continent_contouring_buffer_and_gap_distance_radians(age)
                    self.continent_contouring_buffer_and_gap_distance_radians_function = continent_contouring_buffer_and_gap_distance_radians_function
            else:
                # Use a delegate function that returns the specified parameter.
                def continent_contouring_buffer_and_gap_distance_radians_function(age, area):
                    return continent_contouring_buffer_and_gap_distance_radians
                self.continent_contouring_buffer_and_gap_distance_radians_function = continent_contouring_buffer_and_gap_distance_radians_function
        else:  # no buffer/gap distance (specified either None or zero)
            self.continent_contouring_buffer_and_gap_distance_radians_function = None

        if continent_exclusion_area_threshold_steradians:
            # Convert area threshold to a function of time, if not already a function.
            if callable(continent_exclusion_area_threshold_steradians):
                # We can call the specified function directly.
                self.continent_exclusion_area_threshold_steradians_function = continent_exclusion_area_threshold_steradians
            else:
                # Use a delegate function that returns the specified parameter.
                def continent_exclusion_area_threshold_steradians_function(age):
                    return continent_exclusion_area_threshold_steradians
                self.continent_exclusion_area_threshold_steradians_function = continent_exclusion_area_threshold_steradians_function
        else:  # no area threshold (specified either None or zero)
            self.continent_exclusion_area_threshold_steradians_function = None

        # The number of latitudes (including -90 and 90).
        self.contouring_grid_num_latitudes = int(math.ceil(180.0 / continent_contouring_point_spacing_degrees)) + 1
        # The number of longitudes (including -180 and 180).
        self.contouring_grid_num_longitudes = 2 * (self.contouring_grid_num_latitudes - 1) + 1

        self.contouring_point_spacing_degrees = 180.0 / (self.contouring_grid_num_latitudes - 1)

        # A point grid to calculate contour polygons representing the boundary of reconstructed static polygons that overlap each other.
        #
        # NOTE: We must generate points on the dateline (ie, at both longitude -180 and 180) since the
        #       contouring alorithm depends on this. We also generate points at the North and South poles
        #       for the same reason.
        lats = np.linspace(-90.0, 90.0, self.contouring_grid_num_latitudes)
        lons = np.linspace(-180.0, 180.0, self.contouring_grid_num_longitudes)

        # Create a multipoint grid of points ordered by longitude first then latitude.
        contouring_longitude_array, contouring_latitude_array = np.meshgrid(lons, lats)
        self.contouring_points = pygplates.MultiPointOnSphere(
                zip(contouring_latitude_array.flatten(), contouring_longitude_array.flatten()))

        # Improve efficiency by re-using spatial tree of contouring points over time (when finding points in polygons and finding points near polygons).
        #
        # First calculate the subdivision depth to avoid doing too many point-in-polygon tests (for example) for each spatial tree node.
        # The lat/lon width of a root quad tree node in the spatial tree is 90 degrees which is 'n/2' points wide (where 'n' is 'self.contouring_grid_num_latitudes').
        # So a leaf node at subdivision depth 'D' is '(n/2) / 2^D' = 'n / 2^(D+1)'. The number of points is the square of that 'N = '(n / 2^(D+1)) ^ 2'.
        # Rearranging that gives the subdivision depth 'D' in terms of the number of points we'd like in a deepest (leaf) node N:
        #   D = log2(n / sqrt(N)) - 1
        max_num_points_per_spatial_tree_node = 32  # N
        points_spatial_tree_subdivision_depth = math.ceil(math.log(self.contouring_grid_num_latitudes / math.sqrt(max_num_points_per_spatial_tree_node), 2) - 1)  # D
        self.contouring_points_spatial_tree = points_spatial_tree.PointsSpatialTree(self.contouring_points, points_spatial_tree_subdivision_depth)
    
    
    def get_fragmentation(
            self,
            age):
        """
        Calculate the continental fragmentation index (global perimeter-to-area ratio) at the specified time.
        """

        # Get the contoured continents representing the boundary(s) of the reconstructed continent polygons that overlap each other.
        contoured_continents = self.get_contoured_continents(age)

        total_perimeter = 0.0
        total_area = 0.0
        
        # Update total perimeter and area.
        for contoured_continent in contoured_continents:
            total_perimeter += contoured_continent.get_perimeter()
            total_area += contoured_continent.get_area()

        # Avoid divide-by-zero.
        if total_area == 0.0:
            return 0.0
        
        #print(' global perimeter/area:', total_perimeter / total_area / pygplates.Earth.equatorial_radius_in_kms, 'km-1'); sys.stdout.flush()
        
        #print('age:', age, 'frag_index (1/km):', total_perimeter / total_area / 6371.0); sys.stdout.flush()
        return total_perimeter / total_area


    def get_continent_mask(
            self,
            age):
        """
        Reconstruct the continents (specified in constructor) and then find the latitude/longitude grid points that are inside them.

        The grid spacing of these grid points was specified in the constructor.

        Returns a 2D boolean numpy array of shape (num_latitudes, num_longitudes).
        Note that when writing to a NetCDF grid file you can convert to floating-point (with "continent_mask.astype('float')").
        """

        contoured_continents = self.get_contoured_continents(age)
        
        return self.calculate_continent_mask(contoured_continents)
    
    
    def get_contoured_continents(
            self,
            age):
        """
        Reconstruct the continents (specified in constructor) and then find their boundaries as contoured continents.

        Returns a list of 'ContouredContinent'.
        """

        reconstructed_continent_polygons = self.get_reconstructed_continent_polygons(age)
        
        return self.calculate_contoured_continents(reconstructed_continent_polygons, age)
    
    
    def get_continent_mask_and_contoured_continents(
            self,
            age):
        """
        Reconstruct the continents (specified in constructor) and then find both their boundaries as contoured continents and
        the latitude/longitude grid points that are inside them.

        Returns a 2-tuple of (a 2D boolean numpy array of shape (num_latitudes, num_longitudes), a list of 'ContouredContinent').
        """

        contoured_continents = self.get_contoured_continents(age)
        
        continent_mask = self.calculate_continent_mask(contoured_continents)
    
        return continent_mask, contoured_continents
    
    
    def get_reconstructed_continent_polygons(
            self,
            age):
        """
        Reconstruct the continents (specified in constructor).

        Note that these are just the original continent polygons (but reconstructed).
        They are NOT contoured, so they can still overlap/abutt each other.

        Returns a list of 'pygplates.PolygonOnSphere'.
        """
        
        # Reconstruct static continental polygons.
        reconstructed_feature_geometries = []
        pygplates.reconstruct(self.continent_features, self.rotation_model, reconstructed_feature_geometries, age)
        
        # Get a list of polygons.
        #
        # We should have polygons (not polylines) but turn into a polygon if happens to be a polyline
        # (but that actually only works if the polyline is a closed loop and not just part of a polygon's boundary).
        return [pygplates.PolygonOnSphere(reconstructed_feature_geometry.get_reconstructed_geometry())
                for reconstructed_feature_geometry in reconstructed_feature_geometries]


    def calculate_continent_mask(
            self,
            contoured_continents):
        """
        Find the latitude/longitude grid points that are inside the specified contoured continents.

        The grid spacing of these grid points was specified in the constructor.

        Returns a 2D boolean numpy array of shape (num_latitudes, num_longitudes).
        Note that when writing to a NetCDF grid file you can convert to floating-point (with "continent_mask.astype('float')").
        """
    
        return self._find_grid_points_inside_contoured_continents(contoured_continents)

    
    def calculate_contoured_continents(
            self,
            continent_polygons,
            age = 0):
        """
        Find the boundaries of the specified (potentially overlapping/abutting) continent polygons as contoured continents.

        Note that small contoured continent islands with area less than the area threshold will NOT get returned.

        The 'age' is only used to look up the time-dependent thresholds (passed into constructor).
        If threshold does not vary with time then 'age' does not need to be specified (defaults to present day).

        Returns a list of 'ContouredContinent'.
        """
        
        # Grid points inside the continent polygons.
        grid_points_inside_continents = self._find_grid_points_inside_polygons(continent_polygons)

        # Contour the grid points that are inside the continent polygons.
        contoured_continents = self._find_contoured_continents(grid_points_inside_continents)

        # If a contoured continent area threshold was specified.
        if self.continent_contouring_area_threshold_steradians_function:
            # If the area threshold is non-zero then exclude those contoured continents with area below the threshold.
            # And also exclude any grid points inside those rejected contoured continents.
            continent_contouring_area_threshold_steradians = self.continent_contouring_area_threshold_steradians_function(age)
            if continent_contouring_area_threshold_steradians > 0:
                # Separate contoured continents above and below the area threshold.
                included_contoured_continents = []
                excluded_contoured_continents = []
                for contoured_continent in contoured_continents:
                    if contoured_continent.get_area() >= continent_contouring_area_threshold_steradians:
                        included_contoured_continents.append(contoured_continent)
                    else:
                        excluded_contoured_continents.append(contoured_continent)
                
                if excluded_contoured_continents:
                    grid_points_inside_excluded_continents = self._find_grid_points_inside_contoured_continents(excluded_contoured_continents)
                    # Exclude grid points inside rejected contoured continents.
                    grid_points_inside_continents[grid_points_inside_excluded_continents] = False
                    # Exclude rejected contoured continents.
                    contoured_continents = included_contoured_continents
        
        # Also include grid points *near* the contoured continents (if a buffer/gap distance was specified).
        if self.continent_contouring_buffer_and_gap_distance_radians_function:
            # These might be outside the contoured continents but within the buffer/gap threshold distance from them.
            grid_points_near_contoured_continents = self._find_grid_points_near_contoured_continents(contoured_continents, age)
            # Only re-generate the contoured continents if there are some grid points near them.
            # It's possible the buffer/gap distance is zero for the current time.
            if np.any(grid_points_near_contoured_continents):
                grid_points_inside_continents[grid_points_near_contoured_continents] = True
                # Contour the grid points that are both inside and near the contoured continents
                # (which in turn contoured the grid points inside the continent polygons).
                contoured_continents = self._find_contoured_continents(grid_points_inside_continents)

        # If an exclusive contoured polygon area threshold was specified.
        if self.continent_exclusion_area_threshold_steradians_function:
            # If the area threshold is non-zero then we now include (instead of exclude) grid points inside any
            # excluding contoured polygons with area below the threshold.
            continent_exclusion_area_threshold_steradians = self.continent_exclusion_area_threshold_steradians_function(age)
            if continent_exclusion_area_threshold_steradians > 0:
                # Find those polygon rings that exclude continental crust and that are below the area threshold.
                # These rings will now include (rather than exclude) continental crust.
                ring_polygons_now_including_continent = []
                for contoured_continent in contoured_continents:
                    # For each polygon that *includes* continental crust we'll use its *interior* rings (which exclude continental crust).
                    for polygon in contoured_continent._polygons_including_continent:
                        for interior_ring_index in range(polygon.get_number_of_interior_rings()):
                            interior_ring_polygon = pygplates.PolygonOnSphere(polygon.get_interior_ring_points(interior_ring_index))
                            if interior_ring_polygon.get_area() < continent_exclusion_area_threshold_steradians:
                                ring_polygons_now_including_continent.append(interior_ring_polygon)
                    # For each polygon that *excludes* continental crust we'll use its *exterior* ring (which excludes continental crust).
                    for polygon in contoured_continent._polygons_excluding_continent:
                        exterior_ring_polygon = pygplates.PolygonOnSphere(polygon.get_exterior_ring_points())
                        if exterior_ring_polygon.get_area() < continent_exclusion_area_threshold_steradians:
                            ring_polygons_now_including_continent.append(exterior_ring_polygon)
                
                if ring_polygons_now_including_continent:
                    # Find the grid points inside the ring polygons that now include (rather than exclude) continental crust.
                    grid_points_inside_ring_polygons_now_including_continent = self._find_grid_points_inside_polygons(ring_polygons_now_including_continent)
                    # Include these grid points.
                    grid_points_inside_continents[grid_points_inside_ring_polygons_now_including_continent] = True
                    # Re-contour the grid points (now that we've included some previously excluded grid points).
                    contoured_continents = self._find_contoured_continents(grid_points_inside_continents)
        
        return contoured_continents


    def _find_grid_points_inside_polygons(
            self,
            polygons):
        """
        Find the latitude/longitude grid points that are inside (one or more of) the specified polygons.

        The grid spacing of these grid points was specified in the constructor.

        Returns a 2D boolean numpy array of shape (num_latitudes, num_longitudes).
        """
        
        # Find the polygon (if any) containing each grid point.
        polygons_containing_points = points_in_polygons.find_polygons_using_points_spatial_tree(
                self.contouring_points,
                self.contouring_points_spatial_tree,
                polygons)

        # Determine which grid points are inside the polygons.
        points_inside_contour = np.full(len(self.contouring_points), False)
        for contouring_point_index in range(len(self.contouring_points)):
            # If the current point is inside any polygon then mark it as such.
            if polygons_containing_points[contouring_point_index] is not None:
                points_inside_contour[contouring_point_index] = True

        # Reshape 1D array as 2D array indexed by (latitude, longitude) - same order as the points.
        return points_inside_contour.reshape((self.contouring_grid_num_latitudes, self.contouring_grid_num_longitudes))


    def _find_grid_points_inside_contoured_continents(
            self,
            contoured_continents):
        """
        Find the latitude/longitude grid points that are inside (one or more of) the specified contoured continents.

        The grid spacing of these grid points was specified in the constructor.

        Returns a 2D boolean numpy array of shape (num_latitudes, num_longitudes).
        """
    
        # Test all grid points against each contoured continent.
        #
        # Note that the original point-in-polygon boolean grid mask (calculated from the input continent polygons before contouring)
        # may not match the result of point-in-contoured-continent tests (done here) since any contoured continents with area below
        # the area threshold would have been removed. So we do our own point-in-contoured-continent tests here.
        points_inside_any_contoured_continent = np.full(len(self.contouring_points), False)
        for contoured_continent in contoured_continents:
            points_inside_contoured_continent = contoured_continent.are_points_inside(self.contouring_points, self.contouring_points_spatial_tree)
            
            # Combine the results of current contoured continent with previous contoured continents.
            #
            # Note that there is typically only a handful of contoured continents in general, so this should not be a bottleneck.
            points_inside_any_contoured_continent[points_inside_contoured_continent] = True

        # Reshape 1D array as 2D array indexed by (latitude, longitude) - same order as the points.
        return points_inside_any_contoured_continent.reshape((self.contouring_grid_num_latitudes, self.contouring_grid_num_longitudes))


    def _find_grid_points_near_contoured_continents(
            self,
            contoured_continents,
            age):
        """
        Find the latitude/longitude grid points that are near (one or more of) the specified contoured continents.

        The grid spacing of these grid points was specified in the constructor.

        The 'age' is only used to look up the time-dependent thresholds (passed into constructor).

        Returns a 2D boolean numpy array of shape (num_latitudes, num_longitudes).
        """

        # Determine which grid points are near the contoured continents.
        points_near_any_contoured_continent = np.full(len(self.contouring_points), False)
        for contoured_continent in contoured_continents:

            # The distance threshold for the current contoured continent.
            distance_threshold_radians = self.continent_contouring_buffer_and_gap_distance_radians_function(age, contoured_continent)
            if distance_threshold_radians > 0:

                # Find the contours (if any) near each point.
                points_near_contoured_continent = proximity_query.find_closest_geometries_to_points_using_points_spatial_tree(
                        self.contouring_points,
                        self.contouring_points_spatial_tree,
                        contoured_continent.get_contours(),
                        distance_threshold_radians=distance_threshold_radians)
                
                for contouring_point_index in range(len(self.contouring_points)):
                    if points_near_contoured_continent[contouring_point_index] is not None:
                        points_near_any_contoured_continent[contouring_point_index] = True

        # Reshape 1D array as 2D array indexed by (latitude, longitude) - same order as the points.
        return points_near_any_contoured_continent.reshape((self.contouring_grid_num_latitudes, self.contouring_grid_num_longitudes))

    
    def _find_contoured_continents(
            self,
            points_inside_contour):
        """
        Find the boundaries of the specified mask of grid points.

        Returns a list of 'ContouredContinent'.
        """
        
        num_latitudes = self.contouring_grid_num_latitudes
        num_longitudes = self.contouring_grid_num_longitudes

        num_latitude_intervals = num_latitudes - 1
        num_longitude_intervals = num_longitudes - 1

        #
        # Use the Marching Squares algorithm.
        #
        # This is a 2D version (surface of the globe) of the 3D Marching Cubes algorithm.
        # However the main difference between this and using skimage.measure.find_contours(),
        # that we used previously and that also uses the Marching Squares algorithm, is we wrap across the
        # dateline and handle the poles. In this way we avoid contour polygons clamped to the dateline.
        #
        # The way we handle wrapping around the dateline is to have grid points on the dateline (ie, at both longitude -180 and 180).
        # This way lat/lon points on the left side of uniform lat/lon grid of points actually map to the same points on the globe
        # as the lat/lon points on the right side of the uniform lat/lon grid of points, and so they will generated the same
        # point-in-continent-polygon and point-near-continent-polygon results. This ensures the Marching Squares algorithm (below)
        # will produce continuous contour segments across the dateline (as we move from a square on one side of the dateline to the
        # adjacent square on the other side).
        #
        # We also handle the poles correctly by having the bottom row of lat/lon points map to the South pole and the top row
        # map to the North pole. Because all points in a (top or bottom) row map to the same point (pole) on the globe they will
        # generate the same point-in-continent-polygon and point-near-continent-polygon results. And because the entire row is either
        # inside or outside a contour the Marching Squares algorithm (below) cannot generate contour segments that penetrate the row.
        # This essentially avoids the problem at the poles.
        #

        #
        # First find those latitude/longitude squares (each square has 4 points from uniform lat/lon grid of points)
        # that have some of its 4 points inside a contour and some outside. These are squares that will contain an edge (segment)
        # of a contour polygon. According to the Marching Squares algorithm there are 16 cases. Two of these have all 4 points
        # either inside or outside (and hence have no segments). Twelve cases have a single segment.
        # And two cases have two segments (because two diagonals points are inside and the other two diagonal points are outside).
        # Here we can choose to either join two separated contour islands or keep them separated.
        # We choose to join them because it makes the algorithm easier - if they were separated then a single square would contain
        # two contours belonging to two separate continents and we'd have to be careful that we visited only the contour belonging
        # to the continent we are currently filling. When we join them then the two contours belong to the same continent.
        #
        # Each segment starts at the middle of one side of the square and ends at the middle of another side.
        # Each side of the square is given an identifier...
        #
        #    ---2---
        #   |       |
        #   1       3
        #   |       |
        #    ---0---
        #
        # ...and each segment records a start and end identifier as a 2-tuple.
        #

        # Records the segments contained by all squares.
        marching_squares = []
        # Records the lat/lon index of only those squares containing one (or two) segments.
        marching_squares_containing_segments = set()
        for latitude_index in range(num_latitude_intervals):
            for longitude_index in range(num_longitude_intervals):

                # See which 4 points of the current square are inside a contour.
                bottom_left_square_inside_contour = points_inside_contour[latitude_index, longitude_index]
                bottom_right_square_inside_contour = points_inside_contour[latitude_index, longitude_index + 1]
                top_left_square_inside_contour = points_inside_contour[latitude_index + 1, longitude_index]
                top_right_square_inside_contour = points_inside_contour[latitude_index + 1, longitude_index + 1]

                # Handle the 16 cases of segments in a square.
                #
                # Store 2 segments (most of the time only 1 is needed).
                # Each segment stores a segment start and end identifier.
                if bottom_left_square_inside_contour:
                    if bottom_right_square_inside_contour:
                        if top_left_square_inside_contour:
                            if top_right_square_inside_contour:
                                segments_in_square = None, None
                            else:
                                segments_in_square = (2,3), None
                        else:
                            if top_right_square_inside_contour:
                                segments_in_square = (1,2), None
                            else:
                                segments_in_square = (1,3), None
                    else:
                        if top_left_square_inside_contour:
                            if top_right_square_inside_contour:
                                segments_in_square = (0,3), None
                            else:
                                segments_in_square = (0,2), None
                        else:
                            if top_right_square_inside_contour:
                                # Choose 2 segments that *do* join two islands.
                                segments_in_square = (0,3), (1,2)
                            else:
                                segments_in_square = (0,1), None
                else:
                    if bottom_right_square_inside_contour:
                        if top_left_square_inside_contour:
                            if top_right_square_inside_contour:
                                segments_in_square = (0,1), None
                            else:
                                # Choose 2 segments that *do* join two islands.
                                segments_in_square = (0,1), (2,3)
                        else:
                            if top_right_square_inside_contour:
                                segments_in_square = (0,2), None
                            else:
                                segments_in_square = (0,3), None
                    else:
                        if top_left_square_inside_contour:
                            if top_right_square_inside_contour:
                                segments_in_square = (1,3), None
                            else:
                                segments_in_square = (1,2), None
                        else:
                            if top_right_square_inside_contour:
                                segments_in_square = (2,3), None
                            else:
                                segments_in_square = None, None
                
                # If current square contains a segment then record that.
                if segments_in_square[0]:
                    marching_squares_containing_segments.add((latitude_index, longitude_index))

                marching_squares.append(segments_in_square)
        
        # Mark those points *inside* any contour as requiring a visit.
        # We also need to remove them once they've been visited.
        points_inside_all_contoured_continents_to_visit = set()
        for latitude_index in range(num_latitudes):
            for longitude_index in range(num_longitudes):
                if points_inside_contour[latitude_index, longitude_index]:
                    points_inside_all_contoured_continents_to_visit.add((latitude_index, longitude_index))
        
        #
        # Generate contoured continents.
        #
        # Each contoured continent is found by picking an arbitrary point inside any contours and expanding around it until we've filled
        # the entire contoured continent. As we expand we detect when we reach a contour that has not yet been generated and generate it
        # for the current contoured continent. This expanding fill can detect more than one contour per contoured continent.
        #
        # This is repeated to find all contoured continents (at which time we will have no more points left to visit inside contours).
        #
        contoured_continents = []
        # Keep visting points *inside* any contour until there are no points left to visit.
        while points_inside_all_contoured_continents_to_visit:
            continent_contours = []

            # Keep a queue of points inside the current ContouredContinent that we will search for contours.
            points_inside_contoured_continent = deque()

            # Get any available point inside any contoured continent.
            lat_lon_indices_of_first_point_inside_contoured_continent = points_inside_all_contoured_continents_to_visit.pop()
            point_index_of_first_point_inside_contoured_continent = (lat_lon_indices_of_first_point_inside_contoured_continent[0] * num_longitudes +
                                                                     lat_lon_indices_of_first_point_inside_contoured_continent[1])
            first_point_inside_contoured_continent = self.contouring_points[point_index_of_first_point_inside_contoured_continent]
            # This will be the first point inside the current ContouredContinent.
            points_inside_contoured_continent.append(lat_lon_indices_of_first_point_inside_contoured_continent)

            # Find the remaining points inside the current ContouredContinent by recursively searching
            # nearbouring points until we reach a contour boundary of the current ContouredContinent.
            while points_inside_contoured_continent:
                # Pop the current point to visit.
                latitude_index, longitude_index = points_inside_contoured_continent.popleft()

                # Search the four squares, adjacent to the current point, for a contour.
                #
                # Note that, for an adjacent square containing a contour, we might already have generated
                # the contour in which case all segments of that contour will have been removed from
                # 'marching_squares_containing_segments' and hence we will be essentially searching for the
                # next contour (if any) of the current contoured continent (eg, an interior hole contour).
                # And note that all contours in the four adjacent squares belong to the current continent because
                # we join continent islands (as opposed to separating them) as described above.
                #
                #  +--+--+
                #  |  |  |
                #  +--o--+
                #  |  |  |
                #  +--+--+
                #

                if latitude_index > 0:
                    if longitude_index > 0:
                        neighbour_square_location = latitude_index - 1, longitude_index - 1
                        if neighbour_square_location in marching_squares_containing_segments:
                            continent_contours.append(
                                self._extract_contour(neighbour_square_location, marching_squares, marching_squares_containing_segments))
                    
                    if longitude_index < num_longitudes - 1:
                        neighbour_square_location = latitude_index - 1, longitude_index
                        if neighbour_square_location in marching_squares_containing_segments:
                            continent_contours.append(
                                self._extract_contour(neighbour_square_location, marching_squares, marching_squares_containing_segments))

                if latitude_index < num_latitude_intervals - 1:
                    if longitude_index > 0:
                        neighbour_square_location = latitude_index, longitude_index - 1
                        if neighbour_square_location in marching_squares_containing_segments:
                            continent_contours.append(
                                self._extract_contour(neighbour_square_location, marching_squares, marching_squares_containing_segments))
                    
                    if longitude_index < num_longitudes - 1:
                        neighbour_square_location = latitude_index, longitude_index
                        if neighbour_square_location in marching_squares_containing_segments:
                            continent_contours.append(
                                self._extract_contour(neighbour_square_location, marching_squares, marching_squares_containing_segments))

                #
                # Propagate outwards from current point to progressively fill the inside of the current contoured continent.
                #
                # This requires visiting up to 8 neighbour points (the '+' in diagram below).
                # Only visit those points that are inside (the contour) and that have not yet been visited.
                #
                #  +--+--+
                #  |  |  |
                #  +--o--+
                #  |  |  |
                #  +--+--+
                #
                # Note that we need to wrap around the dateline (longitude) because we need to visit (and remove) ALL points that
                # are inside the *current* contoured continent (before we move onto the next contoured continent).
                # However we don't need to traverse beyond the poles (latitude) in the same way.
                #

                if latitude_index > 0:
                    neighbour_point_location = latitude_index - 1, longitude_index
                    if neighbour_point_location in points_inside_all_contoured_continents_to_visit:
                        points_inside_contoured_continent.append(neighbour_point_location)
                        points_inside_all_contoured_continents_to_visit.remove(neighbour_point_location)

                    if longitude_index > 0:
                        neighbour_point_location = latitude_index - 1, longitude_index - 1
                        if neighbour_point_location in points_inside_all_contoured_continents_to_visit:
                            points_inside_contoured_continent.append(neighbour_point_location)
                            points_inside_all_contoured_continents_to_visit.remove(neighbour_point_location)
                    else:
                        # Wrap around the dateline.
                        neighbour_point_location = latitude_index - 1, num_longitudes - 1
                        if neighbour_point_location in points_inside_all_contoured_continents_to_visit:
                            points_inside_contoured_continent.append(neighbour_point_location)
                            points_inside_all_contoured_continents_to_visit.remove(neighbour_point_location)
                    
                    if longitude_index < num_longitudes - 1:
                        neighbour_point_location = latitude_index - 1, longitude_index + 1
                        if neighbour_point_location in points_inside_all_contoured_continents_to_visit:
                            points_inside_contoured_continent.append(neighbour_point_location)
                            points_inside_all_contoured_continents_to_visit.remove(neighbour_point_location)
                    else:
                        # Wrap around the dateline.
                        neighbour_point_location = latitude_index - 1, 0
                        if neighbour_point_location in points_inside_all_contoured_continents_to_visit:
                            points_inside_contoured_continent.append(neighbour_point_location)
                            points_inside_all_contoured_continents_to_visit.remove(neighbour_point_location)
                
                if latitude_index < num_latitudes - 1:
                    neighbour_point_location = latitude_index + 1, longitude_index
                    if neighbour_point_location in points_inside_all_contoured_continents_to_visit:
                        points_inside_contoured_continent.append(neighbour_point_location)
                        points_inside_all_contoured_continents_to_visit.remove(neighbour_point_location)

                    if longitude_index > 0:
                        neighbour_point_location = latitude_index + 1, longitude_index - 1
                        if neighbour_point_location in points_inside_all_contoured_continents_to_visit:
                            points_inside_contoured_continent.append(neighbour_point_location)
                            points_inside_all_contoured_continents_to_visit.remove(neighbour_point_location)
                    else:
                        # Wrap around the dateline.
                        neighbour_point_location = latitude_index + 1, num_longitudes - 1
                        if neighbour_point_location in points_inside_all_contoured_continents_to_visit:
                            points_inside_contoured_continent.append(neighbour_point_location)
                            points_inside_all_contoured_continents_to_visit.remove(neighbour_point_location)

                    if longitude_index < num_longitudes - 1:
                        neighbour_point_location = latitude_index + 1, longitude_index + 1
                        if neighbour_point_location in points_inside_all_contoured_continents_to_visit:
                            points_inside_contoured_continent.append(neighbour_point_location)
                            points_inside_all_contoured_continents_to_visit.remove(neighbour_point_location)
                    else:
                        # Wrap around the dateline.
                        neighbour_point_location = latitude_index + 1, 0
                        if neighbour_point_location in points_inside_all_contoured_continents_to_visit:
                            points_inside_contoured_continent.append(neighbour_point_location)
                            points_inside_all_contoured_continents_to_visit.remove(neighbour_point_location)

                if longitude_index > 0:
                    neighbour_point_location = latitude_index, longitude_index - 1
                    if neighbour_point_location in points_inside_all_contoured_continents_to_visit:
                        points_inside_contoured_continent.append(neighbour_point_location)
                        points_inside_all_contoured_continents_to_visit.remove(neighbour_point_location)
                else:
                    # Wrap around the dateline.
                    neighbour_point_location = latitude_index, num_longitudes - 1
                    if neighbour_point_location in points_inside_all_contoured_continents_to_visit:
                        points_inside_contoured_continent.append(neighbour_point_location)
                        points_inside_all_contoured_continents_to_visit.remove(neighbour_point_location)

                if longitude_index < num_longitudes - 1:
                    neighbour_point_location = latitude_index, longitude_index + 1
                    if neighbour_point_location in points_inside_all_contoured_continents_to_visit:
                        points_inside_contoured_continent.append(neighbour_point_location)
                        points_inside_all_contoured_continents_to_visit.remove(neighbour_point_location)
                else:
                    # Wrap around the dateline.
                    neighbour_point_location = latitude_index, 0
                    if neighbour_point_location in points_inside_all_contoured_continents_to_visit:
                        points_inside_contoured_continent.append(neighbour_point_location)
                        points_inside_all_contoured_continents_to_visit.remove(neighbour_point_location)

            # The current contoured continent should have encountered one or more contours since
            # it was filled until it reached a boundary (contour) between continent and ocean.
            if not continent_contours:
                # However it is potentially possible for contintental crust to cover the entire globe (ie, no contours).
                # In this case there must be only one continent, which means no continents so far and all points have been visited.
                if contoured_continents or points_inside_all_contoured_continents_to_visit:
                    raise AssertionError('A single continent covering entire globe must be the only continent')

            # Create a ContouredContinent from the contours.
            # This uses an arbitrary point inside the continent to determine which contours are exterior rings and which are interior rings.
            contoured_continent = self._create_contoured_continent_from_contours(continent_contours, first_point_inside_contoured_continent)
            
            contoured_continents.append(contoured_continent)
        
        return contoured_continents


    def _create_contoured_continent_from_contours(
            self,
            contours,
            any_point_inside_contoured_continent):
        """
        Create a ContouredContinent from the specified contours.
        """

        # For each contour create a ring (a polygon with only an exterior ring).
        contour_rings = [pygplates.PolygonOnSphere(contour) for contour in contours]

        def _insert_contour_into_a_landmass_polygon(contour_ring, landmass_polygons):
            # See if new contour ring contains, or is contained by, any of the existing landmass polygons.
            for polygon_index, polygon in enumerate(landmass_polygons):
                exterior_ring, interior_rings = polygon
                # See if exterior ring *contains* the new contour ring (arbitrarily choose first point on ring).
                if exterior_ring.is_point_in_polygon(contour_ring[0]):
                    # Add the new contour ring as an interior ring.
                    interior_rings.append(contour_ring)
                    return
                # See if the new contour ring *contains* the exterior ring (arbitrarily choose first point on ring).
                elif contour_ring.is_point_in_polygon(exterior_ring[0]):
                    # Current exterior ring shouldn't have any interior rings - the contouring algorithm should guarantee this.
                    if interior_rings:
                        raise AssertionError("Cannot make polygon an interior ring of another polygon if it has an interior ring")
                    # The new contour ring becomes an exterior ring and the current exterior ring becomes an interior ring.
                    new_landmass_polygon = contour_ring, [exterior_ring]  # exterior ring, list of interior rings
                    landmass_polygons[polygon_index] = new_landmass_polygon
                    return
            
            # The new contour ring is not inside (and also does not contain) any of the existing landmass polygons.
            # So create a new landmass polygon.
            new_landmass_polygon = contour_ring, []  # exterior ring, list of interior rings
            landmass_polygons.append(new_landmass_polygon)

        # Arrange the contour rings into a landmass polygon.
        #
        # Note that we can get multiple landmass polygons if there is no landmass that is an exterior ring.
        # In this case the continent covers the entire globe except for a few oceanic holes.
        # However you can't have a single polygon with only interior holes (and no exterior ring).
        # So intead we allow for multiple exterior ring polygons to represent these holes (and treat this as a special case).
        landmass_polygons = []
        for contour_ring in contour_rings:
            _insert_contour_into_a_landmass_polygon(contour_ring, landmass_polygons)
        
        contoured_continent = ContouredContinent()
        
        # Add each landmass polygon to the contoured continent.
        for polygon in landmass_polygons:
            # Convert rings (lists of points) into an actual polygon.
            exterior_ring, interior_rings = polygon
            polygon = pygplates.PolygonOnSphere(exterior_ring, interior_rings)
            # A point inside the continent might actually be outside the current landmass polygon.
            # So determine whether this is the case or not (since a contoured continent needs to know this when it's asked to do point-in-polygon tests).
            polygon_interior_contains_continent = polygon.is_point_in_polygon(any_point_inside_contoured_continent)
            # If the polygon *excludes* continental crust then the contouring algorithm should ensure it doesn't have any interior rings (which should be unreachable).
            if (not polygon_interior_contains_continent) and interior_rings:
                raise AssertionError("Contour polygons excluding continental crust should not have interior rings")
            # Add the polygon, and its contours (rings) as polylines.
            # Note: Converting contours directly from polygons (with no interior rings) ensures the first and last points in each polyline are the same.
            contours_as_polylines = [pygplates.PolylineOnSphere(contour_ring) for contour_ring in contour_rings]
            contoured_continent._add_landmass(polygon, polygon_interior_contains_continent, contours_as_polylines)

        return contoured_continent


    def _extract_contour(
            self,
            first_segment_lat_lon_indices,
            marching_squares,
            marching_squares_containing_segments):
        """
        Follow the segment in marching square at specified lat/lon index around contour back to that first segment.

        Note that a marching square can have two segments, in which case that square represents a thin connecting region between
        two larger islands of the contour (but still all just one contour). That square will get traversed twice (once in one direction
        through one segment and once in another the opposite direction through the second segment of that square).
        """

        contour_points = []

        interval_spacing = self.contouring_point_spacing_degrees
        num_latitude_intervals = self.contouring_grid_num_latitudes - 1
        num_longitude_intervals = self.contouring_grid_num_longitudes - 1

        #
        # When a contour is first encountered (during the caller's expanding fill) a contour ring is generated by starting at the first square
        # found that contains one (or two) segments, which represents the start of that contour. We then pick one of that square's segments
        # (in most cases it'll only have one segment) and generate the first contour point at that segment's start. Note that it doesn't matter
        # which segment we pick (if there's two segments) because the contour ring will traverse back to the second segment (since both segments
        # are part of the same contour because their containing square represents a thin connecting region between two larger islands of the ring).
        # We then find the adjacent square to the segment's end (since a segment ends in the middle of a side of the square we can find the adjacent square).
        # We then find the segment in the adjacent square that starts (or ends) at the that point (the previous segment end). The adjacent
        # square may contain two segments in which case we need to find the correct segment (that continues the previous segment).
        # We generate the next contour point at the segment start and continue this process following the contour through segments of squares
        # until we return to the first segment (thus closing the contour loop).
        #

        #
        # Starting at the first segment, follow the segments in a loop until they return to the first segment (thus forming a contour ring).
        #
        latitude_index, longitude_index = first_segment_lat_lon_indices
        prev_segment_end = None
        while True:

            # Get a segment from the current square.
            segment1, segment2 = marching_squares[latitude_index * num_longitude_intervals + longitude_index]
            # If a square has only one segment then it will be in 'segment1' (not 'segment2').
            if segment1 is None:
                # Shouldn't be able to reach a square that doesn't have any segments (or previously had segments but now has none).
                raise AssertionError('Square has no segments')
            
            # Find a segment in the current square such that the segment start matches the
            # end of the previous segment (in the previous square).
            segment_start, segment_end = segment1
            if prev_segment_end is None:  # first segment of current contour...
                # Mark the start of the contour so that later we know when we've completed the contour.
                first_segment_start = segment_start
            else:
                # Continue from the side of the previous square that contains the end point of the previous segment to the side of
                # the current square that should contain the start point of the current segment (that continues previous segment).
                #
                # The right side of previous square continues to left side of current square.
                # The left side of previous square continues to right side of current square.
                # The top side of previous square continues to bottom side of current square.
                # The bottom side of previous square continues to top side of current square.
                #
                # The adjacency relation means 2->0, 0->2, 3->1 and 1->3...
                #
                #    ---2---
                #   |       |
                #   1       3
                #   |       |
                #    ---0---
                #
                # ...which is satisfied by 2^2->0, 0^2->2, 3^2->1 and 1^2->3 (where '^' is exclusive-or).
                # 
                curr_segment_start = prev_segment_end ^ 0b10

                #
                # Find the right segment (if there's two segments) and reverse the segment if necessary
                # so that previous segment end matches current segment start.
                #
                if curr_segment_start == segment_start:
                    # We're traversing segment in the correct direction.
                    pass
                elif curr_segment_start == segment_end:
                    # Reverse segment direction (swap segment start and end points).
                    segment_start, segment_end = segment_end, segment_start
                else:
                    if segment2:
                        # Segment 1 didn't match so swap it with segment 2 (so it can be used later).
                        segment1, segment2 = segment2, segment1

                        segment_start, segment_end = segment1
                        if curr_segment_start == segment_start:
                            # We're traversing segment in the correct direction.
                            pass
                        elif curr_segment_start == segment_end:
                            # Reverse segment direction (swap segment start and end points).
                            segment_start, segment_end = segment_end, segment_start
                        else:
                            raise AssertionError('Unable to find connecting segment')
                    else:
                        raise AssertionError('Unable to find connecting segment')
            
            # The start position of 'segment'.
            # It will be at the midpoint of a side of the square.
            if segment_start == 0:
                segment_start_latitude = -90.0 + latitude_index * interval_spacing
                segment_start_longitude = -180.0 + (longitude_index + 0.5) * interval_spacing
            elif segment_start == 1:
                segment_start_latitude = -90.0 + (latitude_index + 0.5) * interval_spacing
                segment_start_longitude = -180.0 + longitude_index * interval_spacing
            elif segment_start == 2:
                segment_start_latitude = -90.0 + (latitude_index + 1) * interval_spacing
                segment_start_longitude = -180.0 + (longitude_index + 0.5) * interval_spacing
            else:  # segment_start == 3
                segment_start_latitude = -90.0 + (latitude_index + 0.5) * interval_spacing
                segment_start_longitude = -180.0 + (longitude_index + 1) * interval_spacing

            # Generate a contour point at the start of the current segment.
            contour_point = pygplates.PointOnSphere(segment_start_latitude, segment_start_longitude)
            contour_points.append(contour_point)

            # We've just used 'segment1', so discard it by moving 'segment2' into its position to be used later.
            # And if 'segment2' is None then there are no more segments in current square so discard the entire square.
            marching_squares[latitude_index * num_longitude_intervals + longitude_index] = segment2, None
            if segment2 is None:
                # There are no segments left in the current square, so we're finished with it.
                # Note: This will raise KeyError if not present in 'set'.
                marching_squares_containing_segments.remove((latitude_index, longitude_index))

            # We're moving onto the next segment in the next square.
            prev_segment_end = segment_end

            # Move to the next square connected by the end of the current segment.
            #
            #    ---2---
            #   |       |
            #   1       3
            #   |       |
            #    ---0---
            #
            # As noted above, at each pole there is an entire row of lat/lon grid points that are all either inside or outside a contour.
            # This means the Marching Squares algorithm cannot generate contour segments that penetrate the row. So we should not be able
            # to move beyond the poles.
            #
            # Also as noted above, the both the leftmost and rightmost columns of the lat/lon grid of points will be on the dateline
            # (ie, at both longitude -180 and 180). This means the Marching Squares algorithm will produce continuous contour segments across
            # the dateline (as we move from a square on one side of the dateline to the adjacent square on the other side).
            if prev_segment_end == 0:
                latitude_index -= 1
                if latitude_index < 0:
                    raise AssertionError('Segment entered South Pole')
            elif prev_segment_end == 1:
                longitude_index -= 1
                if longitude_index < 0:
                    # Wrap around the dateline.
                    longitude_index += num_longitude_intervals
            elif prev_segment_end == 2:
                latitude_index += 1
                if latitude_index == num_latitude_intervals:
                    raise AssertionError('Segment entered North Pole')
            else:  # prev_segment_end == 3
                longitude_index += 1
                if longitude_index == num_longitude_intervals:
                    # Wrap around the dateline.
                    longitude_index -= num_longitude_intervals
            
            # See if we're returned to the first square (containing the first segment).
            if first_segment_lat_lon_indices == (latitude_index, longitude_index):
                # And make sure the end of the previous segment matches the start of the first segment.
                # See comment above about adjacency relation for explanatation of exclusive-or.
                if first_segment_start == (prev_segment_end ^ 0b10):
                    # Break out of current contour loop (we've completed the contour).
                    break

        # Return the ring of contour points.
        return contour_points

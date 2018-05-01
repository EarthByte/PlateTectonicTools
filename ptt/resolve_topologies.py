
"""
    Copyright (C) 2015 The University of Sydney, Australia
    
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

from __future__ import print_function
import argparse
import math
import sys
import os.path
import pygplates

DEFAULT_OUTPUT_FILENAME_PREFIX = 'topology_'
DEFAULT_OUTPUT_FILENAME_EXTENSION = 'shp'


def resolve_topologies(rotation_model, topological_features, reconstruction_time, output_filename_prefix, output_filename_extension, anchor_plate_id):
    
    # FIXME: Temporary fix to avoid getting OGR GMT/Shapefile error "Mismatch in field names..." and
    # missing geometries when saving resolved topologies/sections to GMT/Shapefile.
    # It's caused by the OGR writer inside pyglates trying to write out features with different
    # shapefiles attribute field (key) names to the same file. We get around this by removing
    # all shapefile attributes.
    topological_features = pygplates.FeaturesFunctionArgument(topological_features).get_features()
    for topological_feature in topological_features:
        topological_feature.remove(pygplates.PropertyName.gpml_shapefile_attributes)
    
    # Resolve our topological plate polygons (and deforming networks) to the current 'reconstruction_time'.
    # We generate both the resolved topology boundaries and the boundary sections between them.
    resolved_topologies = []
    shared_boundary_sections = []
    pygplates.resolve_topologies(
            topological_features, rotation_model, resolved_topologies, reconstruction_time, shared_boundary_sections, anchor_plate_id)

    # We'll create a feature for each boundary polygon feature and each type of
    # resolved topological section feature we find.
    resolved_topology_features = []
    ridge_transform_boundary_section_features = []
    subduction_boundary_section_features = []
    left_subduction_boundary_section_features = []
    right_subduction_boundary_section_features = []

    # Iterate over the resolved topologies.
    for resolved_topology in resolved_topologies:
        resolved_topology_features.append(resolved_topology.get_resolved_feature())

    # Iterate over the shared boundary sections.
    for shared_boundary_section in shared_boundary_sections:
        
        # Get all the geometries of the current boundary section.
        boundary_section_features = [shared_sub_segment.get_resolved_feature()
                for shared_sub_segment in shared_boundary_section.get_shared_sub_segments()]
        
        # Add the feature to the correct list depending on feature type, etc.
        if shared_boundary_section.get_feature().get_feature_type() == pygplates.FeatureType.create_gpml('SubductionZone'):
            
            # Put all subduction zones in one collection/file.
            subduction_boundary_section_features.extend(boundary_section_features)
            
            # Also put subduction zones in left/right collection/file.
            polarity_property = shared_boundary_section.get_feature().get(
                    pygplates.PropertyName.create_gpml('subductionPolarity'))
            if polarity_property:
                polarity = polarity_property.get_value().get_content()
                if polarity == 'Left':
                    left_subduction_boundary_section_features.extend(boundary_section_features)
                elif polarity == 'Right':
                    right_subduction_boundary_section_features.extend(boundary_section_features)
        else:
            # Put all ridges in one collection/file.
            ridge_transform_boundary_section_features.extend(boundary_section_features)

    if resolved_topology_features:
        # Put the features in a feature collection so we can write them to a file.
        resolved_topology_feature_collection = pygplates.FeatureCollection(resolved_topology_features)
        resolved_topology_features_filename = '{0}boundary_polygons_{1:0.2f}Ma.{2}'.format(
                output_filename_prefix, reconstruction_time, output_filename_extension)
        resolved_topology_feature_collection.write(resolved_topology_features_filename)
        
    if ridge_transform_boundary_section_features:
        # Put the features in a feature collection so we can write them to a file.
        ridge_transform_boundary_section_feature_collection = pygplates.FeatureCollection(ridge_transform_boundary_section_features)
        ridge_transform_boundary_section_features_filename = '{0}ridge_transform_boundaries_{1:0.2f}Ma.{2}'.format(
                output_filename_prefix, reconstruction_time, output_filename_extension)
        ridge_transform_boundary_section_feature_collection.write(ridge_transform_boundary_section_features_filename)
        
    if subduction_boundary_section_features:
        # Put the features in a feature collection so we can write them to a file.
        subduction_boundary_section_feature_collection = pygplates.FeatureCollection(subduction_boundary_section_features)
        subduction_boundary_section_features_filename = '{0}subduction_boundaries_{1:0.2f}Ma.{2}'.format(
                output_filename_prefix, reconstruction_time, output_filename_extension)
        subduction_boundary_section_feature_collection.write(subduction_boundary_section_features_filename)
        
    if left_subduction_boundary_section_features:
        # Put the features in a feature collection so we can write them to a file.
        left_subduction_boundary_section_feature_collection = pygplates.FeatureCollection(left_subduction_boundary_section_features)
        left_subduction_boundary_section_features_filename = '{0}subduction_boundaries_sL_{1:0.2f}Ma.{2}'.format(
                output_filename_prefix, reconstruction_time, output_filename_extension)
        left_subduction_boundary_section_feature_collection.write(left_subduction_boundary_section_features_filename)
        
    if right_subduction_boundary_section_features:
        # Put the features in a feature collection so we can write them to a file.
        right_subduction_boundary_section_feature_collection = pygplates.FeatureCollection(right_subduction_boundary_section_features)
        right_subduction_boundary_section_features_filename = '{0}subduction_boundaries_sR_{1:0.2f}Ma.{2}'.format(
                output_filename_prefix, reconstruction_time, output_filename_extension)
        right_subduction_boundary_section_feature_collection.write(right_subduction_boundary_section_features_filename)


if __name__ == "__main__":

    # Check the imported pygplates version.
    required_version = pygplates.Version(9)
    if not hasattr(pygplates, 'Version') or pygplates.Version.get_imported_version() < required_version:
        print('{0}: Error - imported pygplates version {1} but version {2} or greater is required'.format(
                os.path.basename(__file__), pygplates.Version.get_imported_version(), required_version),
            file=sys.stderr)
        sys.exit(1)


    __description__ = \
    """Resolve topological plate polygons (and deforming networks).

    NOTE: Separate the positional and optional arguments with '--' (workaround for bug in argparse module).
    For example...

    python %(prog)s -r rotations1.rot rotations2.rot -m topologies1.gpml topologies2.gpml -t 10 -- topology_"""

    # The command-line parser.
    parser = argparse.ArgumentParser(description = __description__, formatter_class=argparse.RawDescriptionHelpFormatter)
    
    parser.add_argument('-r', '--rotation_filenames', type=str, nargs='+', required=True,
            metavar='rotation_filename', help='One or more rotation files.')
    parser.add_argument('-m', '--topology_filenames', type=str, nargs='+', required=True,
            metavar='topology_filename', help='One or more topology files.')
    parser.add_argument('-a', '--anchor', type=int, default=0,
            dest='anchor_plate_id',
            help='Anchor plate id used for reconstructing. Defaults to zero.')
    
    parser.add_argument('-t', '--reconstruction_times', type=float, nargs='+', required=True,
            metavar='reconstruction_time',
            help='One or more times at which to reconstruct/resolve topologies.')
    
    parser.add_argument('-e', '--output_filename_extension', type=str,
            default='{0}'.format(DEFAULT_OUTPUT_FILENAME_EXTENSION),
            help="The filename extension of the output files containing the resolved topological boundaries and sections "
                "- the default extension is '{0}' - supported extensions include 'shp', 'gmt' and 'xy'."
                .format(DEFAULT_OUTPUT_FILENAME_EXTENSION))
    
    parser.add_argument('output_filename_prefix', type=str, nargs='?',
            default='{0}'.format(DEFAULT_OUTPUT_FILENAME_PREFIX),
            help="The prefix of the output files containing the resolved topological boundaries and sections "
                "- the default prefix is '{0}'".format(DEFAULT_OUTPUT_FILENAME_PREFIX))
    
    
    # Parse command-line options.
    args = parser.parse_args()
    
    rotation_model = pygplates.RotationModel(args.rotation_filenames)
    
    topological_features = [pygplates.FeatureCollection(topology_filename)
            for topology_filename in args.topology_filenames]
    
    for reconstruction_time in args.reconstruction_times:
        resolve_topologies(
                rotation_model,
                topological_features,
                reconstruction_time,
                args.output_filename_prefix,
                args.output_filename_extension,
                args.anchor_plate_id)

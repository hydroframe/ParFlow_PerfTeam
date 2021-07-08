#!/usr/bin/env python3
####################################################################
## Diff two or more Parflow binary files
## Author Ian J. Bertolacci (ianbertolacci@email.arizona.edu)
## University of Arizona
## Adapted from Parflow's pfb2nc tool
## (github.com/parflow/parflow/blob/master/pftools/prepostproc/pfb2nc.py)
## written by Ketan B. Kulkarni (k.kulkarni@fz-juelich.de) at
## SimLab TerrSys Juelich Supercomputing Centre
## Usage: python pfbdiff.py file1 file2 ...
####################################################################

import sys, itertools, argparse
from struct import *
from pprint import pprint
from collections import OrderedDict

class SelfInterrupt(Exception):
  pass

def read_field( file, bytes, format ):
  raw = file.read(bytes)
  data = list(unpack( format, raw ))[0]
  return data

def read_field_across_files( files, bytes, format ):
  return [ read_field(file, bytes, format) for file in files ]

def is_same( data ):
  data_set = set(data)
  return len(data_set) == 1

def is_different( data ):
  return not is_same( data )

def is_same_within_epsilon( data, epsilon, comparison=lambda a,b, e: abs(a-b) <= e ):
  data_set = sorted(list(set(data)))
  return comparison(data[0], data[-1], epsilon)

def is_different_within_epsilon( data, epsilon, comparison=lambda a,b, e: abs(a-b) <= e ):
  return not is_same_within_epsilon( data=data, epsilon=epsilon, comparison=comparison )

def max_difference( values ):
  sorted_values = sorted(values)
  return abs( sorted_values[0] - sorted_values[-1] )

def min_difference( values ):
  sorted_values = sorted(values)
  differences = [ abs(sorted_values[i] - sorted_values[i+1]) for i in range(len(sorted_values)-1) ]
  return min( differences )

# File metadata fields
file_property_fields = OrderedDict([
  # Start index of the domain
  ("x1", { "name" : "x1", "bytes" : 8, "format" : '>d' } ),
  ("y1", { "name" : "y1", "bytes" : 8, "format" : '>d' } ),
  ("z1", { "name" : "z1", "bytes" : 8, "format" : '>d' } ),
  # Number of points in x, y and z direction
  ("nx", { "name" : "nx", "bytes" : 4, "format" : '>i' } ),
  ("ny", { "name" : "ny", "bytes" : 4, "format" : '>i' } ),
  ("nz", { "name" : "nz", "bytes" : 4, "format" : '>i' } ),
  # dx, dy and dz
  ("dx", { "name" : "dx", "bytes" : 8, "format" : '>d' } ),
  ("dy", { "name" : "dy", "bytes" : 8, "format" : '>d' } ),
  ("dz", { "name" : "dz", "bytes" : 8, "format" : '>d' } ),
  # Number of subdomains
  ("nSubGrid", { "name" : "nSubGrid", "bytes" : 4, "format" : '>i' } ),
])

# Subgrid metadata fields
subgrid_property_fields =  OrderedDict([
  # Subgrid indices and counters
  ("ix", { "name" : "ix", "bytes" : 4, "format" : '>i' }),
  ("iy", { "name" : "iy", "bytes" : 4, "format" : '>i' }),
  ("iz", { "name" : "iz", "bytes" : 4, "format" : '>i' }),
  ("nnx", { "name" : "nnx", "bytes" : 4, "format" : '>i' }),
  ("nny", { "name" : "nny", "bytes" : 4, "format" : '>i' }),
  ("nnz", { "name" : "nnz", "bytes" : 4, "format" : '>i' }),
  ("rx", { "name" : "rx", "bytes" : 4, "format" : '>i' }),
  ("ry", { "name" : "ry", "bytes" : 4, "format" : '>i' }),
  ("rz", { "name" : "rz", "bytes" : 4, "format" : '>i' }),
])

print_levels = {
  "verbose" : 2,
  "normal" : 1,
  "silent" : 0
}


diff_status = {
  "same" : 0,
  "different" : 1,
  "help_message" : 2,
  "argument_error" : 3,
}

def main():

  parser = argparse.ArgumentParser(description="Find difference between two or more pfb files. Returns 0 if no difference between all files.")
  parser.add_argument( 'files', type=str, nargs=argparse.REMAINDER )
  parser.add_argument( '--epsilon', '-e', type=float, default=0.0, help="Minimum absolute difference between to floating point numbers to be considered different." )
  parser.add_argument( '--verbose', '-v', default=False, action='store_const', const=True, help="Print all information." )
  parser.add_argument( '--silent', '-s', default=False, action='store_const', const=True, help="Print no information, including difference messages." )
  parser.add_argument( '--quick', '-q', default=False, action='store_const', const=True, help="Return on first found difference" )
  parser.add_argument( '--status_codes', default=False, action='store_const', const=True, help="Print status codes returned by the program.")

  args = parser.parse_args()

  if args.status_codes:
    for name, code in diff_status.items():
      print(f"{name}: {code}")
    return diff_status["help_message"]

  print_level = print_levels["normal"]

  if args.verbose and args.silent:
    print("Error: cannot use both --silent and --verbose options")
    return diff_status["argument_error"]
  elif args.silent:
    print_level = print_levels["silent"]
  elif args.verbose:
    print_level = print_levels["verbose"]

  fileNames = args.files

  if len(fileNames) < 2:
    print("Must list at least 2 files to compare across")
    return diff_status["argument_error"]

  if print_level >= print_levels["verbose"]:
    print(f"Files: {', '.join(fileNames)}")

  file_handles = [ open(filename, "rb") for filename in fileNames ]

  try:
    detected_difference = False

    # Read all files properties
    file_properties = {
      name : { **field, "values" : read_field_across_files( file_handles, field["bytes"], field["format"] ) }
      for name, field in file_property_fields.items()
    }

    if print_level >= print_levels["verbose"]:
      pprint( file_properties )

    # Compare all fields of all files
    for name, database_property in file_properties.items():
      possibly_different = is_different( database_property["values"] )
      if possibly_different:
        detected_difference = True
        if print_level >= print_levels["normal"]:
          print(f"Difference found for {database_property['name']}:")
          for file_index in range(len(database_property['values'])):
            print(f"{fileNames[file_index]}: {database_property['values'][file_index]}")
        if args.quick:
          raise SelfInterrupt

    nSubGrid = set(file_properties["nSubGrid"]["values"]).pop()

    for gridCounter in range(0, nSubGrid):
      if print_level >= print_levels["verbose"]:
        print(f"Grid #{gridCounter}")

      # Read this subgrid's metadata for all files
      subgrid_properties = {
        name : { **field, "values" : read_field_across_files( file_handles, field["bytes"], field["format"] ) }
        for name, field in subgrid_property_fields.items()
      }

      if print_level >= print_levels["verbose"]:
        pprint( subgrid_properties )

      # Compare all subgrid fields of all files
      for name, subgrid_property in file_properties.items():
        possibly_different = is_different( subgrid_property["values"] )
        if possibly_different:
          detected_difference = True
          if print_level >= print_levels["normal"]:
            print(f"Difference found for {subgrid_property['name']}:")
            for file_index in range(len(subgrid_property['values'])):
              print(f"{fileNames[file_index]}: {subgrid_property['values'][file_index]}")
          if args.quick:
            raise SelfInterrupt

      ix = set(subgrid_properties["ix"]["values"]).pop()
      iy = set(subgrid_properties["iy"]["values"]).pop()
      iz = set(subgrid_properties["iz"]["values"]).pop()
      nnx = set(subgrid_properties["nnx"]["values"]).pop()
      nny = set(subgrid_properties["nny"]["values"]).pop()
      nnz = set(subgrid_properties["nnz"]["values"]).pop()

      # Read all data from subgrids and compare
      # Probably ineffecient, original reads nearly an entire subgrid before comparison
      # TODO Make faster by reading a row or plane
      for (z,y,x) in itertools.product( range(iz, iz+nnz), range(iy, iy+nny), range(ix, ix+nnx) ):
        variable_values = read_field_across_files( file_handles, 8, '>d' )
        if print_level >= print_levels["verbose"]:
          print(f"{(x,y,z)}: {variable_values}")
        possibly_different = is_different_within_epsilon( variable_values, args.epsilon )
        if possibly_different:
          detected_difference = True
          if print_level >= print_levels["normal"]:
            print(f"Difference in subgrid data at {(x,y,z)}")
            for file_index in range(len(variable_values)):
              print(f"{fileNames[file_index]}: {variable_values[file_index]}")
            if len(variable_values) > 2:
              print(f"Maximum Difference: {max_difference(variable_values)}\nMinimum Difference: {min_difference(variable_values)}")
            else:
              print(f"Difference: {max_difference(variable_values)}")
          if args.quick:
            raise SelfInterrupt
  except KeyboardInterrupt:
    pass
  except SelfInterrupt:
    pass
  finally:
    for file_handle in file_handles:
      file_handle.close()

  if detected_difference:
    return diff_status["different"]
  else:
    return diff_status["same"]

if __name__ == "__main__":
  exit( main() )

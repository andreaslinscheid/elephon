#!/bin/bash

parent_path=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )
elephon_bin=$parent_path/../release/elephon

$elephon_bin infile

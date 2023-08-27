#!/bin/bash

# Find the directory this script is in
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

# Add that directory to the PATH
echo "export PATH=\$PATH:$DIR" >> ~/.bashrc

# Source .bashrc to apply changes
source ~/.bashrc

# Start sis
./sis


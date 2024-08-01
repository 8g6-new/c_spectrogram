#!/bin/bash

# Function to list files in a folder
list_files() {
    local folder=$1
    echo "Files in $folder:"
    for file in "$folder"/*; do
        [ -f "$file" ] && echo "  $(basename "$file")"
    done
}

# Check if the first argument is 'run'
run_targets() {
    local folder=$1
    echo "Running all executables in $folder:"
    for file in "$folder"/*; do
        if [ -x "$file" ]; then
            echo "  Running $(basename "$file")"
            ./"$file"
        fi
    done
}

# Check if the first argument is 'run'
if [ "$1" == "run" ]; then
    shift
    if [ "$1" == "all" ]; then
        echo "Running all targets..."
        # List all executables in the Makefiles directory
        run_targets "Makefiles/bins"
    else
        # Run a specific executable
        echo "Running specific target: $1"
        ./"Makefiles/bins/$1"
    fi



# Check if the first argument is 'build'
elif [ "$1" == "build" ]; then
    shift
    if [ "$1" == "all" ]; then
        echo "Building all targets..."
        # List all Makefiles in the Makefiles directory
        list_files "Makefiles"
        # Iterate through each Makefile and build
        for makefile in Makefiles/*; do
            [ -f "$makefile" ] && make -f "$makefile"
        done
    else
        # Build a specific Makefile
        echo "Building specific target: $1"
        make -f "Makefiles/$1"
    fi

# Check if the first argument is 'clean'
elif [ "$1" == "clean" ]; then
    shift
    if [ "$1" == "all" ]; then
        echo "Cleaning all targets..."
        # List all Makefiles in the Makefiles directory
        list_files "Makefiles"
        # Iterate through each Makefile and clean
        for makefile in Makefiles/*; do
            [ -f "$makefile" ] && make -f "$makefile" clean
        done
    else
        # Clean a specific Makefile
        echo "Cleaning specific target: $1"
        make -f "Makefiles/$1" clean
    fi

else
    echo "Unknown command. Use 'run', 'build all', 'build <specific-file>', 'clean all', or 'clean <specific-file>'."
fi

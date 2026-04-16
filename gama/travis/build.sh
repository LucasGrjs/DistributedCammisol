#!/bin/bash

echo "Compiling gama.annotations"
cd ../gama.annotations
mvn clean install "$@"

echo "Compiling gama.processor"
cd ../gama.processor 
mvn clean install "$@"

echo "Compiling gama.parent"
cd ../gama.parent 
mvn clean install "$@"

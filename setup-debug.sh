# Configure package and build debug libraries
meson --prefix=$PWD/debug --buildtype debug build-debug
./build-cc-debug.sh
cd python
meson --prefix=$PWD/debug --buildtype debug build-debug
cd ..
./build-py-debug.sh


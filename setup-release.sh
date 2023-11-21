# Configure package and build release libraries
meson setup "$@" --buildtype release build-release
./build-cc-release.sh
cd python
meson setup "$@" --buildtype release build-release
cd ..
./build-py-release.sh

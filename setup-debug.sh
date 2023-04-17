export PKG_CONFIG_PATH=$PKG_CONFIG_PATH:$CONDA_PREFIX/.local/lib64/pkgconfig
meson setup --prefix=$PWD/debug --buildtype debug build-debug
./build-cc-debug.sh
cd python
meson setup --prefix=$PWD/debug --buildtype debug build-debug
cd ..
./build-py-debug.sh


export PKG_CONFIG_PATH=$PKG_CONFIG_PATH:$CONDA_PREFIX/.local/lib64/pkgconfig
prefix="--prefix=$CONDA_PREFIX/.local"
meson $prefix builddir
sh build-cc.sh
meson $prefix python/builddir
sh build-py.sh


# Compile python debug libraries
meson compile -C python/build-debug && meson install -C python/build-debug && meson test -C python/build-debug

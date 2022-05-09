CXXFLAGS="$CXXFLAGS -O3" meson compile -C python/build && meson install -C python/build && meson test -C python/build

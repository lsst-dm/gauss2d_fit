CXXFLAGS="$CXXFLAGS -O3" meson compile -C python/builddir && meson install -C python/builddir && meson test -C python/builddir


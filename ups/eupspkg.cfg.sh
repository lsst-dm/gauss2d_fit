# EupsPkg config file. Sourced by 'eupspkg'

build_cc()
{
    ./build-cc-release.sh
}

build_py()
{
    ./build-py-release.sh
}

build()
{
    (build_cc && build_py)
}

clean()
{
    ./clean.sh
}

config_cc()
{
    ([ -d "$GAUSS2D_FIT_DIR" ] && ./clean-cc.sh && meson setup --prefix="$GAUSS2D_FIT_DIR/build-release" \
     --buildtype release build-release)
}

config_py()
{
    ([ -d "$GAUSS2D_FIT_DIR" ] && ./clean-py.sh && build_cc \
     && cd python && meson setup --prefix="$GAUSS2D_FIT_DIR/python/build-release" --buildtype release build-release)
}

config()
{
    config_cc && build_cc && config_py
}

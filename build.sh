#!/bin/bash
DIST_DIR="./dist"
mkdir -p "$DIST_DIR"
echo "Building Wasm to: $DIST_DIR"

emcc main.cpp \
    -o "$DIST_DIR/circle-splatting.js" \
    -O3 \
    -msimd128 \
    -s WASM=1 \
    -s ALLOW_MEMORY_GROWTH=1 \
    -s MODULARIZE=1 \
    -s EXPORT_NAME="loadCircleSplattingModule" \
    -s ASSERTIONS=0 \
    -s EXPORT_ES6=1 \
    -s ENVIRONMENT='web' \
    --bind \
    -std=c++20

if [ $? -eq 0 ]; then
    echo "Build sucessful!"
    echo "Files located in: $DIST_DIR"
    if [ -n "$1" ]; then
        DEST_DIR="$1"
        mkdir -p "$DEST_DIR"
        cp "$DIST_DIR/circle-splatting.js" "$DEST_DIR/"
        cp "$DIST_DIR/circle-splatting.wasm" "$DEST_DIR/"
        echo "Copied files to: $DEST_DIR"
        echo "Deployment complete."
    fi
else
    echo "Build failed."
    exit 1
fi

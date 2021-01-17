function build(){
    echo "Building $1.cc..."
    clang++ -o $1 -std=c++2a $1.cc -O3 -ferror-limit=1 -I/usr/local/include/eigen3/ -lglfw -lglew -framework OpenGL -DGL_SILENCE_DEPRECATION
}

if [ "$#" -eq 0 ]; then
    build main && echo "Done Building. Running..." && ./main $@
elif [ "$#" -eq 1 ] && [ "$1" == "test" ]; then
    shift # Remove the first arg (which is "test")
    build test && echo "Done Building. Testing..." && ./test $@
fi

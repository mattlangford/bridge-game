clang++ -o test -std=c++2a test.cc -O3 -ferror-limit=1 -lglfw -lglew -framework OpenGL -DGL_SILENCE_DEPRECATION && ./test $@

clang++ -o main -std=c++2a main.cc -O3 -ferror-limit=1 -lglfw -lglew -framework OpenGL -DGL_SILENCE_DEPRECATION && ./main $@

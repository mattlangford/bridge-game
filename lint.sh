find . -name '*.hh' -o -name '*.cc' | xargs clang-format --style WebKit -i

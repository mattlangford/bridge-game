#include "mesh.hh"

#include <iostream>
#include <string>
#include <unordered_map>

void test_mesh()
{
}

int main(int argc, char* argv[])
{
    std::unordered_map<std::string, std::function<void(void)>> tests;
    tests["mesh"] = &test_mesh;

    if (argc == 2) {
        const std::string name = argv[1];
        auto it = tests.find(name);
        if (it == std::end(tests)) {
            std::cout << "No tests found named '" << name << "'!\n";
            return 1;
        }

        std::cout << "Test: '" << it->first << "'\n";
        try {
            it->second();
        } catch (const std::exception& ex) {
            std::cout << "FAILED: '" << ex.what() << "'\n";
            return 1;
        }
        std::cout << "PASSED\n";
    } else {
        std::cout << "Running all tests!\n";
        for (const auto& [name, func] : tests) {
            std::cout << "Test: '" << name << "'\n";
            try {
                func();
            } catch (const std::exception& ex) {
                std::cout << "FAILED: '" << ex.what() << "'\n";
                continue;
            }
            std::cout << "PASSED\n";
        }
    }

    return 0;
}

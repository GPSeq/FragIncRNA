#include <filesystem>
#include <fstream>
#include <string>
#include <stdexcept>
#include <vector>
#include <iostream>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
namespace fs = std::filesystem;
#include "trusted_archive.hpp"

int main() {
    std::string pattern = "/tmp/frag-trusted-XXXXXX";
    std::vector<char> name(pattern.begin(),pattern.end()); name.push_back('\0');
    if (!mkdtemp(name.data())) return 1;
    fs::path root{name.data()};
    try {
        auto write = [](fs::path const & p,std::string const & value) { std::ofstream out(p,std::ios::binary); out << value; };
        auto reject = [](auto action) {
            try { action(); } catch (std::exception const &) { return; }
            throw std::runtime_error("Expected identity-change rejection");
        };
        for (int scenario = 0; scenario < 4; ++scenario) {
            auto file = root/("archive"+std::to_string(scenario)); write(file,"verified archive bytes");
            TrustedArchive trusted{file}; trusted.check();
            std::ifstream pinned(trusted.pinned_path()); std::string text; std::getline(pinned,text);
            if (text != "verified archive bytes") throw std::runtime_error("Pinned descriptor mismatch");
            if (scenario == 0) { // Same-size edit with restored mtime: ctime still detects it.
                struct stat before{}; stat(file.c_str(),&before);
                write(file,"Verified archive bytes");
                struct timespec stamps[2]{before.st_atim,before.st_mtim}; utimensat(AT_FDCWD,file.c_str(),stamps,0);
            }
            if (scenario == 1) fs::resize_file(file,3);
            if (scenario == 2) { auto replacement = root/"replacement"; write(replacement,"verified archive bytes"); fs::rename(replacement,file); }
            if (scenario == 3) { auto original = root/"original"; fs::rename(file,original); fs::create_symlink(original,file); }
            reject([&] { trusted.check(); });
        }
        auto symlink = root/"symlink"; fs::create_symlink(root/"archive0",symlink);
        reject([&] { TrustedArchive invalid{symlink}; });
        auto fifo = root/"fifo"; mkfifo(fifo.c_str(),0600);
        reject([&] { TrustedArchive invalid{fifo}; });
        fs::remove_all(root); std::cout << "Pinned archive identity checks passed\n";
    } catch (...) { fs::remove_all(root); throw; }
}

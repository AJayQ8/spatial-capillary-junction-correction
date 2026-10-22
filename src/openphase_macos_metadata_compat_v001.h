// Build-only macOS compatibility shim for the pinned OpenPhase source.
// The upstream MetaData.cpp uses Linux's /proc/self/exe in its non-Windows
// path. Keep that source unchanged and map only its readlink("/proc/self/exe")
// calls to the equivalent Darwin executable-path API.
#pragma once

#if defined(__APPLE__)
#include <algorithm>
#include <cstdint>
#include <cstring>
#include <mach-o/dyld.h>
#include <string>
#include <sys/types.h>
#include <unistd.h>

inline ssize_t op_macos_executable_readlink(const char* path, char* buffer,
                                            size_t buffer_size) {
    if (std::strcmp(path, "/proc/self/exe") != 0) {
        return ::readlink(path, buffer, buffer_size);
    }

    uint32_t required = 0;
    if (_NSGetExecutablePath(nullptr, &required) != -1 || required == 0) {
        return -1;
    }
    std::string executable_path(required, '\0');
    if (_NSGetExecutablePath(executable_path.data(), &required) != 0) {
        return -1;
    }

    const size_t path_size = std::strlen(executable_path.c_str());
    const size_t copied = std::min(buffer_size, path_size);
    std::memcpy(buffer, executable_path.data(), copied);
    return static_cast<ssize_t>(copied);
}

#define readlink op_macos_executable_readlink
#endif

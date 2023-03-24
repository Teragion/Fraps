#ifndef LOGGER_H
#define LOGGER_H

#include <cstdio>

class FLogger {
public:
    static void logInfo(const char* s) {
        printf("INFO: %s", s);
    }

    static void logErr(const char* s) {

    }
};

#endif // LOGGER_H
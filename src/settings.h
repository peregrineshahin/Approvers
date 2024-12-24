#ifndef SETTINGS_H
#define SETTINGS_H

#include <stddef.h>
struct settings {
    size_t ttSize;
    bool   clear;
};

extern struct settings settings, delayedSettings;

void process_delayed_settings(void);

#endif

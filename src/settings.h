#ifndef SETTINGS_H
#define SETTINGS_H

struct settings {
  size_t ttSize;
  size_t numThreads;
  bool largePages;
  bool clear;
};

extern struct settings settings, delayedSettings;

void process_delayed_settings(void);

#endif


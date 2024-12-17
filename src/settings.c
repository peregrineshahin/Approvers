#include "search.h"
#include "settings.h"
#include "thread.h"
#include "tt.h"
#include "types.h"

struct settings settings, delayedSettings;

// Process Hash, Threads and LargePages settings.

void process_delayed_settings(void) {
    bool ttChange = delayedSettings.ttSize != settings.ttSize;

    if (settings.numThreads != delayedSettings.numThreads)
    {
        settings.numThreads = delayedSettings.numThreads;
        threads_set_number(settings.numThreads);
    }

    if (ttChange)
    {
        tt_free();
        settings.ttSize     = delayedSettings.ttSize;
        tt_allocate(settings.ttSize);
    }

    if (delayedSettings.clear)
    {
        delayedSettings.clear = false;
        search_clear();
    }
}

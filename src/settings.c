#include "search.h"
#include "settings.h"

struct settings delayedSettings;

void process_delayed_settings(void) {
    if (delayedSettings.clear)
    {
        delayedSettings.clear = false;
        search_clear();
    }
}
